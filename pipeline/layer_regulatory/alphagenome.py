"""Track B — AlphaGenome regulatory divergence scan.

Uses the Google DeepMind AlphaGenome Python SDK to predict regulatory effects
of promoter sequence differences between human and resilient species.

Scientific approach:
    For each Tier1/Tier2 gene × resilient species pair:
    1. Fetch ±8 kb around the TSS from NCBI for both human and ortholog.
    2. Pad both sequences to exactly SEQUENCE_LENGTH_16KB (16,384 bp).
    3. Run AlphaGenome predict_sequence() for CAGE (TSS activity) and DNase
       (chromatin accessibility) tracks using brain tissue ontology, which is
       broadly relevant for cancer resistance biology.
    4. Compare predicted tracks at the ±500 bp window around the TSS center.
    5. Divergence = normalised absolute difference in CAGE + DNase signals.

The 16 KB sequence length is the smallest supported by AlphaGenome and is
well-matched to our 16 kb fetch window, capturing promoter + proximal enhancers
without requiring large genomic downloads per species.

API registration: https://deepmind.google.com/science/alphagenome
(Free for non-commercial research. Set ALPHAGENOME_API_KEY env var.)
"""

import logging
import os
import time
from typing import Optional

import numpy as np
import requests

from db.models import CandidateScore, EvolutionScore, Gene, RegulatoryDivergence, Species
from db.session import get_session
from pipeline.config import get_ncbi_api_key, get_ncbi_email

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# AlphaGenome SDK smallest supported length — 16,384 bp.
# We fetch ±FLANKING_BP around the TSS and pad to this length.
SEQUENCE_LENGTH_16KB = 16_384
FLANKING_BP = SEQUENCE_LENGTH_16KB // 2  # ±8192 bp around TSS

# Window around TSS center where we compare predicted tracks (±500 bp)
TSS_WINDOW_HALF = 500

# Significance threshold: divergence score ≥ this → divergent promoter
EFFECT_THRESHOLD = 0.3

# Tissue ontology for predictions: brain (UBERON:0000955) is broadly expressed
# and captures transcriptional regulation relevant to cancer resistance.
ONTOLOGY_TERMS = ["UBERON:0000955"]

# NCBI eutils
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

# Map species_id → lineage group
LINEAGE_MAP = {
    "naked_mole_rat":      "Rodents",
    "blind_mole_rat":      "Rodents",
    "damaraland_mole_rat": "Rodents",
    "ground_squirrel":     "Rodents",
    "spiny_mouse":         "Rodents",
    "bowhead_whale":       "Cetaceans",
    "right_whale":         "Cetaceans",
    "little_brown_bat":    "Bats",
    "brandts_bat":         "Bats",
    "greenland_shark":     "Sharks",
    "rougheye_rockfish":   "Fishes",
    "amazon_parrot":       "Birds",
    "budgerigar":          "Birds",
    "african_elephant":    "Proboscideans",
    "mouse_lemur":         "Primates",
    "human":               "Primates",
    "axolotl":             "Salamanders",
}


# ---------------------------------------------------------------------------
# API key helpers
# ---------------------------------------------------------------------------

def _get_api_key() -> str:
    key = os.environ.get("ALPHAGENOME_API_KEY", "")
    if not key:
        try:
            from pipeline.config import get_config
            cfg = get_config()
            key = (cfg.get("alphagenome") or {}).get("api_key", "")
            if key and "${" in key:
                key = ""
        except Exception:
            pass
    return key.strip()


# ---------------------------------------------------------------------------
# Symbol helpers
# ---------------------------------------------------------------------------

_SPECIES_SUFFIXES = ("_HUMAN", "_MOUSE", "_RAT", "_DOG", "_CAT", "_WHALE", "_BAT", "_SHARK")


def _hgnc_symbol(full_symbol: str) -> str:
    """Strip UniProt species suffix to get clean HGNC gene symbol.

    E.g. 'AKT1_HUMAN' → 'AKT1', 'AKT1' → 'AKT1'
    """
    for suffix in _SPECIES_SUFFIXES:
        if full_symbol.upper().endswith(suffix):
            return full_symbol[: -len(suffix)]
    return full_symbol


# ---------------------------------------------------------------------------
# NCBI helpers
# ---------------------------------------------------------------------------

def _ncbi_base_params() -> dict:
    return {
        "api_key": get_ncbi_api_key(),
        "email": get_ncbi_email(),
    }


def _get_ncbi_gene_id(gene_symbol: str, taxid: int = 9606) -> Optional[str]:
    """Return NCBI Gene ID for a gene symbol + taxonomy ID.

    Accepts both clean symbols (AKT1) and UniProt entries (AKT1_HUMAN).
    """
    clean_symbol = _hgnc_symbol(gene_symbol)
    try:
        r = requests.get(
            NCBI_ESEARCH,
            params={
                **_ncbi_base_params(),
                "db": "gene",
                "term": f"{clean_symbol}[Gene Name] AND {taxid}[Taxonomy ID]",
                "retmax": 1,
                "retmode": "json",
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        return ids[0] if ids else None
    except Exception as exc:
        log.debug("NCBI esearch %s (taxid=%d): %s", clean_symbol, taxid, exc)
        return None


def _fetch_promoter_sequence_16kb(ncbi_gene_id: str) -> Optional[str]:
    """Fetch ±FLANKING_BP bp around TSS from NCBI, padded to SEQUENCE_LENGTH_16KB.

    Returns a DNA string of exactly SEQUENCE_LENGTH_16KB characters, or None on failure.
    """
    try:
        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={
                **_ncbi_base_params(),
                "db": "gene",
                "id": ncbi_gene_id,
                "retmode": "json",
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        result = r.json().get("result", {}).get(ncbi_gene_id, {})
        genomicinfo = result.get("genomicinfo", [])
        if not genomicinfo:
            return None
        gi = genomicinfo[0]
        chrom_accn = gi.get("chraccver")
        chrstart = gi.get("chrstart")
        if not chrom_accn or chrstart is None:
            return None

        tss = int(chrstart)
        seq_start = max(0, tss - FLANKING_BP)
        seq_end = tss + FLANKING_BP

        time.sleep(0.12)
        r2 = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params={
                **_ncbi_base_params(),
                "db": "nuccore",
                "id": chrom_accn,
                "rettype": "fasta",
                "retmode": "text",
                "seq_start": seq_start,
                "seq_stop": seq_end,
            },
            timeout=60,
        )
        if r2.status_code != 200:
            return None
        lines = r2.text.strip().splitlines()
        seq = "".join(l for l in lines if not l.startswith(">")).upper()
        if not seq:
            return None
        # Pad or trim to exactly SEQUENCE_LENGTH_16KB
        if len(seq) < SEQUENCE_LENGTH_16KB:
            seq = seq.ljust(SEQUENCE_LENGTH_16KB, "N")
        else:
            seq = seq[:SEQUENCE_LENGTH_16KB]
        return seq
    except Exception as exc:
        log.debug("Promoter fetch for gene %s: %s", ncbi_gene_id, exc)
        return None


# ---------------------------------------------------------------------------
# AlphaGenome SDK call
# ---------------------------------------------------------------------------

def _compute_divergence(
    human_seq: str,
    alt_seq: str,
    dna_model,
) -> float:
    """Compare predicted regulatory tracks between human and ortholog sequences.

    Runs AlphaGenome on both sequences and computes normalised divergence
    at the ±TSS_WINDOW_HALF window around the sequence center.

    Returns a divergence score in [0, 1].
    """
    try:
        from alphagenome.models import dna_client as _dna_client

        out_h = dna_model.predict_sequence(
            sequence=human_seq,
            requested_outputs=[
                _dna_client.OutputType.CAGE,
                _dna_client.OutputType.DNASE,
            ],
            ontology_terms=ONTOLOGY_TERMS,
        )
        out_a = dna_model.predict_sequence(
            sequence=alt_seq,
            requested_outputs=[
                _dna_client.OutputType.CAGE,
                _dna_client.OutputType.DNASE,
            ],
            ontology_terms=ONTOLOGY_TERMS,
        )

        center = SEQUENCE_LENGTH_16KB // 2
        w = slice(center - TSS_WINDOW_HALF, center + TSS_WINDOW_HALF)

        # CAGE: average over strands and tracks at TSS window
        h_cage = float(np.mean(out_h.cage.values[w]))
        a_cage = float(np.mean(out_a.cage.values[w]))

        # DNase: chromatin accessibility at TSS window
        h_dnase = float(np.mean(out_h.dnase.values[w]))
        a_dnase = float(np.mean(out_a.dnase.values[w]))

        # Normalised absolute divergence (0 = identical, 1 = maximally different)
        cage_div = abs(h_cage - a_cage) / (max(h_cage, a_cage, 1e-6))
        dnase_div = abs(h_dnase - a_dnase) / (max(h_dnase, a_dnase, 1e-6))

        return round(min((cage_div + dnase_div) / 2.0, 1.0), 4)

    except Exception as exc:
        log.debug("AlphaGenome divergence compute: %s", exc)
        return 0.0


# ---------------------------------------------------------------------------
# Per-gene aggregation
# ---------------------------------------------------------------------------

def _compute_lineage_count(gene_id: str, session) -> int:
    """Count independent lineages with regulatory_score ≥ EFFECT_THRESHOLD."""
    rows = (
        session.query(RegulatoryDivergence)
        .filter(
            RegulatoryDivergence.gene_id == gene_id,
            RegulatoryDivergence.regulatory_score >= EFFECT_THRESHOLD,
        )
        .all()
    )
    return len({LINEAGE_MAP.get(r.species_id) for r in rows if LINEAGE_MAP.get(r.species_id)})


def _gene_regulatory_score(gene_id: str, session) -> float:
    """Gene-level score: max divergence + 0.1 × lineage_count, capped at 1.0."""
    rows = session.query(RegulatoryDivergence).filter_by(gene_id=gene_id).all()
    if not rows:
        return 0.0
    effects = [r.regulatory_score for r in rows if r.regulatory_score is not None]
    if not effects:
        return 0.0
    max_effect = max(effects)
    lineage_count = _compute_lineage_count(gene_id, session)
    return round(min(max_effect + 0.1 * lineage_count, 1.0), 4)


# ---------------------------------------------------------------------------
# Main pipeline function
# ---------------------------------------------------------------------------

def run_alphagenome_track() -> int:
    """Run Track B: regulatory divergence for Tier1 + Tier2 genes.

    Scope: Tier1 + Tier2 genes only (from Phase 1 EvolutionScore tiers).
    These are the ~60 genes we care about — no point running AlphaGenome on
    all 818 expression candidates when only ~60 pass Phase 1.

    Returns number of RegulatoryDivergence rows written.
    """
    api_key = _get_api_key()
    if not api_key:
        log.warning(
            "ALPHAGENOME_API_KEY not set. "
            "Set the env var or add alphagenome.api_key to config/environment.yml. "
            "Skipping Track B regulatory divergence."
        )
        return 0

    # Import SDK here so the container fails fast if alphagenome is not installed
    try:
        from alphagenome.models import dna_client as _dna_client
    except ImportError:
        log.error(
            "alphagenome package not installed. "
            "Add 'alphagenome' to the Docker image requirements and rebuild."
        )
        return 0

    # Create the model client once (avoids re-authenticating per gene)
    log.info("Initialising AlphaGenome client (SEQUENCE_LENGTH_16KB)...")
    dna_model = _dna_client.create(api_key)

    with get_session() as session:
        # Only process Tier1 + Tier2 genes — tier lives on CandidateScore, not EvolutionScore
        tier_gene_ids = {
            r.gene_id
            for r in session.query(CandidateScore.gene_id)
            .filter(CandidateScore.tier.in_(["Tier1", "Tier2"]))
            .all()
        }
        all_species = session.query(Species).filter(Species.id != "human").all()
        species_list = [s for s in all_species if s.id in LINEAGE_MAP]
        genes = session.query(Gene).filter(Gene.id.in_(tier_gene_ids)).all()

    log.info(
        "AlphaGenome Track B: %d Tier1/2 genes × %d species = %d pairs to score.",
        len(genes),
        len(species_list),
        len(genes) * len(species_list),
    )

    # Pre-fetch human promoter sequences (cache by gene_id)
    human_seqs: dict[str, str] = {}
    for gene in genes:
        ncbi_id = gene.human_gene_id if (gene.human_gene_id and gene.human_gene_id.isdigit()) else None
        if not ncbi_id:
            ncbi_id = _get_ncbi_gene_id(gene.gene_symbol, taxid=9606)
            time.sleep(0.12)
        if not ncbi_id:
            log.debug("No NCBI gene ID for %s — skipping.", gene.gene_symbol)
            continue
        seq = _fetch_promoter_sequence_16kb(ncbi_id)
        if seq:
            human_seqs[gene.id] = seq
        time.sleep(0.12)

    log.info("Fetched human promoter sequences for %d / %d genes.", len(human_seqs), len(genes))

    # Score each gene × species pair
    written = 0
    for gene in genes:
        human_seq = human_seqs.get(gene.id)
        if not human_seq:
            continue

        for species in species_list:
            ortholog_ncbi_id = _get_ncbi_gene_id(gene.gene_symbol, taxid=species.taxid)
            if not ortholog_ncbi_id:
                continue
            time.sleep(0.12)

            alt_seq = _fetch_promoter_sequence_16kb(ortholog_ncbi_id)
            if not alt_seq:
                continue
            time.sleep(0.12)

            divergence = _compute_divergence(human_seq, alt_seq, dna_model)

            with get_session() as session:
                row = (
                    session.query(RegulatoryDivergence)
                    .filter_by(gene_id=gene.id, species_id=species.id)
                    .first()
                )
                if row is None:
                    row = RegulatoryDivergence(gene_id=gene.id, species_id=species.id)
                    session.add(row)
                row.promoter_divergence = divergence
                row.regulatory_score = divergence
                written += 1

        log.info(
            "AlphaGenome: %s done (%d rows so far).",
            gene.gene_symbol,
            written,
        )

    # Second pass: compute lineage counts
    with get_session() as session:
        gene_ids_done = {r.gene_id for r in session.query(RegulatoryDivergence.gene_id).distinct()}
        for gid in gene_ids_done:
            lc = _compute_lineage_count(gid, session)
            gene_score = _gene_regulatory_score(gid, session)
            for row in session.query(RegulatoryDivergence).filter_by(gene_id=gid).all():
                row.lineage_count = lc
                row.regulatory_score = gene_score  # persist per-gene aggregate score
        log.debug("Lineage count pass complete for %d genes.", len(gene_ids_done))

    log.info("AlphaGenome Track B complete: %d regulatory divergence rows written.", written)
    return written
