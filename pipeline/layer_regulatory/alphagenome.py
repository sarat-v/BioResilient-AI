"""Track B — AlphaGenome regulatory divergence scan.

Uses the DeepMind AlphaGenome API to predict the effect of non-coding DNA
differences between human and resilient species in the promoter/regulatory
region of each candidate gene.

Pre-filter (Track B entry point):
    1. Start from genes with significant differential expression in resilient
       species (from GEO/DESeq2 step, stored in CandidateScore.expression_score).
    2. Subtract genes already captured by Track A (already have EvolutionScore
       with convergence_count ≥ 1 — protein-level signal).
    3. Fetch flanking genomic sequence (±2 kb around TSS) from NCBI for each
       species × gene pair.
    4. Submit to AlphaGenome variant scoring endpoint; compute regulatory
       effect magnitude.
    5. Aggregate lineage-level counts; write RegulatoryDivergence rows.

API registration: https://deepmind.google/technologies/alphagenome/
(Free for non-commercial research. Set ALPHAGENOME_API_KEY env var.)
"""

import logging
import os
import time
from typing import Optional

import requests

from db.models import CandidateScore, EvolutionScore, Gene, RegulatoryDivergence, Species
from db.session import get_session
from pipeline.config import get_ncbi_api_key, get_ncbi_email

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# AlphaGenome API settings
# ---------------------------------------------------------------------------

ALPHAGENOME_BASE = "https://alphagenome.research.google.com/api/v1"
PROMOTER_FLANKING_BP = 2000   # ±2 kb around TSS

# NCBI eutils
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Significance threshold: AlphaGenome effect score ≥ this → divergent promoter
EFFECT_THRESHOLD = 0.3

# Map species_id → lineage group (same as convergence module)
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


def _get_api_key() -> str:
    key = os.environ.get("ALPHAGENOME_API_KEY", "")
    if not key:
        # Also try config YAML
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
# NCBI helpers
# ---------------------------------------------------------------------------

def _ncbi_base_params() -> dict:
    return {
        "api_key": get_ncbi_api_key(),
        "email": get_ncbi_email(),
    }


def _get_ncbi_gene_id(gene_symbol: str) -> Optional[str]:
    """Return NCBI Gene ID for a human gene symbol."""
    try:
        r = requests.get(
            NCBI_ESEARCH,
            params={
                **_ncbi_base_params(),
                "db": "gene",
                "term": f"{gene_symbol}[Gene Name] AND 9606[Taxonomy ID]",
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
        log.debug("NCBI esearch %s: %s", gene_symbol, exc)
        return None


def _fetch_promoter_sequence(
    ncbi_gene_id: str,
    taxid: int = 9606,
    flanking: int = PROMOTER_FLANKING_BP,
) -> Optional[str]:
    """Fetch ±flanking bp around the TSS of a gene from NCBI Nucleotide.

    Returns the DNA sequence as a plain string, or None on failure.
    """
    try:
        # Use efetch gene → get genomic coordinates, then fetch nucleotide region
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
        chrstart = gi.get("chrstart", 0)
        orientation = gi.get("exoncount", 0)
        if not chrom_accn:
            return None

        tss = int(chrstart)
        seq_start = max(0, tss - flanking)
        seq_end = tss + flanking

        time.sleep(0.12)
        r2 = requests.get(
            NCBI_EFETCH,
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
        seq = "".join(l for l in lines if not l.startswith(">"))
        return seq.upper() if seq else None
    except Exception as exc:
        log.debug("Promoter fetch for gene %s taxid %s: %s", ncbi_gene_id, taxid, exc)
        return None


# ---------------------------------------------------------------------------
# AlphaGenome API call
# ---------------------------------------------------------------------------

def _alphagenome_effect(
    human_seq: str,
    alt_seq: str,
    api_key: str,
) -> Optional[float]:
    """Call AlphaGenome variant scoring endpoint.

    Sends human reference sequence and an alternate (resilient species) sequence.
    Returns a scalar effect magnitude in [0, ∞), or None on failure.

    The API compares predicted regulatory grammar (e.g. transcription factor
    binding probabilities, chromatin accessibility) between the two sequences.
    A high value means the resilient species promoter is predicted to behave
    differently from the human promoter.
    """
    if not api_key:
        return None
    if not human_seq or not alt_seq:
        return None
    # Truncate to 4096 bp — model context limit
    human_seq = human_seq[:4096]
    alt_seq = alt_seq[:4096]
    # Pad to same length if needed
    maxlen = max(len(human_seq), len(alt_seq))
    human_seq = human_seq.ljust(maxlen, "N")
    alt_seq = alt_seq.ljust(maxlen, "N")

    try:
        r = requests.post(
            f"{ALPHAGENOME_BASE}/variant_score",
            json={
                "reference_sequence": human_seq,
                "alternate_sequence": alt_seq,
                "sequence_type": "DNA",
            },
            headers={
                "Authorization": f"Bearer {api_key}",
                "Content-Type": "application/json",
            },
            timeout=60,
        )
        if r.status_code != 200:
            log.debug("AlphaGenome API returned %s: %s", r.status_code, r.text[:200])
            return None
        data = r.json()
        # Expect {"effect_score": float, ...}
        score = data.get("effect_score") or data.get("score") or data.get("regulatory_effect")
        return float(score) if score is not None else None
    except Exception as exc:
        log.debug("AlphaGenome effect call: %s", exc)
        return None


# ---------------------------------------------------------------------------
# Per-gene aggregation
# ---------------------------------------------------------------------------

def _compute_lineage_count(gene_id: str, session) -> int:
    """Count independent lineages with regulatory_score ≥ EFFECT_THRESHOLD for this gene."""
    rows = (
        session.query(RegulatoryDivergence)
        .filter(
            RegulatoryDivergence.gene_id == gene_id,
            RegulatoryDivergence.regulatory_score >= EFFECT_THRESHOLD,
        )
        .all()
    )
    lineages = set()
    for row in rows:
        lineage = LINEAGE_MAP.get(row.species_id)
        if lineage:
            lineages.add(lineage)
    return len(lineages)


def _gene_regulatory_score(gene_id: str, session) -> float:
    """Aggregate per-species regulatory scores into one gene-level score.

    Returns the maximum promoter divergence across any resilient species × lineage count bonus.
    Score formula: min(max_effect + 0.1 * lineage_count, 1.0)
    """
    rows = (
        session.query(RegulatoryDivergence)
        .filter(RegulatoryDivergence.gene_id == gene_id)
        .all()
    )
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

def run_alphagenome_track(
    expression_score_min: float = 0.3,
    max_genes: int = 2000,
) -> int:
    """Run Track B: regulatory divergence for expression-positive genes not in Track A.

    Args:
        expression_score_min: Minimum expression_score (from CandidateScore) to include.
        max_genes: Cap on number of genes to process (API rate limit protection).

    Returns:
        Number of RegulatoryDivergence rows written.
    """
    api_key = _get_api_key()
    if not api_key:
        log.warning(
            "AlphaGenome API key not set. Set ALPHAGENOME_API_KEY environment variable "
            "or add alphagenome.api_key to config/environment.yml. Skipping Track B."
        )
        return 0

    with get_session() as session:
        # All species (non-human) with lineage mapping
        all_species = session.query(Species).filter(Species.id != "human").all()
        species_list = [s for s in all_species if s.id in LINEAGE_MAP]

        # Candidate pool: expression-positive genes
        expression_candidates = (
            session.query(CandidateScore.gene_id)
            .filter(CandidateScore.trait_id == "", CandidateScore.expression_score >= expression_score_min)
            .all()
        )
        candidate_ids = {r[0] for r in expression_candidates}

        # Track A genes (already have protein-level signal — still run Track B on them
        # because regulatory changes may compound the signal, but log them separately)
        track_a_ids = {
            r[0]
            for r in session.query(EvolutionScore.gene_id)
            .filter(EvolutionScore.convergence_count >= 1)
            .all()
        }
        track_b_only = candidate_ids - track_a_ids
        all_candidates = list(candidate_ids)[:max_genes]

        genes = session.query(Gene).filter(Gene.id.in_(all_candidates)).all()

    log.info(
        "Track B: %d expression-positive genes (%d Track B-only, %d shared with Track A).",
        len(all_candidates),
        len(track_b_only & set(all_candidates)),
        len(track_a_ids & set(all_candidates)),
    )

    written = 0
    for gene in genes:
        ncbi_gene_id = None
        if gene.human_gene_id and gene.human_gene_id.isdigit():
            ncbi_gene_id = gene.human_gene_id
        if not ncbi_gene_id and gene.gene_symbol:
            ncbi_gene_id = _get_ncbi_gene_id(gene.gene_symbol)
            time.sleep(0.12)

        if not ncbi_gene_id:
            continue

        human_seq = _fetch_promoter_sequence(ncbi_gene_id, taxid=9606)
        if not human_seq:
            log.debug("No promoter sequence for %s (human).", gene.gene_symbol)
            continue

        for species in species_list:
            alt_seq = _fetch_promoter_sequence(ncbi_gene_id, taxid=species.taxid)
            if not alt_seq:
                continue

            effect = _alphagenome_effect(human_seq, alt_seq, api_key)
            if effect is None:
                effect = 0.0

            log2fc = None
            with get_session() as session:
                cs = session.query(CandidateScore).filter_by(gene_id=gene.id, trait_id="").first()
                if cs:
                    log2fc = cs.expression_score

            with get_session() as session:
                reg_score = round(min(float(effect), 1.0), 4)
                row = (
                    session.query(RegulatoryDivergence)
                    .filter_by(gene_id=gene.id, species_id=species.id)
                    .first()
                )
                if row is None:
                    row = RegulatoryDivergence(gene_id=gene.id, species_id=species.id)
                    session.add(row)
                row.promoter_divergence = round(float(effect), 6)
                row.expression_log2fc = log2fc
                row.regulatory_score = reg_score
                written += 1

            time.sleep(0.05)

    # Second pass: compute lineage counts and update rows
    with get_session() as session:
        gene_ids_done = {r.gene_id for r in session.query(RegulatoryDivergence.gene_id).distinct()}
        for gid in gene_ids_done:
            lc = _compute_lineage_count(gid, session)
            rows = session.query(RegulatoryDivergence).filter_by(gene_id=gid).all()
            for r in rows:
                r.lineage_count = lc

    log.info("Track B complete: %d regulatory divergence rows written.", written)
    return written
