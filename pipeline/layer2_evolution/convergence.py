"""Step 7 — Convergence detection.

For each divergent motif:
  1. Query the UCSC REST API for PhyloP conservation scores at each position.
  2. Count how many independent evolutionary lineages show the same directional change.
  3. Flag a motif as convergent if ≥ 3 independent lineages carry it.

Lineage groups (each treated as phylogenetically independent):
  Rodents, Cetaceans, Bats, Sharks, Primates, Salamanders
"""

import logging
import re
import time
from collections import defaultdict
from typing import Optional

import requests

from db.models import EvolutionScore, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_ncbi_api_key, get_ncbi_email, get_thresholds

log = logging.getLogger(__name__)

UCSC_API = "https://api.genome.ucsc.edu/getData/track"
PHYLOP_TRACK = "phyloP100way"   # hg38 100-way vertebrate PhyloP
PHYLOP_GENOME = "hg38"

# Map species_id → lineage group
LINEAGE_MAP = {
    "naked_mole_rat":     "Rodents",
    "damaraland_mole_rat": "Rodents",
    "ground_squirrel":    "Rodents",
    "spiny_mouse":        "Rodents",
    "bowhead_whale":      "Cetaceans",
    "little_brown_bat":   "Bats",
    "greenland_shark":    "Sharks",
    "bowhead_rockfish":   "Sharks",
    "african_elephant":   "Proboscideans",
    "mouse_lemur":        "Primates",
    "axolotl":            "Salamanders",
    "human":              "Primates",
}


def get_phylop_score(chrom: str, start: int, end: int) -> Optional[float]:
    """Query UCSC REST API for mean PhyloP score in a genomic window.

    Note: In Phase 1 we use the human genome coordinates as proxy. A proper
    implementation would lift-over motif positions to hg38. Here we return a
    mean score for a representative locus if coordinates are unavailable.

    Returns None if the API is unreachable or the region has no data.
    """
    try:
        params = {
            "genome": PHYLOP_GENOME,
            "track": PHYLOP_TRACK,
            "chrom": chrom,
            "start": start,
            "end": end,
        }
        r = requests.get(UCSC_API, params=params, timeout=15)
        if r.status_code != 200:
            return None
        data = r.json()
        values = data.get(PHYLOP_TRACK, [])
        if not values:
            return None
        # Average the scores across positions
        scores = [v.get("value", 0) for v in values if "value" in v]
        return round(sum(scores) / len(scores), 4) if scores else None
    except Exception as exc:
        log.debug("PhyloP API error: %s", exc)
        return None


def count_convergent_lineages(motif_species: list[str]) -> int:
    """Count how many independent lineage groups are represented.

    Args:
        motif_species: List of species_ids that carry this motif.

    Returns:
        Number of distinct lineage groups.
    """
    lineages = {LINEAGE_MAP.get(sid) for sid in motif_species if sid in LINEAGE_MAP}
    lineages.discard(None)
    lineages.discard("Primates")   # Human is Primates — exclude from count
    return len(lineages)


def compute_convergence_for_gene(gene_id: str) -> dict:
    """Compute convergence score for one gene.

    Loads all DivergentMotifs for this gene's orthologs, groups by
    window position, and counts independent lineages per window.

    Returns:
      {
        "convergence_count": int,
        "phylop_score": float | None,
        "lineages": [str, ...],
      }
    """
    with get_session() as session:
        orthologs = (
            session.query(Ortholog)
            .filter_by(gene_id=gene_id)
            .all()
        )

        if not orthologs:
            return {"convergence_count": 0, "phylop_score": None, "lineages": []}

        # Group motifs by (start_pos, end_pos)
        window_species: dict[tuple, list[str]] = defaultdict(list)
        for orth in orthologs:
            for motif in orth.motifs:
                window = (motif.start_pos, motif.end_pos)
                window_species[window].append(orth.species_id)

    if not window_species:
        return {"convergence_count": 0, "phylop_score": None, "lineages": []}

    # Find the window with maximum independent lineages
    best_window = max(window_species, key=lambda w: count_convergent_lineages(window_species[w]))
    best_species = window_species[best_window]
    max_lineages = count_convergent_lineages(best_species)
    lineage_names = list({LINEAGE_MAP.get(s) for s in best_species} - {None, "Primates"})

    return {
        "convergence_count": max_lineages,
        "phylop_score": None,   # PhyloP added separately if coordinates available
        "lineages": lineage_names,
    }


def run_convergence_pipeline() -> None:
    """Compute and store convergence scores for all genes in EvolutionScore."""
    thresholds = get_thresholds()
    min_lineages = thresholds.get("convergence_min_lineages", 3)

    with get_session() as session:
        gene_ids = [ev.gene_id for ev in session.query(EvolutionScore).all()]

    if not gene_ids:
        # Fall back to all genes with motifs
        with get_session() as session:
            gene_ids = list({
                o.gene_id
                for o in session.query(Ortholog)
                .join(Ortholog.motifs)
                .all()
            })

    log.info("Computing convergence for %d genes...", len(gene_ids))

    with get_session() as session:
        for gene_id in gene_ids:
            result = compute_convergence_for_gene(gene_id)

            ev = session.get(EvolutionScore, gene_id)
            if ev is None:
                ev = EvolutionScore(gene_id=gene_id)
                session.add(ev)

            ev.convergence_count = result["convergence_count"]
            if result["phylop_score"] is not None:
                ev.phylop_score = result["phylop_score"]

            if result["convergence_count"] >= min_lineages:
                log.info("  CONVERGENT: gene=%s lineages=%d (%s)",
                         gene_id[:8], result["convergence_count"],
                         ", ".join(result["lineages"]))

    log.info("Convergence computation complete.")


NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def build_chrom_map(gene_ids: Optional[list[str]] = None) -> dict[str, tuple[str, int]]:
    """Build {gene_id: (chrom, start)} for hg38 using NCBI Gene API.

    Args:
        gene_ids: Optional list of our Gene.id (UUID). If None, uses all genes with EvolutionScore.

    Returns:
        chrom_map suitable for enrich_phylop_scores(). Chromosome as str (e.g. '1'), start as 0-based.
    """
    with get_session() as session:
        if gene_ids is not None:
            genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()
        else:
            ev_genes = session.query(EvolutionScore.gene_id).distinct().all()
            ids = [r[0] for r in ev_genes]
            if not ids:
                ids = [r[0] for r in session.query(Gene.id).all()]
            genes = session.query(Gene).filter(Gene.id.in_(ids)).all()

    chrom_map: dict[str, tuple[str, int]] = {}
    params_base = {
        "api_key": get_ncbi_api_key(),
        "email": get_ncbi_email(),
    }

    for gene in genes:
        ncbi_gene_id = None
        if gene.human_gene_id and re.match(r"^\d+$", str(gene.human_gene_id)):
            ncbi_gene_id = gene.human_gene_id
        if not ncbi_gene_id and gene.gene_symbol:
            try:
                r = requests.get(
                    NCBI_ESEARCH,
                    params={
                        **params_base,
                        "db": "gene",
                        "term": f"{gene.gene_symbol}[Gene Name] AND 9606[Taxonomy ID]",
                        "retmax": 1,
                        "retmode": "json",
                    },
                    timeout=15,
                )
                if r.status_code == 200:
                    data = r.json()
                    id_list = data.get("esearchresult", {}).get("idlist", [])
                    if id_list:
                        ncbi_gene_id = id_list[0]
            except Exception as exc:
                log.debug("NCBI esearch for %s: %s", gene.gene_symbol, exc)
            time.sleep(0.12)

        if not ncbi_gene_id:
            continue
        try:
            r = requests.get(
                NCBI_ESUMMARY,
                params={
                    **params_base,
                    "db": "gene",
                    "id": ncbi_gene_id,
                    "retmode": "json",
                },
                timeout=15,
            )
            if r.status_code != 200:
                continue
            data = r.json()
            result = data.get("result", {}).get(ncbi_gene_id, {})
            chrom = result.get("chromosome")
            start = result.get("genomicinfo", [{}])[0].get("chrstart") if result.get("genomicinfo") else None
            if chrom is not None and start is not None:
                chrom_map[gene.id] = (str(chrom), int(start))
        except Exception as exc:
            log.debug("NCBI esummary for gene %s: %s", ncbi_gene_id, exc)
        time.sleep(0.12)

    log.info("Built chrom_map for %d genes (hg38).", len(chrom_map))
    return chrom_map


def enrich_phylop_scores(chrom_map: dict[str, tuple[str, int]]) -> None:
    """Optionally enrich EvolutionScore with PhyloP scores.

    Args:
        chrom_map: {gene_id: (chrom, start_pos)} — hg38 coordinates.
                   Typically populated from an NCBI gene annotation lookup.
    """
    with get_session() as session:
        for gene_id, (chrom, start) in chrom_map.items():
            score = get_phylop_score(chrom, start, start + 100)
            if score is None:
                continue
            ev = session.get(EvolutionScore, gene_id)
            if ev:
                ev.phylop_score = score
            time.sleep(0.1)   # Respect UCSC rate limits

    log.info("PhyloP scores enriched for %d genes.", len(chrom_map))
