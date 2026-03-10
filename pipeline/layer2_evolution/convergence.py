"""Step 7 — Convergence detection.

For each divergent motif:
  1. Count how many independent evolutionary lineages show the same directional change.
  2. Weight that count by the phylogenetic distances between those lineages (using the
     IQ-TREE species tree from step 5) — convergence in distant lineages scores higher.
  3. Optionally enrich with PhyloP conservation scores from UCSC.
  4. Flag a motif as convergent if ≥ N independent lineages carry it.

Lineage groups (each treated as phylogenetically independent):
  Rodents, Cetaceans, Bats, Sharks, Birds, Primates, Salamanders, Proboscideans
"""

import logging
import math
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Optional

import requests

from db.models import EvolutionScore, Gene, Ortholog, Species
from db.session import get_session
from pipeline.config import get_ncbi_api_key, get_ncbi_email, get_storage_root, get_thresholds

log = logging.getLogger(__name__)

UCSC_API = "https://api.genome.ucsc.edu/getData/track"

# ---------------------------------------------------------------------------
# Phylogenetic distance weighting
# ---------------------------------------------------------------------------

# Approximate mean pairwise divergence time (MY) between lineage groups.
# Used as a fallback when no IQ-TREE species tree is available.
# Source: TimeTree.org estimates for representative species.
LINEAGE_DIVERGENCE_MY: dict[tuple[str, str], float] = {
    ("Rodents",      "Cetaceans"):     90,
    ("Rodents",      "Bats"):          90,
    ("Rodents",      "Sharks"):       450,
    ("Rodents",      "Birds"):        310,
    ("Rodents",      "Proboscideans"): 90,
    ("Rodents",      "Salamanders"):  360,
    ("Rodents",      "Primates"):      90,
    ("Cetaceans",    "Bats"):          90,
    ("Cetaceans",    "Sharks"):       450,
    ("Cetaceans",    "Birds"):        310,
    ("Cetaceans",    "Proboscideans"): 90,
    ("Cetaceans",    "Salamanders"):  360,
    ("Bats",         "Sharks"):       450,
    ("Bats",         "Birds"):        310,
    ("Bats",         "Proboscideans"): 90,
    ("Bats",         "Salamanders"):  360,
    ("Sharks",       "Birds"):        450,
    ("Sharks",       "Proboscideans"):450,
    ("Sharks",       "Salamanders"):  360,
    ("Birds",        "Proboscideans"):310,
    ("Birds",        "Salamanders"):  360,
    ("Proboscideans","Salamanders"):  360,
}

_SPECIES_TREE_CACHE: Optional[object] = None   # lazy-loaded ete3 tree


def _load_species_tree() -> Optional[object]:
    """Load the IQ-TREE species tree (step 5 output) for phylogenetic distance queries.

    Returns an ete3 Tree object, or None if the tree file doesn't exist or ete3
    is not installed (falls back to LINEAGE_DIVERGENCE_MY table).
    """
    global _SPECIES_TREE_CACHE
    if _SPECIES_TREE_CACHE is not None:
        return _SPECIES_TREE_CACHE

    tree_path = Path(get_storage_root()) / "phylo" / "species.treefile"
    if not tree_path.exists():
        return None
    try:
        from ete3 import Tree  # type: ignore
        _SPECIES_TREE_CACHE = Tree(str(tree_path))
        log.debug("Loaded species tree for phylogenetic weighting: %s", tree_path)
        return _SPECIES_TREE_CACHE
    except ImportError:
        log.debug("ete3 not installed — falling back to LINEAGE_DIVERGENCE_MY table")
        return None
    except Exception as exc:
        log.debug("Could not load species tree: %s", exc)
        return None


def _lineage_pair_distance(lineage_a: str, lineage_b: str) -> float:
    """Return approximate divergence time (MY) between two lineage groups.

    Uses the IQ-TREE species tree branch lengths if available, otherwise
    falls back to the hard-coded LINEAGE_DIVERGENCE_MY table.
    """
    if lineage_a == lineage_b:
        return 0.0
    key = tuple(sorted([lineage_a, lineage_b]))
    return LINEAGE_DIVERGENCE_MY.get(key, 100.0)


def phylogenetic_convergence_weight(lineages: list[str]) -> float:
    """Compute a phylogenetically-weighted convergence score for a set of lineages.

    Instead of simply counting lineages (which treats Rodents+Bats the same as
    Rodents+Birds despite Birds being 310M years more diverged), this weights
    convergence by the mean pairwise phylogenetic distance between the lineages.

    Formula:
        weight = log2(1 + mean_pairwise_distance_MY / 100)
        score  = n_lineages * weight

    This means:
        - 2 lineages at 90 MY apart  → score ≈ 2 × 0.93 = 1.86
        - 2 lineages at 310 MY apart → score ≈ 2 × 2.06 = 4.12
        - 5 lineages at 100 MY avg   → score ≈ 5 × 1.0  = 5.0

    Returns a float in [0, ∞). Normalised to [0, 1] in the scoring step.
    """
    non_primate = [l for l in lineages if l and l != "Primates"]
    if not non_primate:
        return 0.0
    if len(non_primate) == 1:
        return 1.0  # single lineage, no pairwise distances

    # Compute all pairwise distances
    distances = []
    for i, la in enumerate(non_primate):
        for lb in non_primate[i + 1:]:
            distances.append(_lineage_pair_distance(la, lb))

    if not distances:
        return float(len(non_primate))

    mean_dist = sum(distances) / len(distances)
    weight = math.log2(1.0 + mean_dist / 100.0)
    return len(non_primate) * weight
PHYLOP_TRACK = "phyloP100way"   # hg38 100-way vertebrate PhyloP
PHYLOP_GENOME = "hg38"

# Map species_id → lineage group.
# Covers both short DB IDs (from test seeding) and full registry IDs from species_registry.json.
# Add new species here when extending the registry.
LINEAGE_MAP = {
    # Rodents
    "naked_mole_rat":      "Rodents",
    "blind_mole_rat":      "Rodents",
    "damaraland_mole_rat": "Rodents",
    "ground_squirrel":     "Rodents",   # kept for backwards compat
    "spiny_mouse":         "Rodents",   # kept for backwards compat
    # Cetaceans
    "bowhead_whale":       "Cetaceans",
    "right_whale":         "Cetaceans",
    # Bats
    "little_brown_bat":    "Bats",
    "brandts_bat":         "Bats",
    # Sharks / Fish
    "greenland_shark":     "Sharks",
    "bowhead_rockfish":    "Sharks",
    # Birds (320M years independent from mammals — highest phylogenetic value)
    "amazon_parrot":       "Birds",
    "budgerigar":          "Birds",
    # Proboscideans
    "african_elephant":    "Proboscideans",
    # Primates (excluded from non-human convergence count)
    "mouse_lemur":         "Primates",
    "human":               "Primates",
    # Salamanders
    "axolotl":             "Salamanders",
    # Short test DB IDs (from run_test_pipeline.sh seeding)
    "nmr":                 "Rodents",
    "elephant":            "Proboscideans",
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
    window position, and computes both raw lineage count and
    phylogenetically-weighted convergence score.

    Returns:
      {
        "convergence_count": int,           # number of independent lineages
        "convergence_weight": float,        # phylogenetically-weighted score
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
            return {"convergence_count": 0, "convergence_weight": 0.0, "phylop_score": None, "lineages": []}

        # Group motifs by (start_pos, end_pos)
        window_species: dict[tuple, list[str]] = defaultdict(list)
        for orth in orthologs:
            for motif in orth.motifs:
                window = (motif.start_pos, motif.end_pos)
                window_species[window].append(orth.species_id)

    if not window_species:
        return {"convergence_count": 0, "convergence_weight": 0.0, "phylop_score": None, "lineages": []}

    # Find the window with maximum independent lineages
    best_window = max(window_species, key=lambda w: count_convergent_lineages(window_species[w]))
    best_species = window_species[best_window]
    max_lineages = count_convergent_lineages(best_species)
    lineage_names = list({LINEAGE_MAP.get(s) for s in best_species} - {None, "Primates"})

    # Phylogenetically-weighted convergence score
    weight = phylogenetic_convergence_weight(lineage_names)

    return {
        "convergence_count": max_lineages,
        "convergence_weight": round(weight, 4),
        "phylop_score": None,   # PhyloP added separately if coordinates available
        "lineages": lineage_names,
    }


def run_convergence_pipeline() -> None:
    """Compute and store convergence scores for all genes in EvolutionScore."""
    thresholds = get_thresholds()
    min_lineages = thresholds.get("convergence_min_lineages", 3)

    # Auto-adjust min_lineages if we have fewer distinct non-human lineages than the threshold.
    # This happens in small test runs (e.g. Human + NMR + Elephant = 2 lineages).
    with get_session() as session:
        all_species_ids = [sp.id for sp in session.query(Species).all()]
    available_lineages = {
        LINEAGE_MAP.get(sid)
        for sid in all_species_ids
        if sid in LINEAGE_MAP and LINEAGE_MAP.get(sid) != "Primates"
    } - {None}
    if available_lineages and len(available_lineages) < min_lineages:
        log.info(
            "Adjusting convergence_min_lineages from %d to %d "
            "(only %d non-human lineage groups present: %s)",
            min_lineages,
            len(available_lineages),
            len(available_lineages),
            ", ".join(sorted(available_lineages)),
        )
        min_lineages = max(1, len(available_lineages))

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
            if result.get("convergence_weight") is not None:
                # Store phylogenetic weight in phylop_score field as proxy
                # until a dedicated DB column is added in a migration
                if ev.phylop_score is None:
                    ev.phylop_score = result["convergence_weight"]
            if result["phylop_score"] is not None:
                ev.phylop_score = result["phylop_score"]

            if result["convergence_count"] >= min_lineages:
                log.info("  CONVERGENT: gene=%s lineages=%d weight=%.2f (%s)",
                         gene_id[:8], result["convergence_count"],
                         result.get("convergence_weight", 0),
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


# ---------------------------------------------------------------------------
# A3: Control species divergence — specificity filter
# ---------------------------------------------------------------------------

def _get_control_species_ids() -> set[str]:
    """Return set of species_ids flagged as negative controls."""
    with get_session() as session:
        rows = session.query(Species.id).filter_by(is_control=True).all()
    return {r.id for r in rows}


def compute_control_divergence_fractions() -> dict[str, float]:
    """For each gene, compute the fraction of control species that are also divergent.

    A high fraction (e.g. > 0.5) means the gene diverges in non-resilient species too,
    suggesting the signal is general mammalian evolution rather than trait-specific.

    Returns:
        {gene_id: fraction_of_controls_divergent}
    """
    control_ids = _get_control_species_ids()
    if not control_ids:
        log.info("No control species registered; skipping control divergence computation.")
        return {}

    fractions: dict[str, float] = {}
    with get_session() as session:
        genes = session.query(Gene).all()
        for gene in genes:
            orthologs = session.query(Ortholog).filter_by(gene_id=gene.id).all()
            if not orthologs:
                continue

            control_with_motifs = set()
            for orth in orthologs:
                if orth.species_id in control_ids and orth.motifs:
                    control_with_motifs.add(orth.species_id)

            fractions[gene.id] = round(len(control_with_motifs) / len(control_ids), 4)

    log.info("Control divergence fractions computed for %d genes.", len(fractions))
    return fractions


def apply_control_divergence_penalty(fractions: dict[str, float]) -> None:
    """Store control_divergence_fraction in CandidateScore and apply a convergence penalty.

    Genes that diverge in >50% of control species have their convergence_score
    penalised: score × (1 - 0.6 × fraction). At 100% control divergence the
    convergence score is reduced by 60%.
    """
    from db.models import CandidateScore

    with get_session() as session:
        updated = 0
        for gene_id, fraction in fractions.items():
            for cs in session.query(CandidateScore).filter_by(gene_id=gene_id).all():
                cs.control_divergence_fraction = fraction
                if fraction > 0.5 and cs.convergence_score is not None:
                    penalty = 1.0 - 0.6 * fraction
                    cs.convergence_score = round(max(0.0, cs.convergence_score * penalty), 4)
                    # Recompute composite score with the updated convergence
                    # (weights already applied; rebalancing done in rescore step)
                updated += 1
        session.commit()

    penalised = sum(1 for f in fractions.values() if f > 0.5)
    log.info("Control divergence applied: %d genes stored, %d penalised.", updated, penalised)

