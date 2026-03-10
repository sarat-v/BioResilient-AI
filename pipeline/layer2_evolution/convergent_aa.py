"""True convergent amino acid substitution detection.

Current limitation of the sliding-window approach:
  Sequence divergence ≠ convergence.
  A motif that differs in 5 lineages might have 5 *different* substitutions
  — that's parallel divergence, not convergence.

This module detects *true convergence*: the **same** (or biochemically
equivalent) amino acid change evolving independently at the same position
in multiple lineages.

Method:
  For each DivergentMotif, across all species that show divergence at the
  same gene position:
    - Extract the exact amino acid substitution (human_ref → animal_alt) at
      each divergent position.
    - Count how many independent *lineages* carry the same substitution
      (identical or Miyata-equivalent alt amino acid) at that position.
    - Store the result as convergent_aa_count on DivergentMotif.

A high convergent_aa_count means: the *same biochemical change* evolved
independently — strong evidence of adaptive convergence rather than drift.

This supplements (not replaces) the existing lineage-count convergence metric.
"""

import logging
from collections import defaultdict
from typing import Optional

from db.models import DivergentMotif, Gene, Ortholog, Species
from db.session import get_session

log = logging.getLogger(__name__)

# Miyata biochemical equivalence groups
_MIYATA_GROUPS: dict[str, int] = {
    "G": 0, "A": 0, "V": 0, "L": 0, "I": 0,   # nonpolar aliphatic
    "F": 1, "W": 1, "Y": 1,                      # aromatic
    "S": 2, "T": 2, "C": 2,                      # polar uncharged small
    "N": 3, "Q": 3, "H": 3,                      # polar uncharged large
    "D": 4, "E": 4,                               # negatively charged
    "K": 5, "R": 5,                               # positively charged
    "P": 6, "M": 6,                               # special
}


def _miyata_group(aa: str) -> Optional[int]:
    return _MIYATA_GROUPS.get(aa.upper())


def _same_biochemical_change(alt1: str, alt2: str) -> bool:
    """Return True if two amino acid substitution destinations are equivalent."""
    if alt1 == alt2:
        return True
    g1 = _miyata_group(alt1)
    g2 = _miyata_group(alt2)
    return g1 is not None and g1 == g2


def compute_convergent_aa_count(
    gene_id: str,
    species_lineage_map: dict[str, str],
) -> dict[str, int]:
    """Count independently-derived convergent amino acid substitutions per motif.

    For a given gene, loads all DivergentMotif rows across all orthologs, then
    for each alignment position:
      1. Collects the (human_aa, animal_aa) substitution per species.
      2. Groups species by lineage.
      3. Counts how many lineages carry the same alt amino acid (or Miyata-
         equivalent) at that position.

    Returns {motif_id: convergent_aa_count}.
    """
    with get_session() as session:
        motifs = (
            session.query(DivergentMotif, Ortholog.species_id)
            .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
            .filter(Ortholog.gene_id == gene_id)
            .all()
        )

    if not motifs:
        return {}

    # Build per-position substitution map:
    # {abs_position: [(lineage, human_aa, animal_aa), ...]}
    position_subs: dict[int, list[tuple[str, str, str]]] = defaultdict(list)

    motif_obj_by_id: dict[str, DivergentMotif] = {}

    for motif, species_id in motifs:
        motif_obj_by_id[motif.id] = motif
        lineage = species_lineage_map.get(species_id, "Other")
        if lineage in ("Primates", "Other") or species_id == "human":
            continue

        # Walk through the motif window and extract per-position substitutions
        for i, (h_aa, a_aa) in enumerate(zip(motif.human_seq, motif.animal_seq)):
            if h_aa == "-" or a_aa == "-" or h_aa == a_aa:
                continue
            abs_pos = motif.start_pos + i
            position_subs[abs_pos].append((lineage, h_aa, a_aa))

    # For each motif, count how many lineages converge at each of its positions
    motif_conv_counts: dict[str, int] = {}

    for motif_id, motif in motif_obj_by_id.items():
        max_lineage_conv = 0

        for i in range(len(motif.human_seq)):
            abs_pos = motif.start_pos + i
            subs = position_subs.get(abs_pos, [])
            if len(subs) < 2:
                continue

            # Group by alt amino acid (Miyata-equivalent)
            # Count how many *distinct* lineages share the same biochemical change
            # Use a representative alt_aa per group and track lineages
            lineage_seen: dict[str, str] = {}  # lineage → alt_aa
            for lineage, h_aa, a_aa in subs:
                if lineage not in lineage_seen:
                    lineage_seen[lineage] = a_aa

            # Find the most common biochemical change across lineages
            # Group lineages by whether their alt_aa is Miyata-equivalent
            best_count = 0
            lineage_alts = list(lineage_seen.values())
            for j, aa_ref in enumerate(lineage_alts):
                count = sum(
                    1 for aa_other in lineage_alts
                    if _same_biochemical_change(aa_ref, aa_other)
                )
                best_count = max(best_count, count)

            max_lineage_conv = max(max_lineage_conv, best_count)

        motif_conv_counts[motif_id] = max_lineage_conv

    return motif_conv_counts


def annotate_convergent_aa(gene_ids: Optional[list[str]] = None) -> int:
    """Annotate all DivergentMotif rows with convergent_aa_count.

    Args:
        gene_ids: Optional list of gene IDs to process (default: all genes).

    Returns:
        Number of motifs with convergent_aa_count >= 2 (true convergence signal).
    """
    with get_session() as session:
        if gene_ids is None:
            gene_ids = [g.id for g in session.query(Gene).all()]

        # Build species → lineage map from DB
        species_lineage_map = {
            s.id: (s.lineage_group or "Other")
            for s in session.query(Species).all()
        }

    log.info(
        "Computing convergent amino acid substitutions for %d genes...",
        len(gene_ids),
    )
    truly_convergent = 0

    with get_session() as session:
        for gene_id in gene_ids:
            counts = compute_convergent_aa_count(gene_id, species_lineage_map)
            for motif_id, count in counts.items():
                motif = session.get(DivergentMotif, motif_id)
                if motif:
                    motif.convergent_aa_count = count
                    if count >= 2:
                        truly_convergent += 1
        session.commit()

    log.info(
        "Convergent AA annotation complete: %d motifs with true convergence (≥2 lineages, same substitution).",
        truly_convergent,
    )
    return truly_convergent


def run_convergent_aa_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step7b."""
    return annotate_convergent_aa(gene_ids)
