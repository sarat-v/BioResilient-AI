"""True convergent amino acid substitution detection.

Detects *true convergence*: the **same** (or biochemically equivalent) amino
acid change evolving independently at the same position in multiple lineages.

Method:
  For each DivergentMotif, across all species that show divergence at the
  same gene position:
    - Extract the exact amino acid substitution (human_ref -> animal_alt) at
      each divergent position.
    - Count how many independent *lineages* carry the same substitution
      (identical or Miyata-equivalent alt amino acid) at that position.
    - Store the result as convergent_aa_count on DivergentMotif.
"""

import logging
from collections import Counter, defaultdict
from typing import Optional

from db.models import Gene, Species
from db.session import get_session

log = logging.getLogger(__name__)

# Miyata biochemical equivalence groups — map each AA to its group ID
_MIYATA: dict[str, int] = {
    "G": 0, "A": 0, "V": 0, "L": 0, "I": 0,
    "F": 1, "W": 1, "Y": 1,
    "S": 2, "T": 2, "C": 2,
    "N": 3, "Q": 3, "H": 3,
    "D": 4, "E": 4,
    "K": 5, "R": 5,
    "P": 6, "M": 6,
}


def annotate_convergent_aa(gene_ids: Optional[list[str]] = None) -> int:
    """Annotate all DivergentMotif rows with convergent_aa_count.

    Optimised path: raw SQL load, numpy-free in-memory computation using
    pre-computed Miyata groups and Counter-based convergence detection,
    then batched bulk UPDATE via unnest.
    """
    from sqlalchemy import text

    with get_session() as session:
        if gene_ids is None:
            gene_ids = [g.id for g in session.query(Gene).all()]

        species_lineage_map: dict[str, str] = {
            s.id: (s.lineage_group or "Other")
            for s in session.query(Species).all()
        }

        log.info("Loading motif rows for %d genes via raw SQL...", len(gene_ids))
        result = session.execute(
            text("""
                SELECT dm.id, dm.start_pos, dm.human_seq, dm.animal_seq,
                       o.gene_id, o.species_id
                FROM divergent_motif dm
                JOIN ortholog o ON o.id = dm.ortholog_id
                WHERE o.gene_id = ANY(:gene_ids)
            """),
            {"gene_ids": gene_ids},
        )
        rows = result.fetchall()

    log.info("Loaded %d motif rows. Grouping by gene...", len(rows))

    # Group by gene_id — store lightweight tuples, not ORM objects
    # Each entry: (motif_id, start_pos, human_seq, animal_seq, species_id)
    gene_motifs: dict[str, list[tuple]] = defaultdict(list)
    for motif_id, start_pos, human_seq, animal_seq, gene_id, species_id in rows:
        gene_motifs[gene_id].append((str(motif_id), start_pos, human_seq, animal_seq, species_id))

    del rows

    log.info("Computing convergent AA counts for %d genes...", len(gene_motifs))

    skip_lineages = {"Primates", "Other"}
    motif_counts: dict[str, int] = {}
    genes_done = 0

    for gene_id, motifs in gene_motifs.items():
        # Phase 1: build position -> {lineage: miyata_group} map
        # position_lineage_groups[abs_pos][lineage] = miyata_group_of_alt_aa
        position_lineage_groups: dict[int, dict[str, int]] = defaultdict(dict)
        motif_ranges: list[tuple[str, int, int]] = []

        for motif_id, start_pos, human_seq, animal_seq, species_id in motifs:
            motif_ranges.append((motif_id, start_pos, min(len(human_seq), len(animal_seq))))

            lineage = species_lineage_map.get(species_id, "Other")
            if lineage in skip_lineages or species_id == "human":
                continue

            seq_len_cmp = min(len(human_seq), len(animal_seq))
            for i in range(seq_len_cmp):
                h_aa = human_seq[i]
                a_aa = animal_seq[i]
                if h_aa == "-" or a_aa == "-" or h_aa == a_aa:
                    continue
                abs_pos = start_pos + i
                if lineage not in position_lineage_groups[abs_pos]:
                    grp = _MIYATA.get(a_aa.upper(), -1)
                    position_lineage_groups[abs_pos][lineage] = grp

        # Phase 2: for each position, find max lineages sharing same Miyata group
        # Pre-compute per-position max convergence count
        pos_max_conv: dict[int, int] = {}
        for abs_pos, lineage_groups in position_lineage_groups.items():
            if len(lineage_groups) < 2:
                continue
            group_counts = Counter(lineage_groups.values())
            group_counts.pop(-1, None)
            if group_counts:
                pos_max_conv[abs_pos] = max(group_counts.values())

        # Phase 3: for each motif, find the max convergence across its positions
        for motif_id, start_pos, seq_len in motif_ranges:
            max_conv = 0
            for i in range(seq_len):
                conv = pos_max_conv.get(start_pos + i, 0)
                if conv > max_conv:
                    max_conv = conv
            motif_counts[motif_id] = max_conv

        genes_done += 1
        if genes_done % 2000 == 0:
            log.info("  genes processed: %d/%d, motifs so far: %d",
                     genes_done, len(gene_motifs), len(motif_counts))

    log.info("Computation complete: %d motifs scored. Writing to DB...", len(motif_counts))

    # Bulk write using batched UPDATE with unnest
    truly_convergent = sum(1 for c in motif_counts.values() if c >= 2)
    items = list(motif_counts.items())
    BATCH = 10000
    updated = 0

    for i in range(0, len(items), BATCH):
        chunk = items[i : i + BATCH]
        ids = [mid for mid, _ in chunk]
        counts = [c for _, c in chunk]
        with get_session() as session:
            session.execute(text("SET LOCAL statement_timeout = '0'"))
            session.execute(
                text("""
                    UPDATE divergent_motif
                    SET convergent_aa_count = data.cnt
                    FROM (SELECT unnest(:ids) AS id, unnest(:counts) AS cnt) AS data
                    WHERE divergent_motif.id::text = data.id
                """),
                {"ids": ids, "counts": counts},
            )
            session.commit()
        updated += len(chunk)
        if updated % 100000 == 0 or updated == len(items):
            log.info("  DB write progress: %d/%d motifs", updated, len(items))

    log.info(
        "Convergent AA annotation complete: %d motifs with true convergence (>=2 lineages, same substitution).",
        truly_convergent,
    )
    return truly_convergent


def run_convergent_aa_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step7b."""
    return annotate_convergent_aa(gene_ids)
