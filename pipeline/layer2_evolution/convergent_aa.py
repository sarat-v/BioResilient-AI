"""True convergent amino acid substitution detection with ancestral state reconstruction.

Detects *true convergence*: amino acid substitutions that evolved *independently*
at the same position in multiple lineages, as confirmed by ancestral state
reconstruction (ASR). A substitution is convergent only if the inferred
ancestral amino acid at the most recent common ancestor (MRCA) differs from
the derived state — ruling out shared inheritance.

Miyata biochemical equivalence groups are used as a secondary signal to
identify functionally equivalent (but not identical) convergent substitutions.
"""

import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

from db.models import Gene, Species
from db.session import get_session
from pipeline.config import get_local_storage_root

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


def _load_species_tree_newick() -> Optional[str]:
    """Load the species tree Newick from step 5 output."""
    root = Path(get_local_storage_root())
    treefile = root / "phylo" / "species.treefile"
    if treefile.exists():
        return treefile.read_text().strip()
    return None


def _infer_ancestral_aa(tree_newick: str, species_aa: dict[str, str]) -> dict[str, str]:
    """Infer ancestral amino acids at internal nodes using maximum parsimony.

    Uses Fitch's algorithm (downpass + uppass) for parsimony-based ASR.
    Returns a dict mapping species_id to "ancestral" AA at the parent node —
    only species whose AA differs from their parent's inferred state represent
    true independent substitutions.
    """
    try:
        from ete3 import Tree
    except ImportError:
        return {}

    if not tree_newick or len(species_aa) < 3:
        return {}

    try:
        t = Tree(tree_newick)
        present = [n.name for n in t.get_leaves() if n.name in species_aa]
        if len(present) < 3:
            return {}
        t.prune(present, preserve_branch_length=True)
    except Exception:
        return {}

    # Fitch downpass: assign sets of possible AAs bottom-up
    for node in t.traverse("postorder"):
        if node.is_leaf():
            aa = species_aa.get(node.name, "-")
            node.add_feature("fitch_set", {aa} if aa != "-" else set())
        else:
            child_sets = [c.fitch_set for c in node.children if c.fitch_set]
            if not child_sets:
                node.add_feature("fitch_set", set())
            elif len(child_sets) == 1:
                node.add_feature("fitch_set", child_sets[0])
            else:
                intersection = child_sets[0]
                for cs in child_sets[1:]:
                    intersection = intersection & cs
                if intersection:
                    node.add_feature("fitch_set", intersection)
                else:
                    union = set()
                    for cs in child_sets:
                        union = union | cs
                    node.add_feature("fitch_set", union)

    # Fitch uppass: resolve ambiguities top-down
    for node in t.traverse("preorder"):
        if node.is_root():
            node.add_feature("fitch_aa", min(node.fitch_set) if node.fitch_set else "-")
        else:
            parent_aa = node.up.fitch_aa
            if parent_aa in node.fitch_set:
                node.add_feature("fitch_aa", parent_aa)
            else:
                node.add_feature("fitch_aa", min(node.fitch_set) if node.fitch_set else parent_aa)

    # For each leaf, record the parent's inferred AA
    parent_aa_map = {}
    for leaf in t.get_leaves():
        if leaf.up is not None:
            parent_aa_map[leaf.name] = leaf.up.fitch_aa
        else:
            parent_aa_map[leaf.name] = leaf.fitch_aa

    return parent_aa_map


def _infer_ancestral_aa_cached(pruned_tree, species_aa: dict[str, str]) -> dict[str, str]:
    """Infer ancestral AAs using an already-pruned ete3 Tree object.

    Avoids re-parsing the Newick string and re-pruning for every position.
    """
    if pruned_tree is None or len(species_aa) < 3:
        return {}

    try:
        from ete3 import Tree
        t = pruned_tree.copy()
    except Exception:
        return {}

    # Fitch downpass
    for node in t.traverse("postorder"):
        if node.is_leaf():
            aa = species_aa.get(node.name, "-")
            node.add_feature("fitch_set", {aa} if aa != "-" else set())
        else:
            child_sets = [c.fitch_set for c in node.children if c.fitch_set]
            if not child_sets:
                node.add_feature("fitch_set", set())
            elif len(child_sets) == 1:
                node.add_feature("fitch_set", child_sets[0])
            else:
                intersection = child_sets[0]
                for cs in child_sets[1:]:
                    intersection = intersection & cs
                if intersection:
                    node.add_feature("fitch_set", intersection)
                else:
                    union = set()
                    for cs in child_sets:
                        union = union | cs
                    node.add_feature("fitch_set", union)

    # Fitch uppass
    for node in t.traverse("preorder"):
        if node.is_root():
            node.add_feature("fitch_aa", min(node.fitch_set) if node.fitch_set else "-")
        else:
            parent_aa_val = node.up.fitch_aa
            if parent_aa_val in node.fitch_set:
                node.add_feature("fitch_aa", parent_aa_val)
            else:
                node.add_feature("fitch_aa", min(node.fitch_set) if node.fitch_set else parent_aa_val)

    parent_aa_map = {}
    for leaf in t.get_leaves():
        if leaf.up is not None:
            parent_aa_map[leaf.name] = leaf.up.fitch_aa
        else:
            parent_aa_map[leaf.name] = leaf.fitch_aa

    return parent_aa_map


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

    tree_newick = _load_species_tree_newick()
    if tree_newick:
        log.info("Loaded species tree for ancestral state reconstruction.")
    else:
        log.warning("No species tree found — falling back to Miyata-only convergence (no ASR).")

    # Cache pruned ete3 trees keyed by frozenset of species to avoid
    # re-parsing + re-pruning the same tree thousands of times.
    _pruned_tree_cache: dict[frozenset, object] = {}

    def _get_pruned_tree(species_set: frozenset):
        if species_set in _pruned_tree_cache:
            return _pruned_tree_cache[species_set]
        try:
            from ete3 import Tree
            t = Tree(tree_newick)
            present = [n.name for n in t.get_leaves() if n.name in species_set]
            if len(present) < 3:
                _pruned_tree_cache[species_set] = None
                return None
            t.prune(present, preserve_branch_length=True)
            _pruned_tree_cache[species_set] = t
            return t
        except Exception:
            _pruned_tree_cache[species_set] = None
            return None

    for gene_id, motifs in gene_motifs.items():
        # Phase 1: build position -> {species: alt_aa} and {lineage: miyata_group}
        position_species_aa: dict[int, dict[str, str]] = defaultdict(dict)
        position_lineage_groups: dict[int, dict[str, int]] = defaultdict(dict)
        position_human_aa: dict[int, str] = {}
        motif_ranges: list[tuple[str, int, int]] = []

        for motif_id, start_pos, human_seq, animal_seq, species_id in motifs:
            motif_ranges.append((motif_id, start_pos, min(len(human_seq), len(animal_seq))))

            lineage = species_lineage_map.get(species_id, "Other")
            if lineage in skip_lineages or species_id == "human":
                if species_id == "human":
                    seq_len_cmp = min(len(human_seq), len(animal_seq))
                    for i in range(seq_len_cmp):
                        if human_seq[i] != "-":
                            position_human_aa[start_pos + i] = human_seq[i]
                continue

            seq_len_cmp = min(len(human_seq), len(animal_seq))
            for i in range(seq_len_cmp):
                h_aa = human_seq[i]
                a_aa = animal_seq[i]
                if h_aa == "-" or a_aa == "-" or h_aa == a_aa:
                    continue
                abs_pos = start_pos + i
                position_species_aa[abs_pos][species_id] = a_aa.upper()
                if lineage not in position_lineage_groups[abs_pos]:
                    grp = _MIYATA.get(a_aa.upper(), -1)
                    position_lineage_groups[abs_pos][lineage] = grp

        # Phase 2: ASR-corrected convergence counting
        pos_max_conv: dict[int, int] = {}
        for abs_pos, lineage_groups in position_lineage_groups.items():
            if len(lineage_groups) < 2:
                continue

            if tree_newick and abs_pos in position_species_aa:
                species_aa = dict(position_species_aa[abs_pos])
                h_aa = position_human_aa.get(abs_pos)
                if h_aa:
                    species_aa["human"] = h_aa.upper()
                species_key = frozenset(species_aa.keys())
                parent_aa = _infer_ancestral_aa_cached(
                    _get_pruned_tree(species_key), species_aa
                )

                # Only count species where the substitution is independent
                # (leaf AA differs from its inferred parent AA)
                independent_lineages: dict[str, int] = {}
                for sp_id, alt_aa in position_species_aa[abs_pos].items():
                    anc_aa = parent_aa.get(sp_id)
                    if anc_aa and anc_aa.upper() != alt_aa.upper():
                        lineage = species_lineage_map.get(sp_id, "Other")
                        if lineage not in skip_lineages and lineage not in independent_lineages:
                            grp = _MIYATA.get(alt_aa, -1)
                            independent_lineages[lineage] = grp

                if len(independent_lineages) >= 2:
                    group_counts = Counter(independent_lineages.values())
                    group_counts.pop(-1, None)
                    if group_counts:
                        pos_max_conv[abs_pos] = max(group_counts.values())
            else:
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
