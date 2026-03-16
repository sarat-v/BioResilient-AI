"""Step 5 — Phylogenetic tree construction with IQ-TREE2.

Builds a species tree from the concatenated alignment of single-copy orthologs.
The tree is required before running HyPhy (Step 6).

Output: species.treefile in Newick format.
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pipeline.config import get_local_storage_root, get_tool_config

log = logging.getLogger(__name__)


def build_concatenated_alignment(
    aligned_orthogroups: dict[str, dict[str, str]],
    single_copy_only: bool = True,
    max_ogs: int = 500,
) -> Optional[dict[str, str]]:
    """Build a concatenated super-alignment from single-copy orthogroups.

    Single-copy: each species appears exactly once in the orthogroup.

    Args:
        aligned_orthogroups: {og_id: {label: aligned_sequence}}
            Labels are in the format "species_id|species_id|protein_id"
            (OrthoFinder reheadered format).
        single_copy_only: If True, only use orthogroups with exactly one protein per species.
        max_ogs: Maximum number of orthogroups to include. Capped to keep memory manageable
            for IQ-TREE — 500 OGs is more than sufficient for 18-species topology.

    Returns:
        {species_id: concatenated_alignment_string} or None if too few orthogroups.
    """
    species_ids = _collect_species_ids(aligned_orthogroups)
    if not species_ids:
        log.error("build_concatenated_alignment: no species found in aligned_orthogroups "
                  "(dict has %d entries)", len(aligned_orthogroups))
        return None

    n_species = len(species_ids)
    min_coverage = n_species * 0.8

    # Diagnostic counters
    n_multi_copy = 0
    n_low_coverage = 0

    valid_ogs = []
    for og_id, seqs in aligned_orthogroups.items():
        present_species: dict[str, int] = {}
        for label in seqs:
            sid = label.split("|")[0]
            present_species[sid] = present_species.get(sid, 0) + 1

        # Single-copy: each species appears at most once
        if single_copy_only and any(v > 1 for v in present_species.values()):
            n_multi_copy += 1
            continue

        if len(present_species) < min_coverage:
            n_low_coverage += 1
            continue

        valid_ogs.append(og_id)

    log.info(
        "OG filter: %d total, %d single-copy & ≥80%% coverage, "
        "%d rejected (multi-copy), %d rejected (low coverage, need ≥%.0f / %d species).",
        len(aligned_orthogroups), len(valid_ogs),
        n_multi_copy, n_low_coverage, min_coverage, n_species,
    )

    if len(valid_ogs) < 10:
        log.warning("Only %d valid orthogroups found — tree may be unreliable.", len(valid_ogs))
        if len(valid_ogs) == 0:
            return None

    # Cap OG count to keep IQ-TREE memory usage manageable.
    # 500 OGs × ~900 aa avg × 18 species ≈ 8M columns — more than enough for topology.
    if len(valid_ogs) > max_ogs:
        # Prefer OGs with highest species coverage (most informative)
        def _coverage(og_id: str) -> int:
            return len({lbl.split("|")[0] for lbl in aligned_orthogroups[og_id]})
        valid_ogs = sorted(valid_ogs, key=_coverage, reverse=True)[:max_ogs]
        log.info("Capped to top %d orthogroups by species coverage.", max_ogs)

    log.info("Building concatenated alignment from %d orthogroups...", len(valid_ogs))

    concat: dict[str, list[str]] = {sid: [] for sid in species_ids}

    for og_id in valid_ogs:
        seqs = aligned_orthogroups[og_id]
        aln_len = len(next(iter(seqs.values())))

        # Map labels to species_id
        species_seqs: dict[str, str] = {}
        for label, seq in seqs.items():
            sid = label.split("|")[0]
            species_seqs[sid] = seq

        for sid in species_ids:
            if sid in species_seqs:
                concat[sid].append(species_seqs[sid])
            else:
                # Missing species in this OG — fill with gaps
                concat[sid].append("-" * aln_len)

    return {sid: "".join(parts) for sid, parts in concat.items()}


def _collect_species_ids(aligned_orthogroups: dict[str, dict[str, str]]) -> set[str]:
    species = set()
    for seqs in aligned_orthogroups.values():
        for label in seqs:
            species.add(label.split("|")[0])
    return species


def _max_ogs() -> int:
    """Max OGs for concatenated alignment — configurable via iqtree_max_ogs in environment.yml."""
    return int(get_tool_config().get("iqtree_max_ogs", 500))


def run_iqtree(concat_alignment: dict[str, str]) -> Path:
    """Run IQ-TREE2 on the concatenated alignment.

    Args:
        concat_alignment: {species_id: aligned_sequence}

    Returns:
        Path to the .treefile (Newick format).
    """
    cfg = get_tool_config()
    threads = cfg.get("iqtree_threads", "AUTO")
    bootstrap = cfg.get("iqtree_bootstrap", 1000)
    # Default model: Q.INSECT+F+I+G4 — best BIC on empirical 18-species panel
    # (vertebrates + invertebrates). Skipping ModelFinder saves ~30 min per run.
    # Override via iqtree_model in environment.yml if needed.
    model = cfg.get("iqtree_model", "Q.INSECT+F+I+G4")

    # IQ-TREE's AUTO thread detection runs a memory-intensive benchmarking phase
    # that can OOM on large alignments. Pin to a safe default when AUTO is set.
    if str(threads).upper() == "AUTO":
        threads = 4

    root = Path(get_local_storage_root())
    tree_dir = root / "phylo"
    tree_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = tree_dir / "concat_alignment.faa"
    with open(fasta_path, "w") as f:
        for species_id, seq in concat_alignment.items():
            f.write(f">{species_id}\n{seq}\n")

    treefile = tree_dir / "species.treefile"
    if treefile.exists():
        log.info("Species tree already exists at %s — skipping IQ-TREE.", treefile)
        return treefile

    # Support both 'iqtree2' (bioconda newer) and 'iqtree' (bioconda older/conda)
    import shutil as _shutil
    iqtree_bin = "iqtree2" if _shutil.which("iqtree2") else "iqtree"

    cmd = [
        iqtree_bin,
        "-s", str(fasta_path),
        "-m", model,
        "-T", str(threads),
        "--mem", "20G",       # hard cap to avoid OOM on 30GB instances
        "-o", "human",
        "--prefix", str(tree_dir / "species"),
        "--redo",
    ]
    # Bootstrap requires ≥4 sequences — skip for small test runs
    if len(concat_alignment) >= 4:
        cmd += ["-bb", str(bootstrap)]

    log.info("Running IQ-TREE2: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=False, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"IQ-TREE2 exited with code {result.returncode}")

    if not treefile.exists():
        raise FileNotFoundError(f"IQ-TREE2 did not produce {treefile}")

    log.info("Species tree written to %s", treefile)
    return treefile


def load_tree(treefile: Path) -> str:
    """Read Newick tree string from file."""
    return treefile.read_text().strip()


def prune_tree_to_species(treefile: Path, species_ids: list[str]) -> str:
    """Prune the species tree to only the species present in a given orthogroup.

    Uses ETE3 if available, otherwise returns the full tree.
    """
    try:
        from ete3 import Tree

        t = Tree(str(treefile))
        # Remove leaves not in species_ids
        leaves_to_remove = [n for n in t.get_leaves() if n.name not in species_ids]
        for node in leaves_to_remove:
            node.detach()
        return t.write(format=1)

    except ImportError:
        log.debug("ete3 not available — using full tree without pruning.")
        return load_tree(treefile)
    except Exception as exc:
        log.warning("Tree pruning failed: %s — using full tree.", exc)
        return load_tree(treefile)
