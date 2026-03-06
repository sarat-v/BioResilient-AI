"""Step 4a — MAFFT alignment wrapper.

For each ortholog group (OrthoGroup), collects all protein sequences and runs
MAFFT to produce a multiple sequence alignment (MSA) in aligned FASTA format.

Gate 1 (funnel): filter_orthogroups_by_global_identity() reduces orthogroups
before MAFFT to those where ≥2 resilient species diverge from human by >15%.
"""

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

from db.models import Ortholog
from db.session import get_session
from pipeline.config import get_storage_root, get_tool_config

log = logging.getLogger(__name__)


def _alignments_dir() -> Path:
    root = Path(get_storage_root())
    d = root / "alignments"
    d.mkdir(parents=True, exist_ok=True)
    return d


def align_orthogroup(og_id: str, sequences: dict[str, str], threads: int = 4) -> Optional[dict[str, str]]:
    """Run MAFFT on a set of sequences for one orthogroup.

    Args:
        og_id: OrthoGroup identifier (used for output file naming).
        sequences: {label: protein_sequence} — label should be unique.
        threads: CPU threads for MAFFT.

    Returns:
        {label: aligned_sequence} or None on failure.
    """
    if len(sequences) < 2:
        log.debug("  Skipping %s — only %d sequence(s)", og_id, len(sequences))
        return None

    cfg = get_tool_config()
    threads = cfg.get("mafft_threads", threads)

    aln_dir = _alignments_dir()
    out_path = aln_dir / f"{og_id}.afa"

    if out_path.exists() and out_path.stat().st_size > 0:
        return _parse_aligned_fasta(out_path)

    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as tmp:
        for label, seq in sequences.items():
            tmp.write(f">{label}\n{seq}\n")
        tmp_path = tmp.name

    cmd = [
        "mafft",
        "--auto",
        "--thread", str(threads),
        "--quiet",
        tmp_path,
    ]

    try:
        with open(out_path, "w") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, check=False)

        if result.returncode != 0:
            log.warning("  MAFFT failed for %s: %s", og_id, result.stderr.decode())
            return None

        return _parse_aligned_fasta(out_path)

    except Exception as exc:
        log.error("  MAFFT error for %s: %s", og_id, exc)
        return None
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _parse_aligned_fasta(path: Path) -> dict[str, str]:
    """Parse aligned FASTA. Returns {seq_id: aligned_sequence}."""
    result = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        result[rec.id] = str(rec.seq)
    return result


def calculate_sequence_identity(seq_a: str, seq_b: str) -> float:
    """Compute pairwise sequence identity between two aligned sequences.

    Returns a value in [0, 100] (percentage).
    Gaps in both sequences are excluded from the denominator.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("Sequences must be aligned (same length).")

    matches = 0
    comparable = 0
    for a, b in zip(seq_a, seq_b):
        if a == "-" and b == "-":
            continue
        comparable += 1
        if a == b:
            matches += 1

    if comparable == 0:
        return 0.0
    return (matches / comparable) * 100.0


def align_all_orthogroups(
    orthogroups: dict[str, dict[str, str]],
) -> dict[str, dict[str, str]]:
    """Align all orthogroups and return a map of og_id → aligned sequences.

    Args:
        orthogroups: {og_id: {label: sequence}} — output from load_orthogroup_sequences().

    Returns:
        {og_id: {label: aligned_sequence}}
    """
    aligned: dict[str, dict[str, str]] = {}
    for og_id, seqs in orthogroups.items():
        result = align_orthogroup(og_id, seqs)
        if result:
            aligned[og_id] = result
    log.info("Aligned %d / %d orthogroups.", len(aligned), len(orthogroups))
    return aligned


def _quick_pairwise_identity(seq_a: str, seq_b: str) -> float:
    """Compute approximate pairwise identity (0–100) using global alignment.
    Used by Gate 1 filter; no need for full MSA.
    """
    if not seq_a or not seq_b:
        return 0.0
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(seq_a, seq_b)
    if not alignments:
        return 0.0
    aln = alignments[0]
    aln1, aln2 = str(aln[0]), str(aln[1])
    matches = sum(1 for a, b in zip(aln1, aln2) if a == b and a != "-")
    comparable = sum(1 for a, b in zip(aln1, aln2) if a != "-" or b != "-")
    if comparable == 0:
        return 0.0
    return 100.0 * matches / comparable


def filter_orthogroups_by_global_identity(
    orthogroups: dict[str, dict[str, str]],
    min_divergent_species: int = 2,
    divergence_pct_min: float = 15.0,
) -> dict[str, dict[str, str]]:
    """Gate 1: keep only orthogroups where ≥ min_divergent_species show divergence > divergence_pct_min from human.

    Uses quick pairwise identity (no MAFFT). Expected reduction: ~15k → 2k–4k orthogroups.
    """
    identity_max = 100.0 - divergence_pct_min  # e.g. 85% identity → 15% divergence
    filtered: dict[str, dict[str, str]] = {}
    for og_id, seqs in orthogroups.items():
        human_label = next((k for k in seqs if "human" in k.lower()), None)
        if not human_label:
            continue
        human_seq = seqs[human_label]
        count_divergent = 0
        for label, seq in seqs.items():
            if "human" in label.lower():
                continue
            identity = _quick_pairwise_identity(human_seq, seq)
            if identity < identity_max:
                count_divergent += 1
        if count_divergent >= min_divergent_species:
            filtered[og_id] = seqs
    log.info(
        "  Gate 1: %d orthogroups pass global identity filter (≥%d species with >%.0f%% divergence).",
        len(filtered),
        min_divergent_species,
        divergence_pct_min,
    )
    return filtered


def load_orthogroup_sequences_from_db() -> dict[str, dict[str, str]]:
    """Load orthogroup sequences from the database.

    Returns {og_id: {"{species_id}|{protein_id}": sequence}}.
    Only includes orthogroups with at least one human + one non-human protein.
    """
    orthogroups: dict[str, dict[str, str]] = {}

    with get_session() as session:
        orthologs = session.query(Ortholog).filter(Ortholog.protein_seq.isnot(None)).all()
        for o in orthologs:
            og = o.orthofinder_og
            if og not in orthogroups:
                orthogroups[og] = {}
            label = f"{o.species_id}|{o.protein_id}"
            orthogroups[og][label] = o.protein_seq

    # Filter: must have human + at least one other species
    filtered = {
        og_id: seqs
        for og_id, seqs in orthogroups.items()
        if any("human" in k for k in seqs)
        and sum(1 for k in seqs if "human" not in k) >= 1
    }
    log.info("  %d orthogroups loaded with human + ≥1 other species.", len(filtered))
    return filtered
