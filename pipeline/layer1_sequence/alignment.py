"""Step 4a — MAFFT alignment wrapper.

For each ortholog group (OrthoGroup), collects all protein sequences and runs
MAFFT to produce a multiple sequence alignment (MSA) in aligned FASTA format.

Gate 1 (funnel): filter_orthogroups_by_global_identity() reduces orthogroups
before MAFFT to those where ≥2 resilient species diverge from human by >15%.
"""

import logging
import os
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

from db.models import Ortholog
from db.session import get_session
from pipeline.config import get_local_storage_root, get_tool_config

log = logging.getLogger(__name__)

try:
    import numpy as np
    _NUMPY_AVAILABLE = True
except ImportError:
    _NUMPY_AVAILABLE = False


def _alignments_dir() -> Path:
    root = Path(get_local_storage_root())
    d = root / "alignments"
    d.mkdir(parents=True, exist_ok=True)
    return d


def align_orthogroup(og_id: str, sequences: dict[str, str], threads: int = 1) -> Optional[dict[str, str]]:
    """Run MAFFT then trimAl-mask on a set of sequences for one orthogroup.

    Args:
        og_id: OrthoGroup identifier (used for output file naming).
        sequences: {label: protein_sequence} — label should be unique.
        threads: CPU threads for MAFFT per worker.

    Returns:
        {label: aligned_sequence} (trimAl-masked when available) or None on failure.
    """
    if len(sequences) < 2:
        log.debug("  Skipping %s — only %d sequence(s)", og_id, len(sequences))
        return None

    cfg = get_tool_config()
    threads = cfg.get("mafft_threads", threads)

    aln_dir = _alignments_dir()
    out_path = aln_dir / f"{og_id}.afa"

    if not (out_path.exists() and out_path.stat().st_size > 0):
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
                out_path.unlink(missing_ok=True)
                return None
        except Exception as exc:
            log.error("  MAFFT error for %s: %s", og_id, exc)
            out_path.unlink(missing_ok=True)
            return None
        finally:
            Path(tmp_path).unlink(missing_ok=True)

    masked_path = mask_alignment_trimal(og_id, out_path)
    if masked_path:
        return _parse_aligned_fasta(masked_path)
    return _parse_aligned_fasta(out_path)


def mask_alignment_trimal(og_id: str, aln_path: Path) -> Optional[Path]:
    """Run trimAl to remove poorly-aligned or gap-heavy columns.

    Uses --gt 0.1 (≥90 % of sequences must have a residue) and --cons 60
    (column conservation threshold 60 %) which together remove the majority
    of misaligned/gapped columns while keeping informative positions.

    Writes the masked alignment alongside the original as <og_id>.masked.afa.
    Returns the masked path on success, or None if trimAl is not available
    or the masked file would be empty.
    """
    masked_path = aln_path.parent / f"{og_id}.masked.afa"
    if masked_path.exists() and masked_path.stat().st_size > 0:
        return masked_path

    cmd = [
        "trimal",
        "-in", str(aln_path),
        "-out", str(masked_path),
        "-gt", "0.1",
        "-cons", "60",
        "-fasta",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0 or not masked_path.exists() or masked_path.stat().st_size == 0:
            log.debug("  trimAl unavailable or empty result for %s; using raw MAFFT output", og_id)
            masked_path.unlink(missing_ok=True)
            return None

        raw = _parse_aligned_fasta(aln_path)
        masked = _parse_aligned_fasta(masked_path)
        if not masked:
            masked_path.unlink(missing_ok=True)
            return None

        raw_len = len(next(iter(raw.values()), ""))
        masked_len = len(next(iter(masked.values()), ""))
        if raw_len > 0:
            frac_removed = (raw_len - masked_len) / raw_len
            log.debug("  trimAl %s: %d → %d cols (%.1f %% removed)", og_id, raw_len, masked_len, frac_removed * 100)
            if frac_removed > 0.80:
                log.warning("  trimAl removed >80%% of columns for %s — falling back to raw alignment", og_id)
                masked_path.unlink(missing_ok=True)
                return None

        return masked_path
    except FileNotFoundError:
        log.debug("  trimAl not found on PATH; skipping alignment masking")
        return None
    except Exception as exc:
        log.warning("  trimAl error for %s: %s", og_id, exc)
        masked_path.unlink(missing_ok=True)
        return None


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
    Uses numpy for fast vectorized comparison.
    """
    if len(seq_a) != len(seq_b):
        raise ValueError("Sequences must be aligned (same length).")

    try:
        if _NUMPY_AVAILABLE:
            a = np.frombuffer(seq_a.encode(), dtype=np.uint8)
            b = np.frombuffer(seq_b.encode(), dtype=np.uint8)
            gap = ord("-")
            both_gap = (a == gap) & (b == gap)
            comparable = int((~both_gap).sum())
            if comparable == 0:
                return 0.0
            matches = int(((a == b) & ~both_gap).sum())
            return (matches / comparable) * 100.0
    except Exception as exc:
        log.debug("numpy identity calculation failed, falling back to pure Python: %s", exc)
    # Pure Python fallback
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


def _align_worker(args: tuple) -> tuple[str, Optional[dict[str, str]]]:
    """Top-level helper for multiprocessing (must be picklable)."""
    og_id, seqs, threads = args
    return og_id, align_orthogroup(og_id, seqs, threads=threads)


def align_all_orthogroups(
    orthogroups: dict[str, dict[str, str]],
) -> dict[str, dict[str, str]]:
    """Align all orthogroups in parallel and return a map of og_id → aligned sequences.

    Uses ProcessPoolExecutor with one worker per logical CPU.  Each worker runs
    MAFFT with 1 thread so the total CPU usage is n_workers × 1 ≈ all cores.

    Args:
        orthogroups: {og_id: {label: sequence}} — output from load_orthogroup_sequences().

    Returns:
        {og_id: {label: aligned_sequence}}
    """
    n_cpu = os.cpu_count() or 8
    # Reserve 2 cores for OS/DB overhead; each MAFFT uses 1 thread
    n_workers = max(1, n_cpu - 2)
    mafft_threads_per_worker = 1

    items = list(orthogroups.items())
    total = len(items)
    aligned: dict[str, dict[str, str]] = {}
    done = 0

    log.info("Aligning %d orthogroups with %d parallel workers (1 MAFFT thread each).", total, n_workers)

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(_align_worker, (og_id, seqs, mafft_threads_per_worker)): og_id
            for og_id, seqs in items
        }
        for future in as_completed(futures, timeout=600):
            try:
                og_id, result = future.result(timeout=120)
            except Exception as exc:
                og_id = futures[future]
                log.warning("Alignment worker failed for %s: %s", og_id, exc)
                done += 1
                continue
            done += 1
            if result:
                aligned[og_id] = result
            if done % 500 == 0 or done == total:
                log.info("  Aligned %d / %d orthogroups (%.0f%%).", done, total, 100 * done / total)

    log.info("Aligned %d / %d orthogroups.", len(aligned), total)
    return aligned


def _quick_pairwise_identity(seq_a: str, seq_b: str) -> float:
    """Compute approximate pairwise identity (0–100) using global alignment.
    Used by Gate 1 filter; no need for full MSA.
    """
    if not seq_a or not seq_b:
        return 0.0
    # Truncate very long sequences to first 300 aa for speed
    seq_a = seq_a[:300]
    seq_b = seq_b[:300]
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.1
    # Use score_only to avoid OverflowError from enumerating all optimal alignments
    score = aligner.score(seq_a, seq_b)
    max_len = max(len(seq_a), len(seq_b))
    if max_len == 0:
        return 0.0
    return 100.0 * score / max_len


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


def load_orthogroup_sequences_from_db(one_to_one_only: bool = True) -> dict[str, dict[str, str]]:
    """Load orthogroup sequences from the database.

    Args:
        one_to_one_only: When True (default), only load orthogroups flagged as 1:1
            (no paralogs). This prevents false convergence signals from paralog
            confusion. Set False only for exploratory analysis.

    Returns {og_id: {"{species_id}|{protein_id}": sequence}}.
    Only includes orthogroups with at least one human + one non-human protein.
    """
    orthogroups: dict[str, dict[str, str]] = {}

    with get_session() as session:
        q = session.query(Ortholog).filter(Ortholog.protein_seq.isnot(None), Ortholog.protein_seq != "")
        if one_to_one_only:
            q = q.filter(Ortholog.is_one_to_one.is_(True))
        orthologs = q.all()
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
    if one_to_one_only:
        log.info("  %d 1:1 orthogroups loaded (paralogs excluded).", len(filtered))
    else:
        log.info("  %d orthogroups loaded with human + ≥1 other species.", len(filtered))
    return filtered
