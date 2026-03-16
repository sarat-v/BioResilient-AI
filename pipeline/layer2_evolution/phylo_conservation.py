"""Step 3d — Phylogenetic conservation scoring with phyloP / PhastCons.

Builds a multi-species nucleotide alignment (MSA) for each gene region from the
NucleotideRegion table, then runs phyloP or PhastCons to compute evolutionary
conservation scores.  Scores indicate:
  - High positive value → strong purifying selection (functionally constrained)
  - Large negative value → accelerated evolution (possible adaptive pressure)
  - Near zero → evolving neutrally

Tools (checked in order):
  1. phyloP (PHAST package) — preferred; per-site log-odds conservation score
  2. phastCons (PHAST package) — alternative; posterior probability of conservation

Both tools require a species tree (the IQ-TREE2 treefile from step 5).
If neither tool is installed, all scores are left NULL and a WARN is logged.

Parallelisation: each gene's three regions (cds, promoter, downstream) are scored
by independent subprocess invocations.  ProcessPoolExecutor runs multiple genes
in parallel so all CPU cores are used simultaneously.

Scope: restricted to genes in tier_gene_ids (Tier1 + Tier2 after step9).
       When called before step9 (run_pipeline order), all genes are scored.

Entry point: run_phylo_conservation(tier_gene_ids, treefile) -> int
"""

import logging
import os
import re
import shutil
import subprocess
import tempfile
import uuid
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from db.models import Gene, NucleotideRegion, PhyloConservationScore
from db.session import get_session
from pipeline.config import get_tool_config

log = logging.getLogger(__name__)

_REGION_TYPES = ("cds", "promoter", "downstream")
_REGION_TO_FIELD = {
    "cds":        "cds_phylo_score",
    "promoter":   "promoter_phylo_score",
    "downstream": "downstream_phylo_score",
}

# Minimum number of species required to run phyloP (tool needs ≥3 for a tree model)
_MIN_SPECIES_FOR_PHYLOP = 3


# ─────────────────────────────────────────────────────────────────────────────
# Tool detection
# ─────────────────────────────────────────────────────────────────────────────

def _detect_phylo_tool() -> str:
    """Return 'phyloP', 'phastCons', or '' if neither is available."""
    for tool in ("phyloP", "phastCons"):
        if shutil.which(tool):
            return tool
    return ""


# ─────────────────────────────────────────────────────────────────────────────
# MSA building
# ─────────────────────────────────────────────────────────────────────────────

def build_msa(gene_id: str, region_type: str, session) -> Optional[Path]:
    """Build a FASTA MSA for a gene region across all available species.

    Sequences come from NucleotideRegion. Sequences of different lengths are
    padded with '-' to the length of the longest sequence (crude padding —
    sufficient for phyloP's per-region average score computation).

    Returns path to a temporary FASTA file, or None if < _MIN_SPECIES_FOR_PHYLOP
    sequences are available.  The caller is responsible for cleaning up.
    """
    regions = (
        session.query(NucleotideRegion)
        .filter_by(gene_id=gene_id, region_type=region_type)
        .all()
    )
    if not regions:
        return None

    # Filter to non-empty sequences
    records = [(r.species_id, r.sequence) for r in regions if r.sequence]
    if len(records) < _MIN_SPECIES_FOR_PHYLOP:
        return None

    max_len = max(len(seq) for _, seq in records)
    fasta_records = [
        SeqRecord(Seq(seq.ljust(max_len, "-")), id=sid, description="")
        for sid, seq in records
    ]

    tmp = tempfile.NamedTemporaryFile(
        suffix=".fa", delete=False, mode="w", prefix=f"msa_{gene_id[:8]}_{region_type}_"
    )
    SeqIO.write(fasta_records, tmp, "fasta")
    tmp.close()
    return Path(tmp.name)


def build_msa_from_records(
    gene_id: str,
    region_type: str,
    records: list[tuple[str, str]],
) -> Optional[Path]:
    """Build a FASTA MSA from pre-fetched (species_id, sequence) tuples.

    Used by workers to avoid opening a DB connection per-thread.
    Returns path to a temporary FASTA file, or None if too few species.
    """
    if len(records) < _MIN_SPECIES_FOR_PHYLOP:
        return None
    max_len = max(len(seq) for _, seq in records)
    fasta_records = [
        SeqRecord(Seq(seq.ljust(max_len, "-")), id=sid, description="")
        for sid, seq in records
    ]
    tmp = tempfile.NamedTemporaryFile(
        suffix=".fa", delete=False, mode="w", prefix=f"msa_{gene_id[:8]}_{region_type}_"
    )
    SeqIO.write(fasta_records, tmp, "fasta")
    tmp.close()
    return Path(tmp.name)


# ─────────────────────────────────────────────────────────────────────────────
# phyloP / phastCons runners
# ─────────────────────────────────────────────────────────────────────────────

def _newick_with_branch_lengths(treefile: Path) -> Optional[Path]:
    """Return a Newick file with branch lengths stripped of bootstrap labels.

    IQ-TREE2 produces a tree with bootstrap labels — phyloFit needs a clean
    Newick with only branch lengths.
    """
    try:
        newick = treefile.read_text().strip()
        # Remove bootstrap labels: digits just before colons (e.g. )85: → ):)
        newick_clean = __import__("re").sub(r"\)(\d+\.?\d*):", "):", newick)
        tmp = tempfile.NamedTemporaryFile(suffix=".nwk", delete=False, mode="w")
        tmp.write(newick_clean + "\n")
        tmp.close()
        return Path(tmp.name)
    except Exception as exc:
        log.debug("Tree cleaning failed: %s", exc)
        return None


def _fit_neutral_model(treefile: Path, concat_aln: Path) -> Optional[Path]:
    """Run phyloFit to produce a neutral .mod file for phyloP.

    phyloP requires a pre-fitted neutral evolutionary model (.mod file),
    not a raw Newick tree.  We fit the model once using the concatenated
    alignment from step 5 and cache it alongside the species tree.

    Args:
        treefile: IQ-TREE2 species treefile (Newick with branch lengths).
        concat_aln: Concatenated protein alignment FASTA from step 5.

    Returns:
        Path to the fitted .mod file, or None if phyloFit fails.
    """
    if not shutil.which("phyloFit"):
        log.warning("phyloFit not found — cannot fit neutral model for phyloP.")
        return None

    mod_path = treefile.parent / "neutral.mod"
    if mod_path.exists():
        log.info("Neutral model already exists at %s — skipping phyloFit.", mod_path)
        return mod_path

    nwk_path = _newick_with_branch_lengths(treefile)
    if nwk_path is None:
        return None

    try:
        log.info("Fitting neutral model with phyloFit (this may take a few minutes)...")
        result = subprocess.run(
            [
                "phyloFit",
                "--tree", str(nwk_path),
                "--msa-format", "FASTA",
                "--subst-mod", "SSREV",
                "--out-root", str(treefile.parent / "neutral"),
                str(concat_aln),
            ],
            capture_output=True, text=True, timeout=1800,
        )
        if result.returncode == 0 and mod_path.exists():
            log.info("Neutral model written to %s", mod_path)
            return mod_path
        log.warning("phyloFit failed (exit %d): %s", result.returncode, result.stderr[:400])
    except Exception as exc:
        log.warning("phyloFit error: %s", exc)
    finally:
        try:
            nwk_path.unlink(missing_ok=True)
        except Exception:
            pass
    return None


def _parse_phylop_output(stdout: str) -> Optional[float]:
    """Parse phyloP --wig-scores output and return the mean conservation score."""
    scores = []
    for line in stdout.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            scores.append(float(line.split()[-1]))
        except ValueError:
            continue
    if not scores:
        return None
    return round(sum(scores) / len(scores), 4)


def _parse_phastcons_output(stdout: str) -> Optional[float]:
    """Parse phastCons output (posterior probabilities) and return mean score."""
    scores = []
    for line in stdout.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        for token in line.split():
            try:
                val = float(token)
                if 0.0 <= val <= 1.0:
                    scores.append(val)
            except ValueError:
                continue
    if not scores:
        return None
    return round(sum(scores) / len(scores), 4)


def _prune_newick_ete3(newick: str, keep: set[str]) -> Optional[str]:
    """Prune Newick using ETE3. Returns pruned string or None if ETE3 unavailable."""
    try:
        from ete3 import Tree  # type: ignore
        t = Tree(newick, format=1)
        t.prune(list(keep), preserve_branch_length=True)
        return t.write(format=1)
    except Exception:
        return None


def _pruned_mod_via_tree_doctor(mod_path: Path, species_in_msa: list[str]) -> Optional[Path]:
    """Use PHAST's tree_doctor to prune a .mod file to only species_in_msa.

    tree_doctor writes the pruned .mod to stdout — capture it and write to a
    temp file for use by phyloP.
    """
    if not shutil.which("tree_doctor"):
        return None
    try:
        result = subprocess.run(
            ["tree_doctor",
             "--prune-all-but", ",".join(species_in_msa),
             str(mod_path)],
            capture_output=True, text=True, timeout=30,
        )
        if result.returncode == 0 and result.stdout.strip():
            tmp = tempfile.NamedTemporaryFile(suffix=".mod", delete=False, mode="w")
            tmp.write(result.stdout)
            tmp.close()
            return Path(tmp.name)
        log.debug("tree_doctor failed (exit %d): %s", result.returncode, result.stderr[:200])
    except Exception as exc:
        log.debug("tree_doctor error: %s", exc)
    return None


def _pruned_mod(mod_path: Path, species_in_msa: list[str]) -> Optional[Path]:
    """Write a copy of the .mod file with the TREE pruned to only species_in_msa.

    phyloP errors with "no match for leaves of tree in alignment" when the MSA
    contains a subset of species but the model tree has all 18 leaves.

    Tries (in order): tree_doctor (PHAST, most reliable) → ETE3 + manual edit.
    """
    # Preferred: tree_doctor handles .mod files natively
    result = _pruned_mod_via_tree_doctor(mod_path, species_in_msa)
    if result is not None:
        return result

    # Fallback: ETE3 + manual TREE line replacement
    try:
        mod_text = mod_path.read_text()
        tree_match = re.search(r"^TREE:\s*(.+)$", mod_text, re.MULTILINE)
        if not tree_match:
            return None
        newick = tree_match.group(1).strip()
        keep = set(species_in_msa)
        pruned_newick = _prune_newick_ete3(newick, keep)
        if pruned_newick is None:
            log.debug("_pruned_mod: both tree_doctor and ETE3 failed for %d species", len(species_in_msa))
            return None
        pruned_mod = (mod_text[:tree_match.start(1)]
                      + pruned_newick
                      + mod_text[tree_match.end(1):])
        tmp = tempfile.NamedTemporaryFile(suffix=".mod", delete=False, mode="w")
        tmp.write(pruned_mod)
        tmp.close()
        return Path(tmp.name)
    except Exception as exc:
        log.debug("_pruned_mod ETE3 fallback failed: %s", exc)
        return None


def run_phylop(msa_path: Path, mod_path: Path) -> Optional[float]:
    """Run phyloP on the MSA using a pre-fitted neutral model (.mod file).

    The neutral model contains all 18 species; each per-gene MSA may have
    fewer.  We prune the model tree to only the species in the MSA before
    calling phyloP to avoid "no match for leaves of tree in alignment" errors.
    """
    # Extract species present in this MSA
    species_in_msa = [
        rec.id for rec in SeqIO.parse(str(msa_path), "fasta")
    ]
    if not species_in_msa:
        return None

    pruned_path: Optional[Path] = None
    effective_mod = mod_path
    try:
        pruned_path = _pruned_mod(mod_path, species_in_msa)
        if pruned_path is not None:
            effective_mod = pruned_path

        result = subprocess.run(
            ["phyloP", "--wig-scores", "--msa-format", "FASTA",
             str(effective_mod), str(msa_path)],
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode == 0:
            return _parse_phylop_output(result.stdout)
        log.warning("phyloP non-zero exit (%d) for %s: %s",
                    result.returncode, msa_path.name, result.stderr[:300])
    except Exception as exc:
        log.debug("phyloP execution error: %s", exc)
    finally:
        if pruned_path is not None:
            try:
                pruned_path.unlink(missing_ok=True)
            except Exception:
                pass
    return None


def run_phastcons(msa_path: Path, tree_path: Path) -> Optional[float]:
    """Run phastCons on the MSA and return the mean score, or None on failure."""
    with tempfile.TemporaryDirectory() as tmpdir:
        score_file = Path(tmpdir) / "scores.wig"
        try:
            result = subprocess.run(
                ["phastCons", str(msa_path), str(tree_path),
                 "--score", "--no-post-probs",
                 "--score-file", str(score_file)],
                capture_output=True, text=True, timeout=300,
            )
            if result.returncode == 0 and score_file.exists():
                return _parse_phastcons_output(score_file.read_text())
            log.debug("phastCons non-zero exit (%d)", result.returncode)
        except Exception as exc:
            log.debug("phastCons execution error: %s", exc)
    return None


def score_region(msa_path: Path, mod_path: Path, tool: str) -> Optional[float]:
    """Run the selected phylogenetic tool and return a mean conservation score."""
    if tool == "phyloP":
        return run_phylop(msa_path, mod_path)
    if tool == "phastCons":
        return run_phastcons(msa_path, mod_path)
    return None


# ─────────────────────────────────────────────────────────────────────────────
# DB upsert
# ─────────────────────────────────────────────────────────────────────────────

def _upsert_phylo_score(session, gene_id: str, scores: dict[str, Optional[float]]) -> None:
    existing = session.get(PhyloConservationScore, gene_id)
    if existing is None:
        existing = PhyloConservationScore(gene_id=gene_id)
        session.add(existing)
    for field, value in scores.items():
        if value is not None:
            setattr(existing, field, value)


def _score_gene_worker(args: tuple) -> tuple[str, dict[str, Optional[float]]]:
    """Thread worker: build MSAs from pre-fetched sequences and run phyloP/phastCons.

    Accepts pre-fetched sequences to avoid DB connection pool exhaustion when
    running 60 concurrent threads.  No DB access in this function.

    args: (gene_id, {region_type: [(species_id, sequence), ...]}, mod_path_str, tool)
    """
    gene_id, region_seqs, mod_path_str, tool = args
    mod_path = Path(mod_path_str)

    region_scores: dict[str, Optional[float]] = {}
    tmp_files: list[Path] = []

    try:
        for region_type in _REGION_TYPES:
            records = region_seqs.get(region_type, [])
            msa_path = build_msa_from_records(gene_id, region_type, records)
            if msa_path is None:
                continue
            tmp_files.append(msa_path)
            score = score_region(msa_path, mod_path, tool)
            field = _REGION_TO_FIELD[region_type]
            region_scores[field] = score
    finally:
        for f in tmp_files:
            try:
                f.unlink(missing_ok=True)
            except Exception:
                pass

    return gene_id, region_scores


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def run_phylo_conservation(
    tier_gene_ids: list[str] | None,
    treefile: Path,
) -> int:
    """Build MSAs and score conservation with phyloP/phastCons for each gene.

    Runs all gene-level subprocess invocations in parallel using
    ProcessPoolExecutor so all CPU cores are utilised simultaneously.

    Args:
        tier_gene_ids: List of gene IDs to score. If None, all genes with
                       NucleotideRegion data are scored.
        treefile: Path to the IQ-TREE2 species treefile (step 5 output).

    Returns count of genes with at least one phylo score written.
    """
    tool = _detect_phylo_tool()
    if not tool:
        log.warning(
            "phyloP and phastCons not found. Phylogenetic conservation scores will be NULL. "
            "Install the PHAST package: conda install -c bioconda phast"
        )
        return 0

    if not treefile.exists():
        log.warning("Treefile not found at %s — phylo conservation skipped.", treefile)
        return 0

    # Fit a neutral evolutionary model with phyloFit — required by phyloP.
    # The concatenated alignment from step 5 lives alongside the treefile.
    concat_aln = treefile.parent / "concat_alignment.faa"
    if not concat_aln.exists():
        log.warning(
            "Concatenated alignment not found at %s — cannot fit neutral model. "
            "Phylogenetic conservation scores will be NULL.",
            concat_aln,
        )
        return 0

    mod_path = _fit_neutral_model(treefile, concat_aln)
    if mod_path is None:
        log.warning("Could not fit neutral model for phyloP — skipped.")
        return 0

    with get_session() as session:
        if tier_gene_ids:
            gene_ids = [g.id for g in session.query(Gene).filter(Gene.id.in_(tier_gene_ids)).all()]
        else:
            gene_ids = [
                r[0] for r in session.query(NucleotideRegion.gene_id).distinct().all()
            ]

        # Bulk-fetch all nucleotide sequences in one query — avoids opening a
        # DB connection per worker thread (which exhausts the connection pool).
        log.info("Pre-fetching nucleotide sequences for %d genes...", len(gene_ids))
        rows = (
            session.query(
                NucleotideRegion.gene_id,
                NucleotideRegion.region_type,
                NucleotideRegion.species_id,
                NucleotideRegion.sequence,
            )
            .filter(NucleotideRegion.gene_id.in_(gene_ids))
            .filter(NucleotideRegion.sequence.isnot(None))
            .all()
        )

        # Build {gene_id: {region_type: [(species_id, sequence)]}}
        gene_region_seqs: dict[str, dict[str, list[tuple[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for gene_id, region_type, species_id, sequence in rows:
        gene_region_seqs[gene_id][region_type].append((species_id, sequence))

    log.info("Phylogenetic conservation scoring with %s for %d genes (parallel)...",
             tool, len(gene_ids))

    n_workers = int(get_tool_config().get("phylop_workers", min(60, max(1, (os.cpu_count() or 4) * 4))))
    work_items = [
        (gid, dict(gene_region_seqs.get(gid, {})), str(mod_path), tool)
        for gid in gene_ids
    ]
    done = 0
    total = len(work_items)

    _DB_FLUSH_INTERVAL = 1000  # write to DB every N genes to show incremental progress

    genes_scored = 0
    errors = 0
    pending: dict[str, dict[str, Optional[float]]] = {}

    # ThreadPoolExecutor: phyloP/phastCons are external subprocesses — they
    # release the GIL and are I/O-bound, so threads outperform processes here.
    # No DB access in workers — sequences pre-fetched above.
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_score_gene_worker, item): item[0] for item in work_items}
        for future in as_completed(futures):
            try:
                gene_id, region_scores = future.result()
                if any(v is not None for v in region_scores.values()):
                    pending[gene_id] = region_scores
            except Exception as exc:
                errors += 1
                gid = futures[future]
                if errors <= 3:
                    log.warning("  phylo_conservation: worker error for gene %s: %s",
                                gid, exc, exc_info=True)
                else:
                    log.debug("  phylo_conservation: worker error for gene %s: %s", gid, exc)
            done += 1
            if done % 100 == 0 or done == total:
                log.info("  phylo_conservation: %d / %d genes (%.0f%%) — %d scored, %d errors.",
                         done, total, 100 * done / total, genes_scored + len(pending), errors)
            # Flush pending results to DB periodically
            if len(pending) >= _DB_FLUSH_INTERVAL:
                with get_session() as session:
                    for gid, region_scores in pending.items():
                        _upsert_phylo_score(session, gid, region_scores)
                    session.commit()
                genes_scored += len(pending)
                pending.clear()

    # Final flush
    if pending:
        with get_session() as session:
            for gid, region_scores in pending.items():
                _upsert_phylo_score(session, gid, region_scores)
            session.commit()
        genes_scored += len(pending)
        pending.clear()

    if errors:
        log.warning("  phylo_conservation: %d genes failed scoring (logged at DEBUG level).", errors)
    log.info("Phylogenetic conservation scoring complete: %d genes scored", genes_scored)
    return genes_scored
