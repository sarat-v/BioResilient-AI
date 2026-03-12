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

Scope: restricted to genes in tier_gene_ids (Tier1 + Tier2 after step9).
       When called before step9 (run_pipeline order), all genes are scored.

Entry point: run_phylo_conservation(tier_gene_ids, treefile) -> int
"""

import logging
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from db.models import Gene, NucleotideRegion, PhyloConservationScore
from db.session import get_session

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


# ─────────────────────────────────────────────────────────────────────────────
# phyloP / phastCons runners
# ─────────────────────────────────────────────────────────────────────────────

def _newick_with_branch_lengths(treefile: Path) -> Optional[Path]:
    """Return a Newick file with branch lengths, as required by PHAST tools.

    IQ-TREE2 produces a tree with bootstrap labels — phyloP needs a tree
    with only branch lengths.  We strip bootstrap values here.
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


def run_phylop(msa_path: Path, tree_path: Path) -> Optional[float]:
    """Run phyloP on the MSA and return the mean score, or None on failure."""
    try:
        result = subprocess.run(
            ["phyloP", "--wig-scores", str(tree_path), str(msa_path)],
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode == 0:
            return _parse_phylop_output(result.stdout)
        log.debug("phyloP non-zero exit (%d): %s", result.returncode, result.stderr[:200])
    except Exception as exc:
        log.debug("phyloP execution error: %s", exc)
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


def score_region(msa_path: Path, tree_path: Path, tool: str) -> Optional[float]:
    """Run the selected phylogenetic tool and return a mean conservation score."""
    if tool == "phyloP":
        return run_phylop(msa_path, tree_path)
    if tool == "phastCons":
        return run_phastcons(msa_path, tree_path)
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


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def run_phylo_conservation(
    tier_gene_ids: list[str] | None,
    treefile: Path,
) -> int:
    """Build MSAs and score conservation with phyloP/phastCons for each gene.

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

    tree_path = _newick_with_branch_lengths(treefile)
    if tree_path is None:
        log.warning("Could not prepare tree for phyloP — skipped.")
        return 0

    log.info("Phylogenetic conservation scoring with %s...", tool)
    genes_scored = 0
    tmp_files: list[Path] = []

    try:
        with get_session() as session:
            if tier_gene_ids:
                genes = session.query(Gene).filter(Gene.id.in_(tier_gene_ids)).all()
            else:
                # Fall back to all genes that have NucleotideRegion data
                gene_ids_with_regions = [
                    r[0] for r in session.query(NucleotideRegion.gene_id).distinct().all()
                ]
                genes = session.query(Gene).filter(Gene.id.in_(gene_ids_with_regions)).all()

            log.info("Scoring %d genes with %s...", len(genes), tool)

            for gene in genes:
                region_scores: dict[str, Optional[float]] = {}

                for region_type in _REGION_TYPES:
                    msa_path = build_msa(gene.id, region_type, session)
                    if msa_path is None:
                        continue
                    tmp_files.append(msa_path)

                    score = score_region(msa_path, tree_path, tool)
                    field = _REGION_TO_FIELD[region_type]
                    region_scores[field] = score

                if any(v is not None for v in region_scores.values()):
                    _upsert_phylo_score(session, gene.id, region_scores)
                    genes_scored += 1

            session.commit()

    finally:
        # Clean up temporary MSA and tree files
        for f in tmp_files:
            try:
                f.unlink(missing_ok=True)
            except Exception:
                pass
        if tree_path:
            try:
                tree_path.unlink(missing_ok=True)
            except Exception:
                pass

    log.info("Phylogenetic conservation scoring complete: %d genes scored", genes_scored)
    return genes_scored
