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

import hashlib
import logging
import os
import re
import shutil
import subprocess
import tempfile
import threading
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


def _build_nucleotide_concat_for_phylofit(treefile: Path) -> Optional[Path]:
    """Build a concatenated nucleotide CDS alignment for phyloFit training.

    phyloFit requires nucleotide data for nucleotide substitution models.
    We concatenate CDS regions from the NucleotideRegion table across genes
    that have coverage in ≥80% of species.
    """
    try:
        species_seqs: dict[str, list[str]] = {}
        with get_session() as session:
            genes_with_cds = (
                session.query(NucleotideRegion.gene_id)
                .filter_by(region_type="cds")
                .distinct()
                .all()
            )
            gene_ids = [g[0] for g in genes_with_cds]

            all_species = set()
            gene_data: list[dict[str, str]] = []
            for gene_id in gene_ids[:500]:
                regions = (
                    session.query(NucleotideRegion)
                    .filter_by(gene_id=gene_id, region_type="cds")
                    .all()
                )
                seqs = {r.species_id: r.sequence for r in regions if r.sequence and len(r.sequence) >= 30}
                if len(seqs) >= 3:
                    gene_data.append(seqs)
                    all_species.update(seqs.keys())

        if not gene_data or not all_species:
            log.warning("No suitable CDS regions for phyloFit training.")
            return None

        n_species = len(all_species)
        min_coverage = n_species * 0.5
        valid_genes = [g for g in gene_data if len(g) >= min_coverage]

        if len(valid_genes) < 5:
            log.warning("Only %d CDS genes with sufficient coverage for phyloFit.", len(valid_genes))
            valid_genes = gene_data[:50] if gene_data else []
            if not valid_genes:
                return None

        concat: dict[str, list[str]] = {sp: [] for sp in all_species}
        for gene_seqs in valid_genes[:200]:
            max_len = max(len(s) for s in gene_seqs.values())
            for sp in all_species:
                seq = gene_seqs.get(sp, "-" * max_len)
                concat[sp].append(seq.ljust(max_len, "-"))

        out_path = treefile.parent / "cds_concat_for_phylofit.fna"
        with open(out_path, "w") as f:
            for sp, parts in concat.items():
                f.write(f">{sp}\n{''.join(parts)}\n")

        log.info("Built nucleotide CDS concat for phyloFit: %d species, %d genes.",
                 len(concat), len(valid_genes))
        return out_path

    except Exception as exc:
        log.warning("Failed to build nucleotide concat for phyloFit: %s", exc)
        return None


def _fit_neutral_model(treefile: Path, concat_aln: Path) -> Optional[Path]:
    """Run phyloFit to produce a neutral .mod file for phyloP.

    phyloP requires a pre-fitted neutral evolutionary model (.mod file).
    We train phyloFit on a concatenated nucleotide CDS alignment using the
    REV (general reversible) nucleotide substitution model.

    Args:
        treefile: IQ-TREE2 species treefile (Newick with branch lengths).
        concat_aln: Concatenated protein alignment FASTA from step 5
                    (used only as fallback species list; nucleotide CDS
                    alignment is built from NucleotideRegion table).

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

    nt_aln = _build_nucleotide_concat_for_phylofit(treefile)
    if nt_aln is None:
        log.warning("Could not build nucleotide alignment for phyloFit — skipping.")
        return None

    try:
        log.info("Fitting neutral model with phyloFit on nucleotide CDS alignment...")
        result = subprocess.run(
            [
                "phyloFit",
                "--tree", str(nwk_path),
                "--msa-format", "FASTA",
                "--subst-mod", "REV",
                "--out-root", str(treefile.parent / "neutral"),
                str(nt_aln),
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
    """Parse phyloP --wig-scores output and return the mean conservation score.

    phyloP emits WIG format which contains header lines like:
      fixedStep chrom=gene1 start=1 step=1
    Taking line.split()[-1] on a score-only line is correct, but on a header
    line it yields tokens like "step=1" (ValueError, skipped) or large integers
    from start= coordinates.  We guard against misparsed header tokens by
    clamping to the physically meaningful phyloP range [-30, 30]; anything
    outside that range is a coordinate or header artefact, not a real score.
    """
    _PHYLOP_MIN, _PHYLOP_MAX = -30.0, 30.0
    scores = []
    for line in stdout.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        # Skip WIG header lines (fixedStep / variableStep / track)
        if line.startswith(("fixedStep", "variableStep", "track")):
            continue
        # WIG variableStep lines: "position score" — take the last token
        try:
            val = float(line.split()[-1])
        except ValueError:
            continue
        if _PHYLOP_MIN <= val <= _PHYLOP_MAX:
            scores.append(val)
        # else: silently discard — it's a misparsed coordinate or header artefact
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


# Cache of pruned .mod file content keyed by frozenset of species names.
# Most genes share the same species subset, so this eliminates ~99% of
# tree_doctor subprocess calls (38K calls → ~50 unique species sets).
_mod_cache: dict[frozenset, str] = {}
_mod_cache_lock = threading.Lock()


def _pruned_mod(mod_path: Path, species_in_msa: list[str]) -> Optional[Path]:
    """Write a pruned copy of the .mod file for the given species subset.

    Results are cached in memory by species set — tree_doctor is only called
    once per unique species combination across all 12K+ gene workers.
    """
    cache_key = frozenset(species_in_msa)

    # Check cache first (read with lock)
    with _mod_cache_lock:
        cached_content = _mod_cache.get(cache_key)

    if cached_content is None:
        # Cache miss — call tree_doctor once and store content
        tmp_path = _pruned_mod_via_tree_doctor(mod_path, sorted(species_in_msa))
        if tmp_path is None:
            # Fallback: ETE3 + manual TREE line replacement
            try:
                mod_text = mod_path.read_text()
                tree_match = re.search(r"^TREE:\s*(.+)$", mod_text, re.MULTILINE)
                if tree_match:
                    newick = tree_match.group(1).strip()
                    pruned_newick = _prune_newick_ete3(newick, set(species_in_msa))
                    if pruned_newick is not None:
                        cached_content = (mod_text[:tree_match.start(1)]
                                          + pruned_newick
                                          + mod_text[tree_match.end(1):])
            except Exception as exc:
                log.debug("_pruned_mod ETE3 fallback failed: %s", exc)
            if cached_content is None:
                log.debug("_pruned_mod: both tree_doctor and ETE3 failed for %d species",
                          len(species_in_msa))
                return None
        else:
            cached_content = tmp_path.read_text()
            tmp_path.unlink(missing_ok=True)

        with _mod_cache_lock:
            _mod_cache[cache_key] = cached_content

    # Write cached content to a fresh temp file for this call
    tmp = tempfile.NamedTemporaryFile(suffix=".mod", delete=False, mode="w")
    tmp.write(cached_content)
    tmp.close()
    return Path(tmp.name)


def run_phylop(msa_path: Path, mod_path: Path) -> Optional[float]:
    """Run phyloP on the MSA using a pre-pruned neutral model (.mod file).

    The mod_path is expected to already be pruned to the species in the MSA
    (done once upfront in run_phylo_conservation). No tree pruning here.
    """
    try:
        result = subprocess.run(
            ["phyloP", "--wig-scores", "--msa-format", "FASTA",
             str(mod_path), str(msa_path)],
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode == 0:
            return _parse_phylop_output(result.stdout)
        log.warning("phyloP non-zero exit (%d) for %s: %s",
                    result.returncode, msa_path.name, result.stderr[:200])
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

    Accepts pre-fetched sequences and pre-built pruned mod paths — no DB access,
    no tree_doctor calls, no temp .mod file writes in this function.

    args: (gene_id, {region_type: [(species_id, sequence)]},
                    {region_type: pruned_mod_path_str | None}, tool)
    """
    gene_id, region_seqs, region_mods, tool = args

    region_scores: dict[str, Optional[float]] = {}
    tmp_files: list[Path] = []

    try:
        for region_type in _REGION_TYPES:
            records = region_seqs.get(region_type, [])
            # Strip ambiguous IUPAC characters phyloP rejects (keep only ACGT-)
            clean_records = [
                (sp, "".join(c if c in "ACGTacgt-" else "N" for c in seq))
                for sp, seq in records
            ]
            msa_path = build_msa_from_records(gene_id, region_type, clean_records)
            if msa_path is None:
                continue
            tmp_files.append(msa_path)

            pruned_mod_str = region_mods.get(region_type)
            if pruned_mod_str is None:
                continue
            score = score_region(msa_path, Path(pruned_mod_str), tool)
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
    # Step 3d depends on the IQ-TREE species treefile produced by Step 5.
    # In the full pipeline the orchestrator guarantees Step 5 runs before Step 3d,
    # but when calling run_phylo_conservation() directly (e.g. in tests) the
    # treefile may be absent.  Fail fast with a clear message rather than silently
    # returning 0 and leaving all scores NULL.
    if not treefile.exists():
        raise FileNotFoundError(
            f"Species treefile not found at {treefile}.\n"
            "Step 3d (phylo_conservation) requires the IQ-TREE tree from Step 5. "
            "Run Step 5 (phylo_tree.run_iqtree) before Step 3d, or pass the "
            "correct treefile path."
        )

    tool = _detect_phylo_tool()
    if not tool:
        log.warning(
            "phyloP and phastCons not found. Phylogenetic conservation scores will be NULL. "
            "Install the PHAST package: conda install -c bioconda phast"
        )
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

        # Skip genes already cleanly scored — safe to resume mid-run.
        #
        # A gene is "already scored" when:
        #   (a) At least one phylo score column is non-NULL (cds OR promoter), AND
        #   (b) The cds score (if present) is not a WIG-header parsing artefact.
        #       Artefacts come from phyloP WIG "start=" coordinates being misread
        #       as scores; legitimate phyloP scores are in [-30, 30].
        #
        # Genes with only a promoter score (cds=NULL) are treated as done — they
        # may have had no CDS coverage in enough species.  Genes whose cds score
        # falls outside [-30, 30] will be rescored to overwrite the artefact.
        already_scored = {
            r[0] for r in session.query(PhyloConservationScore.gene_id)
            .filter(
                # Condition (a): at least one region scored
                (PhyloConservationScore.cds_phylo_score.isnot(None)) |
                (PhyloConservationScore.promoter_phylo_score.isnot(None))
            )
            .filter(
                # Condition (b): cds score absent (only promoter) OR in valid range
                PhyloConservationScore.cds_phylo_score.is_(None) |
                PhyloConservationScore.cds_phylo_score.between(-30.0, 30.0)
            )
            .all()
        }
        if already_scored:
            before = len(gene_ids)
            gene_ids = [g for g in gene_ids if g not in already_scored]
            log.info("Skipping %d already-scored genes; %d remaining.",
                     before - len(gene_ids), len(gene_ids))

        if not gene_ids:
            log.info("All genes already scored — step3d complete.")
            return len(already_scored)

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

    # Pre-build pruned .mod files for every unique species combination.
    # With 18 species and variable coverage, there are ~50-200 unique subsets.
    # Building them upfront (single-threaded, cached on disk) means workers
    # never call tree_doctor — they just pass a stable path to phyloP.
    log.info("Pre-building pruned neutral models for unique species subsets...")
    species_set_to_mod: dict[frozenset, Optional[Path]] = {}
    pruned_mod_dir = mod_path.parent / "pruned_mods"
    pruned_mod_dir.mkdir(exist_ok=True)

    for gene_id, region_map in gene_region_seqs.items():
        for region_type, records in region_map.items():
            if len(records) < _MIN_SPECIES_FOR_PHYLOP:
                continue
            key = frozenset(sp for sp, _ in records)
            if key not in species_set_to_mod:
                species_list = sorted(key)
                # Use a stable filename so we can skip re-building on reruns
                subset_hash = hashlib.md5(",".join(species_list).encode()).hexdigest()[:8]
                cached_path = pruned_mod_dir / f"pruned_{subset_hash}_{len(species_list)}sp.mod"
                if cached_path.exists():
                    species_set_to_mod[key] = cached_path
                else:
                    tmp_path = _pruned_mod_via_tree_doctor(mod_path, species_list)
                    if tmp_path is not None:
                        tmp_path.rename(cached_path)
                        species_set_to_mod[key] = cached_path
                    else:
                        species_set_to_mod[key] = None

    n_unique = sum(1 for v in species_set_to_mod.values() if v is not None)
    log.info("Pre-built %d pruned models for %d unique species subsets.",
             n_unique, len(species_set_to_mod))

    log.info("Phylogenetic conservation scoring with %s for %d genes (parallel)...",
             tool, len(gene_ids))

    # Build work items: pass stable pruned mod path per gene/region directly
    # so workers never call tree_doctor or write temp .mod files.
    # phyloP spawns external subprocesses that release the GIL — threads are
    # fine. Each gene runs 3 phyloP calls (cds/promoter/downstream), so with
    # N_cpu cores we can saturate all cores with N_cpu/3 * 2 concurrent genes.
    # Use all available CPUs; oversubscription is acceptable since phyloP
    # subprocesses are short-lived (~0.1–0.5s each).
    default_workers = max(1, os.cpu_count() or 8)
    n_workers = int(get_tool_config().get("phylop_workers", default_workers))

    def _make_work_item(gid: str) -> tuple:
        region_map = dict(gene_region_seqs.get(gid, {}))
        # Attach the pre-built mod path for each region's species set
        region_mods: dict[str, Optional[str]] = {}
        for rt, records in region_map.items():
            key = frozenset(sp for sp, _ in records)
            p = species_set_to_mod.get(key)
            region_mods[rt] = str(p) if p is not None else None
        return (gid, region_map, region_mods, tool)

    work_items = [_make_work_item(gid) for gid in gene_ids]
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
