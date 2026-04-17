"""Step 6 — PAML codeml branch-site model A for positive selection analysis.

Statistical framework:
  Branch-site model A (Yang & Nielsen 2002 Genetics 161:1727-1737)
  H0 (null):        ω fixed at 1 in foreground branches (no positive selection)
  H1 (alternative): ω free  in foreground branches (positive selection allowed)
  LRT = 2 × (lnL_H1 − lnL_H0)

  The null hypothesis has a boundary constraint (ω=1 is on the boundary of
  the parameter space), so the LRT statistic follows a 50:50 mixture of
  χ²(df=1) and a point mass at zero.  Conservative approximation df=1 is used
  per Zhang et al. (2005) Mol Biol Evol 22:2472-2479 and Yang (2007) Mol
  Evolution p.130.  This is the standard for publication.

Output columns written to evolution_score:
  dnds_pvalue      ← LRT p-value (branch-site positive selection test)
  dnds_ratio       ← foreground ω in the positive selection site class (ω_fg)
  selection_model  ← "paml_branch_site"
  busted_pvalue    ← NULL  (HyPhy BUSTED not run; field intentionally empty)
  relax_k          ← NULL  (HyPhy RELAX not run; field intentionally empty)
  relax_pvalue     ← NULL  (HyPhy RELAX not run; field intentionally empty)
  fel_sites        ← NULL  (HyPhy FEL not run; field intentionally empty)

Rationale for PAML-only approach:
  FEL, BUSTED, and RELAX test different (and for this study less appropriate)
  hypotheses.  FEL detects pervasive selection across ALL branches — diluting
  the foreground signal.  BUSTED is gene-wide and non-directional.  RELAX
  tests constraint relaxation, not positive selection.  PAML branch-site model
  A directly and specifically tests the biological question: did positive
  selection act in the cancer-resistant (foreground) lineages?
"""

from __future__ import annotations

import json
import logging
import os
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Test-species cache (loaded once per worker process)
# ---------------------------------------------------------------------------

_test_species_cache: Optional[frozenset[str]] = None


def _get_paml_test_species() -> frozenset[str]:
    """Return cancer-resistant (foreground) species IDs from the DB.

    Cached for the lifetime of the process so the DB is only hit once per
    Batch job regardless of how many OGs are in the batch.
    """
    global _test_species_cache
    if _test_species_cache is not None:
        return _test_species_cache

    try:
        from db.models import Species
        from db.session import get_session

        with get_session() as session:
            rows = session.query(Species.id, Species.phenotypes, Species.is_control).all()

        _test_species_cache = frozenset(
            r.id
            for r in rows
            if (r.phenotypes or [])
            and "outgroup" not in (r.phenotypes or [])
            and not r.is_control
            and "baseline" not in (r.phenotypes or [])
        )
        log.info("PAML: loaded %d foreground (test) species", len(_test_species_cache))
    except Exception as exc:
        log.warning("PAML: could not load test species from DB (%s) — proceeding without labels", exc)
        _test_species_cache = frozenset()

    return _test_species_cache


# ---------------------------------------------------------------------------
# File-writing helpers
# ---------------------------------------------------------------------------


def _write_paml_phylip(codon_aln: dict[str, str], path: Path) -> None:
    """Write codon alignment in PHYLIP sequential format.

    PAML requires: first line is ' N  L', then each sequence on its own line
    as 'name  SEQUENCE' (name separated from sequence by ≥2 spaces).
    """
    seqs = list(codon_aln.items())
    if not seqs:
        raise ValueError("Empty codon alignment")

    n = len(seqs)
    aln_len = len(seqs[0][1])
    if aln_len % 3 != 0:
        raise ValueError(f"Alignment length {aln_len} is not a multiple of 3 (codon data required)")

    with open(path, "w") as f:
        f.write(f" {n}  {aln_len}\n")
        for name, seq in seqs:
            # Replace any gaps/ambiguous characters that might trip up PAML
            seq_clean = seq.upper().replace(".", "-")
            f.write(f"{name}  {seq_clean}\n")


def _label_tree_for_paml(tree_newick: str, test_species: frozenset[str]) -> str:
    """Add PAML foreground markers (#1) to test-species branches.

    PAML branch-site notation: 'species_name #1:branch_length'
    The marker goes between the name and the colon.

    Also strips internal node labels (bootstrap values) since PAML does not
    handle them in the branch-site model; it only needs leaf names.
    """
    # 1. Strip internal node numeric/string labels (keep only branch lengths)
    #    Pattern: ')something:' → '):' (internal nodes)
    labeled = re.sub(r'\)([\w.]+):', r'):', tree_newick)

    # 2. Mark each test species: 'name:' → 'name #1:'
    #    Use word-boundary match to avoid partial replacements.
    for sp in sorted(test_species, key=len, reverse=True):
        # Only replace if the species is present in the tree (exact word match)
        labeled = re.sub(
            r'\b(' + re.escape(sp) + r')(:)',
            r'\1 #1\2',
            labeled,
        )

    return labeled


def _write_paml_ctl(
    seq_file: str,
    tree_file: str,
    out_file: str,
    ctl_path: Path,
    fix_omega: int = 0,
    omega: float = 1.5,
) -> None:
    """Write a PAML codeml control file for branch-site model A.

    model=2 (branch models) + NSsites=2 (positive selection) = branch-site model A.

    H1 (alt):  fix_omega=0, omega=1.5  → free foreground ω, can exceed 1
    H0 (null): fix_omega=1, omega=1.0  → foreground ω fixed at 1 (no PS)
    """
    ctl = (
        f"      seqfile = {seq_file}\n"
        f"     treefile = {tree_file}\n"
        f"      outfile = {out_file}\n"
        "\n"
        "        noisy = 0\n"
        "      verbose = 0\n"
        "      runmode = 0\n"
        "      seqtype = 1\n"
        "    CodonFreq = 2\n"
        "        clock = 0\n"
        "        model = 2\n"
        "      NSsites = 2\n"
        "        icode = 0\n"
        "    fix_kappa = 0\n"
        "        kappa = 2\n"
        f"    fix_omega = {fix_omega}\n"
        f"        omega = {omega:.2f}\n"
        "    fix_alpha = 1\n"
        "        alpha = 0.\n"
        "       Malpha = 0\n"
        "        ncatG = 3\n"
        "        getSE = 0\n"
        " RateAncestor = 0\n"
        "   Small_Diff = 5e-7\n"
        "    cleandata = 0\n"
        "       method = 0\n"
    )
    ctl_path.write_text(ctl)


# ---------------------------------------------------------------------------
# Output parsers
# ---------------------------------------------------------------------------


def _parse_lnl(output_file: Path) -> Optional[float]:
    """Parse log-likelihood from codeml output file.

    Looks for the pattern: 'lnL(ntime: N  np: M):  -1234.567890'
    """
    try:
        text = output_file.read_text(errors="replace")
        m = re.search(r'lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*([-\d.]+)', text)
        if m:
            return float(m.group(1))
    except Exception as exc:
        log.debug("PAML lnL parse failed (%s): %s", output_file, exc)
    return None


def _parse_site_class_omegas(output_file: Path) -> tuple[Optional[float], Optional[float]]:
    """Parse foreground and background omega from branch-site model A output.

    Returns (omega_fg, omega_bg) where omega_fg is the ω in site class 2b
    (foreground positive selection class) and omega_bg is the background ω.

    The output section looks like:
        Parameters in branch-site model A, Model 2:
        p0 = ...  p1 = ...
        proportion   0.xxx   0.xxx   0.xxx   0.xxx
        background w  0.xxx   1.000   0.xxx   1.000
        foreground w  0.xxx   1.000   0.xxx  99.999    ← we want last value
    """
    try:
        text = output_file.read_text(errors="replace")
        fg_match = re.search(r'foreground w\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', text)
        bg_match = re.search(r'background w\s+([\d.]+)', text)

        omega_fg = float(fg_match.group(4)) if fg_match else None
        omega_bg = float(bg_match.group(1)) if bg_match else None

        # Cap extreme values (PAML reports 999 when synonymous substitutions ≈ 0)
        if omega_fg is not None:
            omega_fg = min(omega_fg, 99.0)

        return omega_fg, omega_bg
    except Exception as exc:
        log.debug("PAML omega parse failed (%s): %s", output_file, exc)
        return None, None


# ---------------------------------------------------------------------------
# Main PAML function
# ---------------------------------------------------------------------------


def run_paml_branch_site(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
    work_dir: Path,
) -> Optional[dict]:
    """Run PAML codeml branch-site model A (H0 vs H1 LRT).

    Args:
        codon_aln:    {species_id: codon_sequence}  (all seqs same length, mult of 3)
        tree_newick:  Newick string (pruned to species in codon_aln)
        og_id:        Orthogroup ID (for logging)
        work_dir:     Directory for PAML input/output files (will be created)

    Returns dict with keys:
        paml_pvalue, lrt_stat, lnl_h0, lnl_h1, omega_fg, omega_bg,
        test_species_count
    or None on failure.
    """
    from scipy.stats import chi2 as _chi2

    work_dir.mkdir(parents=True, exist_ok=True)

    test_species = _get_paml_test_species()
    present_test = [sp for sp in codon_aln if sp in test_species]

    if not present_test:
        log.debug("PAML %s: no test species in codon alignment — skipping", og_id)
        return None

    seq_file = "codon_aln.phy"
    tree_file = "tree.nwk"
    h1_out    = "h1_output.txt"
    h0_out    = "h0_output.txt"
    h1_ctl    = work_dir / "h1.ctl"
    h0_ctl    = work_dir / "h0.ctl"

    # Write alignment
    _write_paml_phylip(codon_aln, work_dir / seq_file)

    # Write tree with foreground markers
    labeled_tree = _label_tree_for_paml(tree_newick, frozenset(present_test))
    (work_dir / tree_file).write_text(labeled_tree)

    # Write control files
    _write_paml_ctl(seq_file, tree_file, h1_out, h1_ctl, fix_omega=0, omega=1.5)
    _write_paml_ctl(seq_file, tree_file, h0_out, h0_ctl, fix_omega=1, omega=1.0)

    paml_bin = "codeml"
    timeout = int(os.environ.get("PAML_TIMEOUT", "600"))

    # Run H1 (alternative: positive selection allowed in foreground)
    try:
        r1 = subprocess.run(
            [paml_bin, "h1.ctl"],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        if r1.returncode != 0:
            log.debug("PAML H1 rc=%d for %s: %s", r1.returncode, og_id, r1.stderr[:200])
    except subprocess.TimeoutExpired:
        log.warning("PAML H1 timed out for %s (>%ds)", og_id, timeout)
        return None
    except FileNotFoundError:
        log.error("codeml binary not found — is PAML installed?")
        return None

    lnl_h1 = _parse_lnl(work_dir / h1_out)
    if lnl_h1 is None:
        log.debug("PAML H1 produced no lnL for %s (stderr: %s)", og_id,
                  r1.stderr[:200] if r1 else "")
        return None

    # Run H0 (null: omega fixed at 1 in foreground)
    try:
        r0 = subprocess.run(
            [paml_bin, "h0.ctl"],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        if r0.returncode != 0:
            log.debug("PAML H0 rc=%d for %s: %s", r0.returncode, og_id, r0.stderr[:200])
    except subprocess.TimeoutExpired:
        log.warning("PAML H0 timed out for %s (>%ds)", og_id, timeout)
        return None

    lnl_h0 = _parse_lnl(work_dir / h0_out)
    if lnl_h0 is None:
        log.debug("PAML H0 produced no lnL for %s", og_id)
        return None

    # LRT statistic and p-value.
    # Branch-site model A has a boundary constraint on the null: under H0 the
    # foreground ω is fixed at 1 (boundary of the parameter space), so the LRT
    # statistic follows a 50:50 mixture of χ²(df=1) and a point mass at zero.
    # The standard conservative approximation (Zhang et al. 2005 Mol Biol Evol
    # 22:2472-2479; Yang 2007 Mol Evolution p.130) uses df=1, which is correct
    # and preferred for publication.  df=2 would be anti-conservative.
    lrt = max(0.0, 2.0 * (lnl_h1 - lnl_h0))
    pvalue = float(_chi2.sf(lrt, df=1)) if lrt > 0 else 1.0

    omega_fg, omega_bg = _parse_site_class_omegas(work_dir / h1_out)

    result = {
        "og_id": og_id,
        "paml_pvalue":       round(pvalue, 8),
        "lrt_stat":          round(lrt, 4),
        "lnl_h0":            round(lnl_h0, 4),
        "lnl_h1":            round(lnl_h1, 4),
        "omega_fg":          round(omega_fg, 4)  if omega_fg  is not None else None,
        "omega_bg":          round(omega_bg, 4)  if omega_bg  is not None else None,
        "test_species_count": len(present_test),
    }

    log.info(
        "PAML %s: p=%.4f  LRT=%.2f  ω_fg=%.2f  ω_bg=%.4f  test_sp=%d",
        og_id, pvalue, lrt,
        omega_fg or 0.0, omega_bg or 0.0, len(present_test),
    )
    return result


def parse_paml_results(paml_result: Optional[dict]) -> dict:
    """Convert PAML branch-site output to the EvolutionScore DB column format.

    Only PAML-native outputs are written.  HyPhy fields (busted_pvalue,
    relax_k, relax_pvalue, fel_sites) are left NULL because HyPhy tests are
    not run in this pipeline — leaving them NULL is scientifically honest and
    prevents downstream confusion if those columns are ever queried directly.
    """
    if paml_result is None:
        return {
            "dnds_pvalue":     1.0,
            "dnds_ratio":      None,
            "selection_model": "paml_no_signal",
            "busted_pvalue":   None,
            "relax_k":         None,
            "relax_pvalue":    None,
            "fel_sites":       None,
        }

    pvalue   = paml_result.get("paml_pvalue", 1.0)
    omega_fg = paml_result.get("omega_fg") or 0.0

    return {
        "dnds_pvalue":     round(pvalue, 8),
        "dnds_ratio":      round(omega_fg, 4),
        "selection_model": "paml_branch_site",
        "busted_pvalue":   None,   # HyPhy BUSTED not run — intentionally NULL
        "relax_k":         None,   # HyPhy RELAX not run — intentionally NULL
        "relax_pvalue":    None,   # HyPhy RELAX not run — intentionally NULL
        "fel_sites":       None,   # HyPhy FEL not run  — intentionally NULL
    }


# ---------------------------------------------------------------------------
# Parallel pipeline runner  (drop-in replacement for run_meme_pipeline)
# ---------------------------------------------------------------------------


def _paml_worker(
    args: tuple[str, dict[str, str], str, dict[str, list]],
) -> tuple[str, Optional[dict]]:
    """Process-pool worker: build codon alignment → run PAML → return parsed result.

    Args:
        args: (og_id, protein_alignment, species_tree_text, motifs_by_og)

    Returns:
        (og_id, parsed_result_dict)  where parsed_result_dict is None on failure.
    """
    og_id, protein_aln, tree_text, motifs_by_og = args

    try:
        from pipeline.layer2_evolution.meme_selection import (
            protein_to_codon_alignment,
            fetch_cds_for_protein,
        )
        from pipeline.config import get_local_storage_root

        # Build cds_seqs dict from the disk cache populated by _prefetch_all_cds
        # (called before the worker pool in run_paml_pipeline — no NCBI calls here).
        cds_seqs: dict[str, str] = {
            label: seq
            for label in protein_aln
            if (seq := fetch_cds_for_protein(label))
        }

        codon_aln = protein_to_codon_alignment(protein_aln, cds_seqs)
        if not codon_aln:
            # Not enough CDS data — use a null result (pval=1.0, no selection)
            return og_id, parse_paml_results(None)

        work_dir = Path(get_local_storage_root()) / "paml" / og_id
        result = run_paml_branch_site(codon_aln, tree_text, og_id, work_dir)
        return og_id, parse_paml_results(result)

    except Exception as exc:
        log.warning("PAML worker failed for %s: %s", og_id, exc)
        return og_id, None


def run_paml_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
    motifs_by_og: dict[str, list],
    species_treefile: Path,
    flush_callback=None,
) -> dict[str, dict]:
    """Run PAML branch-site model A in parallel for all candidate orthogroups.

    Drop-in replacement for run_meme_pipeline():
      - Same call signature
      - Same return type: {og_id: parsed_selection_result}
      - Writes selection_model="paml_branch_site" to every EvolutionScore row

    Args:
        aligned_orthogroups: {og_id: {species_id: protein_aligned_seq}}
        motifs_by_og:        {og_id: [motif_dict, ...]}  — divergence gate
        species_treefile:    Path to Newick species tree
        flush_callback:      optional callable(batch) for incremental DB writes
    """
    from pipeline.layer2_evolution.meme_selection import _prefetch_all_cds
    from pipeline.layer2_evolution.selection import _should_run_hyphy

    candidates = [og for og in aligned_orthogroups if _should_run_hyphy(og, motifs_by_og)]
    log.info(
        "PAML branch-site: running on %d candidate orthogroups...",
        len(candidates),
    )

    # Pre-fetch all CDS sequences to disk cache before spawning workers.
    # This avoids workers racing against NCBI's rate limit.
    _prefetch_all_cds(aligned_orthogroups)

    tree_text = Path(species_treefile).read_text().strip()
    n_cpu = os.cpu_count() or 4
    n_workers = max(1, n_cpu - 2)

    work_items = [
        (og_id, aligned_orthogroups[og_id], tree_text, motifs_by_og)
        for og_id in candidates
    ]

    all_results: dict[str, dict] = {}
    pending_flush: dict[str, dict] = {}
    success = 0
    no_cds = 0
    failed = 0
    done = 0
    total = len(candidates)
    _FLUSH_INTERVAL = 500

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_paml_worker, item): item[0] for item in work_items}
        for future in as_completed(futures):
            og_id, parsed = future.result()
            done += 1

            if parsed is not None:
                model = parsed.get("selection_model", "")
                if model == "paml_branch_site":
                    success += 1
                elif model == "paml_no_signal":
                    no_cds += 1
                else:
                    no_cds += 1
                all_results[og_id] = parsed
                pending_flush[og_id] = parsed
            else:
                failed += 1

            if done % 100 == 0 or done == total:
                log.info(
                    "  PAML: %d / %d (%.0f%%) — %d p-values, %d no-signal, %d failed",
                    done, total, 100 * done / total, success, no_cds, failed,
                )

            if flush_callback and len(pending_flush) >= _FLUSH_INTERVAL:
                flush_callback(dict(pending_flush))
                pending_flush.clear()

    if flush_callback and pending_flush:
        flush_callback(pending_flush)
        pending_flush.clear()

    log.info(
        "PAML complete: %d p-values, %d no-signal, %d failed",
        success, no_cds, failed,
    )
    return all_results
