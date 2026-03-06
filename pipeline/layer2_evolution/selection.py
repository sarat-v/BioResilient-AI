"""Step 6 — HyPhy aBSREL evolutionary selection analysis.

Runs HyPhy aBSREL on each candidate orthogroup that passes the divergence
threshold from Step 4. aBSREL (adaptive Branch-Site Random Effects Likelihood)
tests for episodic positive selection on each individual branch.

Threshold for running HyPhy: sequence identity < 85% in ≥ 2 protective species.
"""

import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from Bio import SeqIO

from db.models import EvolutionScore, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_storage_root, get_thresholds, get_tool_config

log = logging.getLogger(__name__)

HYPHY_SIGNIFICANCE_THRESHOLD = 0.05


def _should_run_hyphy(og_id: str, motifs_by_og: dict[str, list]) -> bool:
    """Gate: only run HyPhy on orthogroups with sufficient divergence signal."""
    thresholds = get_thresholds()
    min_species = thresholds.get("divergence_min_species", 2)
    species_with_motifs = {m["species_id"] for m in motifs_by_og.get(og_id, [])}
    # Exclude human from the count
    species_with_motifs.discard("human")
    return len(species_with_motifs) >= min_species


def write_hyphy_input(
    og_id: str,
    aligned_seqs: dict[str, str],
    tree_newick: str,
) -> tuple[Path, Path]:
    """Write alignment and tree files for HyPhy input.

    Returns (alignment_path, tree_path).
    """
    root = Path(get_storage_root())
    hyphy_dir = root / "hyphy" / og_id
    hyphy_dir.mkdir(parents=True, exist_ok=True)

    aln_path = hyphy_dir / "aligned.faa"
    with open(aln_path, "w") as f:
        for label, seq in aligned_seqs.items():
            f.write(f">{label}\n{seq}\n")

    tree_path = hyphy_dir / "species.treefile"
    tree_path.write_text(tree_newick)

    return aln_path, tree_path


def run_absrel(aln_path: Path, tree_path: Path, og_id: str) -> Optional[dict]:
    """Run HyPhy aBSREL and return parsed JSON results.

    Returns the parsed results dict, or None on failure.
    """
    root = Path(get_storage_root())
    out_path = root / "hyphy" / og_id / "results.json"

    if out_path.exists():
        try:
            with open(out_path) as f:
                return json.load(f)
        except json.JSONDecodeError:
            pass   # Re-run if corrupted

    cfg = get_tool_config()
    threads = cfg.get("hyphy_threads", 4)

    cmd = [
        "hyphy",
        "aBSREL",
        "--alignment", str(aln_path),
        "--tree", str(tree_path),
        "--output", str(out_path),
        "--branches", "All",
        "--srv", "Yes",
    ]

    log.info("  Running HyPhy aBSREL for %s...", og_id)
    result = subprocess.run(
        cmd,
        capture_output=True,
        check=False,
        env={**os.environ, "OMP_NUM_THREADS": str(threads)},
    )

    if result.returncode != 0:
        log.warning("  HyPhy failed for %s: %s", og_id, result.stderr.decode()[:300])
        return None

    if not out_path.exists():
        return None

    try:
        with open(out_path) as f:
            return json.load(f)
    except json.JSONDecodeError as exc:
        log.warning("  HyPhy output JSON invalid for %s: %s", og_id, exc)
        return None


def parse_absrel_results(hyphy_json: dict) -> dict:
    """Extract dN/dS, p-values, and selected branches from HyPhy aBSREL output.

    Returns:
      {
        "dnds_ratio": float,
        "dnds_pvalue": float,
        "selection_model": "aBSREL",
        "branches_under_selection": [species_label, ...],
      }
    """
    branch_results = hyphy_json.get("branch attributes", {}).get("0", {})

    min_pvalue = 1.0
    max_dnds = 0.0
    selected_branches = []

    for branch_name, attrs in branch_results.items():
        pvalue = attrs.get("Corrected P-value", attrs.get("Uncorrected P-value", 1.0))
        # aBSREL reports omega (dN/dS) per branch
        omega = attrs.get("Rate Distributions", [[1.0, 1.0]])[0][0]

        if pvalue is None:
            pvalue = 1.0

        if pvalue < HYPHY_SIGNIFICANCE_THRESHOLD:
            selected_branches.append(branch_name)
            if omega > max_dnds:
                max_dnds = omega
            if pvalue < min_pvalue:
                min_pvalue = pvalue

    # If no branch is significant, report the minimum observed p-value
    if not selected_branches:
        for attrs in branch_results.values():
            pv = attrs.get("Corrected P-value", 1.0) or 1.0
            if pv < min_pvalue:
                min_pvalue = pv
        max_dnds = 1.0   # neutral

    return {
        "dnds_ratio": round(max_dnds, 4),
        "dnds_pvalue": round(min_pvalue, 6),
        "selection_model": "aBSREL",
        "branches_under_selection": selected_branches,
    }


def load_selection_scores(
    selection_results: dict[str, dict],
    gene_by_og: dict[str, str],
) -> int:
    """Save EvolutionScore rows from HyPhy results.

    Args:
        selection_results: {og_id: parsed_absrel_result}
        gene_by_og: {og_id: gene_id}

    Returns:
        Number of rows inserted/updated.
    """
    saved = 0
    with get_session() as session:
        for og_id, result in selection_results.items():
            gene_id = gene_by_og.get(og_id)
            if not gene_id:
                continue

            ev = session.get(EvolutionScore, gene_id)
            if ev is None:
                ev = EvolutionScore(gene_id=gene_id)
                session.add(ev)

            ev.dnds_ratio = result["dnds_ratio"]
            ev.dnds_pvalue = result["dnds_pvalue"]
            ev.selection_model = result["selection_model"]
            ev.branches_under_selection = result["branches_under_selection"]
            saved += 1

    log.info("Saved evolution scores for %d genes.", saved)
    return saved


def build_gene_og_map() -> dict[str, str]:
    """Build {og_id: gene_id} from the Ortholog table."""
    og_map: dict[str, str] = {}
    with get_session() as session:
        for o in session.query(Ortholog).filter(Ortholog.orthofinder_og.isnot(None)):
            if o.orthofinder_og not in og_map:
                og_map[o.orthofinder_og] = o.gene_id
    return og_map


def run_selection_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
    motifs_by_og: dict[str, list],
    species_treefile: Path,
) -> dict[str, dict]:
    """Run HyPhy aBSREL for all candidate orthogroups.

    Args:
        aligned_orthogroups: {og_id: {label: aligned_seq}}
        motifs_by_og: {og_id: [motif_dict, ...]} — from divergence step
        species_treefile: Path to the IQ-TREE Newick file

    Returns:
        {og_id: parsed_absrel_result}
    """
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species

    all_results: dict[str, dict] = {}
    candidates = [og for og in aligned_orthogroups if _should_run_hyphy(og, motifs_by_og)]
    log.info("Running HyPhy aBSREL on %d candidate orthogroups...", len(candidates))

    for og_id in candidates:
        aligned_seqs = aligned_orthogroups[og_id]
        species_ids = [label.split("|")[0] for label in aligned_seqs]
        pruned_tree = prune_tree_to_species(species_treefile, species_ids)

        aln_path, tree_path = write_hyphy_input(og_id, aligned_seqs, pruned_tree)
        raw_result = run_absrel(aln_path, tree_path, og_id)

        if raw_result is None:
            continue

        parsed = parse_absrel_results(raw_result)
        all_results[og_id] = parsed
        log.debug("  %s: dN/dS=%.3f p=%.4f branches=%s",
                  og_id, parsed["dnds_ratio"], parsed["dnds_pvalue"],
                  parsed["branches_under_selection"])

    log.info("HyPhy complete: %d / %d orthogroups processed.", len(all_results), len(candidates))
    return all_results
