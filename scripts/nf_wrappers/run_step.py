#!/usr/bin/env python3
"""
Thin CLI wrapper for running individual BioResilient pipeline steps.
Called by Nextflow processes inside Docker containers.

Usage:
  python -m scripts.nf_wrappers.run_step --step step2 --phenotype cancer_resistance
  python -m scripts.nf_wrappers.run_step --step step6 --phenotype cancer_resistance --og-id OG0001234
"""
import argparse
import json
import logging
import os
import pickle
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-5s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("nf_wrapper")


def _setup_env(args):
    if args.db_url:
        os.environ["DATABASE_URL"] = args.db_url
    if args.storage_root:
        os.environ["BIORESILIENT_STORAGE_ROOT"] = args.storage_root
    if args.ncbi_api_key:
        os.environ["NCBI_API_KEY"] = args.ncbi_api_key


def _load_species(phenotype: str) -> list[dict]:
    from pipeline.config import load_species_registry
    return load_species_registry(phenotype)


def _load_aligned(path: str) -> dict:
    with open(path, "rb") as f:
        return pickle.load(f)


def run_step1(args):
    from pipeline.orchestrator import step1_validate_environment
    step1_validate_environment()


def run_step2(args):
    from pipeline.orchestrator import step2_download_proteomes
    species = _load_species(args.phenotype)
    step2_download_proteomes(species)


def run_step3(args):
    from pipeline.orchestrator import step3_run_orthofinder
    step3_run_orthofinder()


def run_step3b(args):
    from pipeline.orchestrator import step3b_load_orthologs
    results_dir = Path(args.input_dir) if args.input_dir else None
    step3b_load_orthologs(results_dir)


def run_step3c(args):
    from pipeline.orchestrator import step3c_nucleotide_conservation
    species = _load_species(args.phenotype)
    step3c_nucleotide_conservation(species)


def run_step3d(args):
    from pipeline.orchestrator import step3d_phylo_conservation
    step3d_phylo_conservation()


def run_step4(args):
    from pipeline.orchestrator import step4_alignment_and_divergence
    species = _load_species(args.phenotype)
    aligned, motifs = step4_alignment_and_divergence(species)
    out = Path(args.output_dir or ".")
    with open(out / "aligned_orthogroups.pkl", "wb") as f:
        pickle.dump({"aligned": aligned, "motifs_by_og": motifs}, f)
    log.info("Wrote aligned_orthogroups.pkl (%d OGs)", len(aligned))


def run_step4b(args):
    from pipeline.orchestrator import step4b_domain_and_consequence
    step4b_domain_and_consequence()


def run_step4c(args):
    from pipeline.orchestrator import step4c_esm1v
    step4c_esm1v()


def run_step4d(args):
    from pipeline.orchestrator import step4d_variant_direction
    step4d_variant_direction()


def run_step5(args):
    from pipeline.orchestrator import step5_phylogenetic_tree
    data = _load_aligned(args.input_pkl)
    aligned = data.get("aligned", data)
    treefile = step5_phylogenetic_tree(aligned)
    log.info("Tree written to %s", treefile)


def run_step6_single_og(args):
    """Per-OG MEME — called by Nextflow scatter."""
    from pipeline.layer2_evolution.meme_selection import (
        _build_codon_alignment,
        _run_hyphy_meme,
        _prefetch_cds_for_og,
    )
    og_id = args.og_id
    data = _load_aligned(args.input_pkl)
    aligned = data.get("aligned", data)

    if og_id not in aligned:
        log.warning("OG %s not in aligned data, skipping", og_id)
        return

    treefile = Path(args.treefile)
    out_dir = Path(args.output_dir or ".")
    out_dir.mkdir(parents=True, exist_ok=True)

    seqs = aligned[og_id]
    _prefetch_cds_for_og(seqs)
    codon_aln = _build_codon_alignment(og_id, seqs, out_dir)
    if codon_aln:
        result = _run_hyphy_meme(og_id, codon_aln, treefile, out_dir)
        with open(out_dir / "result.json", "w") as f:
            json.dump(result or {}, f)
    else:
        log.warning("No codon alignment for %s", og_id)
        with open(out_dir / "result.json", "w") as f:
            json.dump({"og_id": og_id, "status": "no_codon_alignment"}, f)


def run_step6_collect(args):
    """Collect per-OG MEME results and write to DB."""
    from pipeline.layer2_evolution.selection import load_selection_scores
    from pipeline.layer2_evolution.meme_selection import build_gene_og_map

    results_dir = Path(args.input_dir)
    gene_by_og = build_gene_og_map()
    batch = []
    for rfile in results_dir.glob("*/result.json"):
        with open(rfile) as f:
            r = json.load(f)
        if r and r.get("status") != "no_codon_alignment":
            batch.append(r)

    if batch:
        load_selection_scores(batch, gene_by_og)
        log.info("Loaded %d MEME results to DB", len(batch))


def run_step6b_single_og(args):
    """Per-OG FEL+BUSTED — called by Nextflow scatter."""
    og_id = args.og_id
    out_dir = Path(args.output_dir or ".")
    out_dir.mkdir(parents=True, exist_ok=True)
    codon_aln = Path(args.codon_aln)
    treefile = Path(args.treefile)

    if not codon_aln.exists():
        log.warning("No codon alignment for %s", og_id)
        json.dump({"og_id": og_id, "status": "skipped"}, open(out_dir / "result.json", "w"))
        return

    from pipeline.layer2_evolution.meme_selection import _run_hyphy_tool
    fel = _run_hyphy_tool("fel", og_id, codon_aln, treefile, out_dir)
    busted = _run_hyphy_tool("busted", og_id, codon_aln, treefile, out_dir)
    with open(out_dir / "result.json", "w") as f:
        json.dump({"og_id": og_id, "fel": fel, "busted": busted}, f)


def run_step6c_single_og(args):
    """Per-OG RELAX — called by Nextflow scatter."""
    og_id = args.og_id
    out_dir = Path(args.output_dir or ".")
    out_dir.mkdir(parents=True, exist_ok=True)
    codon_aln = Path(args.codon_aln)
    treefile = Path(args.treefile)

    if not codon_aln.exists():
        log.warning("No codon alignment for %s", og_id)
        json.dump({"og_id": og_id, "status": "skipped"}, open(out_dir / "result.json", "w"))
        return

    from pipeline.layer2_evolution.meme_selection import _run_hyphy_tool
    relax = _run_hyphy_tool("relax", og_id, codon_aln, treefile, out_dir)
    with open(out_dir / "result.json", "w") as f:
        json.dump({"og_id": og_id, "relax": relax}, f)


def run_step7(args):
    from pipeline.orchestrator import step7_convergence
    step7_convergence()


def run_step7b(args):
    from pipeline.orchestrator import step7b_convergent_aa
    step7b_convergent_aa()


def run_step8(args):
    from pipeline.orchestrator import step8_expression
    species = _load_species(args.phenotype)
    step8_expression(species)


def run_step8b(args):
    from pipeline.orchestrator import step8b_bgee
    step8b_bgee()


def run_step9(args):
    from pipeline.orchestrator import step9_composite_score
    step9_composite_score()


def run_step10b(args):
    from pipeline.orchestrator import step10b_alphagenome
    step10b_alphagenome()


def run_step11(args):
    from pipeline.orchestrator import step11_disease_annotation
    step11_disease_annotation(trait_id=args.phenotype)


def run_step11b(args):
    from pipeline.orchestrator import step11b_rare_variants
    step11b_rare_variants(trait_id=args.phenotype)


def run_step11c(args):
    from pipeline.orchestrator import step11c_literature
    step11c_literature(trait_id=args.phenotype)


def run_step11d(args):
    from pipeline.orchestrator import step11d_pathway_convergence
    step11d_pathway_convergence()


def run_step12(args):
    from pipeline.orchestrator import step12_druggability
    step12_druggability(trait_id=args.phenotype)


def run_step12b(args):
    from pipeline.orchestrator import step12b_p2rank
    step12b_p2rank(trait_id=args.phenotype)


def run_step13(args):
    from pipeline.orchestrator import step13_gene_therapy
    step13_gene_therapy(trait_id=args.phenotype)


def run_step14(args):
    from pipeline.orchestrator import step14_safety
    step14_safety(trait_id=args.phenotype)


def run_step14b(args):
    from pipeline.orchestrator import step14b_depmap_gtex
    step14b_depmap_gtex(trait_id=args.phenotype)


def run_step15(args):
    from pipeline.orchestrator import step15_rescore
    step15_rescore()


STEP_MAP = {
    "step1": run_step1,
    "step2": run_step2,
    "step3": run_step3,
    "step3b": run_step3b,
    "step3c": run_step3c,
    "step3d": run_step3d,
    "step4": run_step4,
    "step4b": run_step4b,
    "step4c": run_step4c,
    "step4d": run_step4d,
    "step5": run_step5,
    "step6_single_og": run_step6_single_og,
    "step6_collect": run_step6_collect,
    "step6b_single_og": run_step6b_single_og,
    "step6c_single_og": run_step6c_single_og,
    "step7": run_step7,
    "step7b": run_step7b,
    "step8": run_step8,
    "step8b": run_step8b,
    "step9": run_step9,
    "step10b": run_step10b,
    "step11": run_step11,
    "step11b": run_step11b,
    "step11c": run_step11c,
    "step11d": run_step11d,
    "step12": run_step12,
    "step12b": run_step12b,
    "step13": run_step13,
    "step14": run_step14,
    "step14b": run_step14b,
    "step15": run_step15,
}


def main():
    parser = argparse.ArgumentParser(description="BioResilient Nextflow step wrapper")
    parser.add_argument("--step", required=True, choices=sorted(STEP_MAP.keys()))
    parser.add_argument("--phenotype", default="cancer_resistance")
    parser.add_argument("--db-url", default=None)
    parser.add_argument("--storage-root", default=None)
    parser.add_argument("--ncbi-api-key", default=None)
    parser.add_argument("--og-id", default=None)
    parser.add_argument("--input-pkl", default=None)
    parser.add_argument("--input-dir", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--treefile", default=None)
    parser.add_argument("--codon-aln", default=None)

    args = parser.parse_args()
    _setup_env(args)

    log.info("Running step: %s", args.step)
    STEP_MAP[args.step](args)
    log.info("Step %s complete", args.step)


if __name__ == "__main__":
    main()
