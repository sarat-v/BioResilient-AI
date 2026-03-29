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
    from pipeline.orchestrator import _load_species_registry
    return _load_species_registry(phenotype)


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
    """Per-OG MEME — called by Nextflow scatter. Mirrors _meme_worker logic."""
    import gc
    from pipeline.layer2_evolution.meme_selection import (
        fetch_cds_for_protein,
        protein_to_codon_alignment,
        run_meme,
        parse_meme_results,
    )
    from pipeline.layer2_evolution.selection import (
        write_hyphy_input, run_absrel, parse_absrel_results,
    )
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species

    og_id = args.og_id
    data = _load_aligned(args.input_pkl)
    aligned = data.get("aligned", data)

    if og_id not in aligned:
        log.warning("OG %s not in aligned data, skipping", og_id)
        return

    treefile = Path(args.treefile)
    out_dir = Path(args.output_dir or ".")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Extract only this OG's data then immediately free the full pickle from
    # memory — the pkl holds all OGs and can be several GB. HyPhy needs the
    # headroom, especially for large OGs like OG0000000.
    seqs = dict(aligned[og_id])
    del data, aligned
    gc.collect()

    species_ids = [label.split("|")[0] for label in seqs]
    pruned_tree = prune_tree_to_species(treefile, species_ids)

    cds_seqs: dict[str, str] = {}
    for label in seqs:
        parts = label.split("|")
        accession = parts[-1] if parts else label
        cds = fetch_cds_for_protein(accession)
        if cds:
            cds_seqs[label] = cds

    # Run MEME if we have CDS for ≥60% of species in the OG (strict 100% caused
    # all OGs to fall through to aBSREL proxy when any one NCBI fetch failed).
    codon_aln = None
    cds_coverage = len(cds_seqs) / len(seqs) if seqs else 0
    log.info("OG %s: CDS coverage %d/%d (%.0f%%)", og_id, len(cds_seqs), len(seqs), cds_coverage * 100)
    if cds_coverage >= 0.60 and len(cds_seqs) >= 3:
        # Subset seqs and pruned tree to only those species with CDS
        seqs_with_cds = {label: seq for label, seq in seqs.items() if label in cds_seqs}
        cds_species = [label.split("|")[0] for label in cds_seqs]
        pruned_tree = prune_tree_to_species(treefile, cds_species)
        codon_aln = protein_to_codon_alignment(seqs_with_cds, cds_seqs)

    result = {"og_id": og_id, "cds_coverage": round(cds_coverage, 3)}
    if codon_aln is not None:
        meme_json = run_meme(codon_aln, pruned_tree, og_id)
        if meme_json is not None:
            result.update(parse_meme_results(meme_json, og_id))
            result["selection_model"] = "meme"
            # Copy meme outputs from local storage root into out_dir so Nextflow
            # stages them to S3 — run_meme writes to get_local_storage_root()/meme/og_id/
            # which is NOT the task work dir that Nextflow syncs.
            from pipeline.config import get_local_storage_root
            import shutil
            meme_local = Path(get_local_storage_root()) / "meme" / og_id
            for fname in ("codon_aln.fna", "meme.json", "species.treefile"):
                src = meme_local / fname
                if src.exists():
                    shutil.copy2(src, out_dir / fname)
        else:
            result["status"] = "meme_failed"
    else:
        aln_path, tree_path = write_hyphy_input(og_id, seqs, pruned_tree)
        raw = run_absrel(aln_path, tree_path, og_id)
        if raw is not None:
            result.update(parse_absrel_results(raw))
            result["selection_model"] = "proxy"
        else:
            result["status"] = "no_codon_alignment"

    with open(out_dir / "result.json", "w") as f:
        json.dump(result, f)


def run_step6_collect(args):
    """Collect per-OG MEME results and write to DB."""
    from pipeline.layer2_evolution.selection import load_selection_scores, build_gene_og_map

    results_dir = Path(args.input_dir)
    gene_by_og = build_gene_og_map()
    batch: dict[str, dict] = {}
    for rfile in results_dir.glob("*/result.json"):
        with open(rfile) as f:
            r = json.load(f)
        og_id = r.get("og_id") or rfile.parent.name
        if r and r.get("status") not in ("no_codon_alignment", "meme_failed", "skipped"):
            batch[og_id] = r

    if batch:
        load_selection_scores(batch, gene_by_og)
        log.info("Loaded %d MEME results to DB", len(batch))


def run_step6b_single_og(args):
    """Per-OG FEL+BUSTED — called by Nextflow scatter."""
    from Bio import SeqIO
    from pipeline.layer2_evolution.meme_selection import (
        run_fel, run_busted,
        parse_fel_results, parse_busted_results,
    )
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species

    og_id = args.og_id
    out_dir = Path(args.output_dir or ".")
    out_dir.mkdir(parents=True, exist_ok=True)
    codon_aln_path = Path(args.codon_aln)
    treefile = Path(args.treefile)

    if not codon_aln_path.exists():
        log.warning("No codon alignment for %s", og_id)
        with open(out_dir / "result.json", "w") as f:
            json.dump({"og_id": og_id, "status": "skipped"}, f)
        return

    codon_aln = {r.id: str(r.seq) for r in SeqIO.parse(str(codon_aln_path), "fasta")}
    species_ids = [label.split("|")[0] for label in codon_aln]
    pruned_tree = prune_tree_to_species(treefile, species_ids)

    fel_json = run_fel(codon_aln, pruned_tree, og_id)
    busted_json = run_busted(codon_aln, pruned_tree, og_id)
    fel_result = parse_fel_results(fel_json) if fel_json else {"fel_sites": 0}
    busted_result = parse_busted_results(busted_json) if busted_json else {"busted_pvalue": 1.0}

    with open(out_dir / "result.json", "w") as f:
        json.dump({"og_id": og_id, **fel_result, **busted_result}, f)


def run_step6c_single_og(args):
    """Per-OG RELAX — called by Nextflow scatter."""
    from Bio import SeqIO
    from pipeline.layer2_evolution.meme_selection import run_relax, parse_relax_results
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species

    og_id = args.og_id
    out_dir = Path(args.output_dir or ".")
    out_dir.mkdir(parents=True, exist_ok=True)
    codon_aln_path = Path(args.codon_aln)
    treefile = Path(args.treefile)

    if not codon_aln_path.exists():
        log.warning("No codon alignment for %s", og_id)
        with open(out_dir / "result.json", "w") as f:
            json.dump({"og_id": og_id, "status": "skipped"}, f)
        return

    codon_aln = {r.id: str(r.seq) for r in SeqIO.parse(str(codon_aln_path), "fasta")}
    species_ids = [label.split("|")[0] for label in codon_aln]
    pruned_tree = prune_tree_to_species(treefile, species_ids)

    relax_json = run_relax(codon_aln, pruned_tree, og_id)
    relax_result = parse_relax_results(relax_json) if relax_json else {"relax_k": None, "relax_pvalue": 1.0}

    with open(out_dir / "result.json", "w") as f:
        json.dump({"og_id": og_id, **relax_result}, f)


def run_step6b_collect(args):
    """Collect per-OG FEL+BUSTED results and write to DB."""
    from pipeline.layer2_evolution.selection import load_fel_busted_scores, build_gene_og_map

    results_dir = Path(args.input_dir)
    gene_by_og = build_gene_og_map()
    batch: dict[str, dict] = {}
    for rfile in results_dir.glob("*/result.json"):
        with open(rfile) as f:
            r = json.load(f)
        og_id = r.get("og_id") or rfile.parent.name
        if r and r.get("status") != "skipped":
            batch[og_id] = r

    if batch:
        load_fel_busted_scores(batch, gene_by_og)
        log.info("Loaded %d FEL+BUSTED results to DB", len(batch))


def run_step6c_collect(args):
    """Collect per-OG RELAX results and write to DB."""
    from pipeline.layer2_evolution.selection import load_relax_scores, build_gene_og_map

    results_dir = Path(args.input_dir)
    gene_by_og = build_gene_og_map()
    batch: dict[str, dict] = {}
    for rfile in results_dir.glob("*/result.json"):
        with open(rfile) as f:
            r = json.load(f)
        og_id = r.get("og_id") or rfile.parent.name
        if r and r.get("status") != "skipped":
            batch[og_id] = r

    if batch:
        load_relax_scores(batch, gene_by_og)
        log.info("Loaded %d RELAX results to DB", len(batch))


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


def run_dedup_motifs(args):
    """Remove duplicate divergent_motif rows created by running step4 more than once.
    Keeps the row with motif_direction populated (gnomAD-aware); falls back to latest id."""
    from db.session import get_session
    from sqlalchemy import text

    with get_session() as s:
        before = s.execute(text("SELECT COUNT(*) FROM divergent_motif")).scalar()
        s.execute(text("""
            DELETE FROM divergent_motif
            WHERE id IN (
                SELECT id FROM (
                    SELECT id,
                           ROW_NUMBER() OVER (
                               PARTITION BY ortholog_id, start_pos
                               ORDER BY motif_direction NULLS LAST, id DESC
                           ) AS rn
                    FROM divergent_motif
                ) ranked
                WHERE rn > 1
            )
        """))
        after = s.execute(text("SELECT COUNT(*) FROM divergent_motif")).scalar()
    log.info("Deduped divergent_motif: %d → %d rows (removed %d duplicates)", before, after, before - after)


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
    "step6b_collect": run_step6b_collect,
    "step6c_single_og": run_step6c_single_og,
    "step6c_collect": run_step6c_collect,
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
    "dedup_motifs": run_dedup_motifs,
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
