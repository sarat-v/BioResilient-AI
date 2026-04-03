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


def _s3_download(storage_root: str, relative_path: str, dst: Path) -> None:
    """Download *relative_path* from storage_root into *dst*.

    Uses boto3 for S3 (no aws-cli required; boto3 is always present).
    Falls back to a plain file copy for local storage roots.
    NOTE: Fusion only mounts *declared* Nextflow inputs as POSIX paths,
    so arbitrary s3:// URIs cannot be accessed via the filesystem —
    boto3 is the correct cross-container approach.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    if storage_root.startswith("s3://"):
        import boto3
        bucket = storage_root[len("s3://"):]          # e.g. "bioresilient-data"
        log.info("Downloading s3://%s/%s → %s", bucket, relative_path, dst)
        boto3.client("s3").download_file(bucket, relative_path, str(dst))
    elif storage_root.startswith("gs://"):
        from google.cloud import storage as gcs
        bucket_name, prefix = storage_root[len("gs://"):].split("/", 1) if "/" in storage_root[len("gs://"):] else (storage_root[len("gs://"):], "")
        blob_name = f"{prefix}/{relative_path}".lstrip("/")
        log.info("Downloading gs://%s/%s → %s", bucket_name, blob_name, dst)
        gcs.Client().bucket(bucket_name).blob(blob_name).download_to_filename(str(dst))
    else:
        import shutil
        src = Path(storage_root) / relative_path
        log.info("Copying %s → %s", src, dst)
        shutil.copy2(str(src), str(dst))


def _s3_upload(storage_root: str, relative_path: str, src: Path) -> None:
    """Upload *src* to *relative_path* under *storage_root*.

    Mirrors _s3_download — uses boto3 for S3, plain copy for local.
    """
    if storage_root.startswith("s3://"):
        import boto3
        bucket = storage_root[len("s3://"):]
        log.info("Uploading %s → s3://%s/%s", src, bucket, relative_path)
        boto3.client("s3").upload_file(str(src), bucket, relative_path)
    elif storage_root.startswith("gs://"):
        from google.cloud import storage as gcs
        bucket_name, prefix = storage_root[len("gs://"):].split("/", 1) if "/" in storage_root[len("gs://"):] else (storage_root[len("gs://"):], "")
        blob_name = f"{prefix}/{relative_path}".lstrip("/")
        gcs.Client().bucket(bucket_name).blob(blob_name).upload_from_filename(str(src))
    else:
        import shutil
        dst_path = Path(storage_root) / relative_path
        dst_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(src), str(dst_path))


def _s3_exists(storage_root: str, relative_path: str) -> bool:
    """Return True if *relative_path* exists under *storage_root*."""
    if storage_root.startswith("s3://"):
        import boto3
        from botocore.exceptions import ClientError
        bucket = storage_root[len("s3://"):]
        try:
            boto3.client("s3").head_object(Bucket=bucket, Key=relative_path)
            return True
        except ClientError:
            return False
    else:
        return (Path(storage_root) / relative_path).exists()


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


def run_step4c_chunk(args):
    gene_id_file = getattr(args, 'gene_id_file', None)
    if not gene_id_file:
        raise ValueError("--gene-id-file required for step4c_chunk")
    with open(gene_id_file) as f:
        gene_ids = [line.strip() for line in f if line.strip()]
    log.info("ESM-1v chunk: %d genes from %s", len(gene_ids), gene_id_file)
    from pipeline.layer1_sequence.esm1v import run_esm1v_pipeline
    n = run_esm1v_pipeline(gene_ids=gene_ids)
    log.info("ESM-1v chunk complete: %d motifs scored.", n)


def run_step4d(args):
    from pipeline.orchestrator import step4d_variant_direction
    step4d_variant_direction()


def run_step4d_scatter(args):
    """Emit N gene-batch files for parallel variant direction classification.

    Queries all gene_ids that have divergent motifs and splits them into
    batches of ~500 genes per file → typically 25–30 batch files for our
    12k-gene dataset.  Each file is named direction_batch_NNN.txt.
    """
    from db.models import Gene, Ortholog, DivergentMotif
    from db.session import get_session

    batch_size = 500  # genes per chunk; each chunk ~1–5 min on spot

    with get_session() as session:
        gene_ids = [
            row[0] for row in (
                session.query(Gene.id)
                .join(Ortholog, Ortholog.gene_id == Gene.id)
                .join(DivergentMotif, DivergentMotif.ortholog_id == Ortholog.id)
                .distinct()
            )
        ]

    log.info("step4d_scatter: %d genes → batches of %d", len(gene_ids), batch_size)

    for idx in range(0, len(gene_ids), batch_size):
        batch = gene_ids[idx : idx + batch_size]
        fname = f"direction_batch_{idx // batch_size:04d}.txt"
        with open(fname, "w") as fh:
            fh.write("\n".join(batch))
        log.info("  Wrote %s (%d genes)", fname, len(batch))

    log.info("step4d_scatter: %d batch files written.", (len(gene_ids) + batch_size - 1) // batch_size)


def run_step4d_chunk(args):
    """Classify variant directions for the gene subset in args.gene_id_file."""
    if not args.gene_id_file:
        raise ValueError("--gene-id-file is required for step4d_chunk")

    gene_ids = [line.strip() for line in open(args.gene_id_file) if line.strip()]
    log.info("step4d_chunk: %d genes from %s", len(gene_ids), args.gene_id_file)

    from pipeline.layer1_sequence.variant_direction import run_variant_direction_pipeline
    n = run_variant_direction_pipeline(gene_ids=gene_ids)
    log.info("step4d_chunk: %d motifs classified.", n)


def run_step5(args):
    from pipeline.orchestrator import step5_phylogenetic_tree
    import shutil as _shutil
    data = _load_aligned(args.input_pkl)
    aligned = data.get("aligned", data)
    treefile = step5_phylogenetic_tree(aligned)
    log.info("Tree written to %s", treefile)
    # Copy treefile into the NF work dir so Nextflow can stage it downstream.
    # treefile is a local path (e.g. /tmp/bioresilient/phylo/species.treefile)
    # — cannot _s3_download it; it was never uploaded by step5_phylogenetic_tree.
    dst = Path("species.treefile")
    if not dst.exists():
        _shutil.copy2(str(treefile), str(dst))
    log.info("species.treefile ready in work dir (%.1f KB)", dst.stat().st_size / 1e3)
    # Upload to S3 cache/ so future runs find it via the done-check.
    storage = os.environ.get("BIORESILIENT_STORAGE_ROOT", ".")
    try:
        _s3_upload(storage, "cache/species.treefile", dst)
        log.info("species.treefile uploaded to %s/cache/", storage)
    except Exception as exc:
        log.warning("Failed to upload species.treefile to storage: %s", exc)


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
        seen_sp = set()
        deduped_aln = {}
        for label, seq in codon_aln.items():
            sp = label.split("|")[0]
            if sp not in seen_sp:
                seen_sp.add(sp)
                deduped_aln[sp] = seq
        codon_aln = deduped_aln
        codon_species = list(codon_aln.keys())
        pruned_tree = prune_tree_to_species(treefile, codon_species)
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


def _find_result_jsons(results_dir: Path) -> list[Path]:
    """Walk results_dir to find result.json files, following symlinks (needed for Fusion).

    pathlib.glob('**') in Python <3.13 doesn't follow symlinks, which breaks
    on Fusion-staged directories.  os.walk(followlinks=True) handles them correctly.
    Falls back to subprocess `find` if os.walk also returns nothing.
    """
    import os
    import subprocess
    found = []
    for dirpath, _dirnames, filenames in os.walk(results_dir, followlinks=True):
        for fname in filenames:
            if fname == "result.json":
                found.append(Path(dirpath) / fname)
    if not found:
        log.warning("os.walk found 0 result.json in %s — trying `find` fallback", results_dir)
        try:
            proc = subprocess.run(
                ["find", str(results_dir), "-name", "result.json", "-type", "f"],
                capture_output=True, text=True, timeout=120,
            )
            found = [Path(p) for p in proc.stdout.strip().split("\n") if p]
        except Exception as exc:
            log.warning("find fallback failed: %s", exc)
    if not found:
        log.warning("pathlib.glob fallback for %s", results_dir)
        found = list(results_dir.glob("**/result.json"))
    log.info("Found %d result.json files in %s", len(found), results_dir)
    return found


def run_step6_collect(args):
    """Collect per-OG MEME results and write to DB. Supports both flat and batch-nested dirs."""
    from pipeline.layer2_evolution.selection import load_selection_scores, build_gene_og_map

    results_dir = Path(args.input_dir)
    gene_by_og = build_gene_og_map()
    all_results = _find_result_jsons(results_dir)
    batch: dict[str, dict] = {}
    for rfile in all_results:
        with open(rfile) as f:
            r = json.load(f)
        og_id = r.get("og_id") or rfile.parent.name
        if r and r.get("status") not in ("no_codon_alignment", "meme_failed", "skipped", "not_in_aligned"):
            batch[og_id] = r

    log.info("MEME collect: %d valid results out of %d total", len(batch), len(all_results))
    if batch:
        load_selection_scores(batch, gene_by_og)
        log.info("Loaded %d MEME results to DB", len(batch))
    else:
        log.warning("No valid MEME results found to load")


def _apply_fdr_correction(
    meme_batch: dict[str, dict],
    fb_batch: dict[str, dict],
) -> None:
    """Apply Benjamini-Hochberg FDR correction to MEME and BUSTED p-values.

    Corrects across all tested OGs (not per-batch) to properly control
    the false discovery rate.
    """
    try:
        import numpy as np

        meme_pvals = [(og, r.get("dnds_pvalue")) for og, r in meme_batch.items()
                      if r.get("selection_model") == "meme" and r.get("dnds_pvalue") is not None]
        if meme_pvals:
            ogs, pvals = zip(*meme_pvals)
            pvals_arr = np.array(pvals, dtype=float)
            n = len(pvals_arr)
            sorted_idx = np.argsort(pvals_arr)
            sorted_pvals = pvals_arr[sorted_idx]
            ranks = np.arange(1, n + 1)
            adjusted = np.minimum(1.0, sorted_pvals * n / ranks)
            for i in range(n - 2, -1, -1):
                adjusted[i] = min(adjusted[i], adjusted[i + 1])
            fdr = np.empty(n)
            fdr[sorted_idx] = adjusted
            for og, adj_p in zip(ogs, fdr):
                meme_batch[og]["meme_pvalue_fdr"] = round(float(adj_p), 6)
            log.info("BH-FDR corrected %d MEME p-values", n)

        busted_pvals = [(og, r.get("busted_pvalue")) for og, r in fb_batch.items()
                        if r.get("busted_pvalue") is not None]
        if busted_pvals:
            ogs, pvals = zip(*busted_pvals)
            pvals_arr = np.array(pvals, dtype=float)
            n = len(pvals_arr)
            sorted_idx = np.argsort(pvals_arr)
            sorted_pvals = pvals_arr[sorted_idx]
            ranks = np.arange(1, n + 1)
            adjusted = np.minimum(1.0, sorted_pvals * n / ranks)
            for i in range(n - 2, -1, -1):
                adjusted[i] = min(adjusted[i], adjusted[i + 1])
            fdr = np.empty(n)
            fdr[sorted_idx] = adjusted
            for og, adj_p in zip(ogs, fdr):
                fb_batch[og]["busted_pvalue_fdr"] = round(float(adj_p), 6)
            log.info("BH-FDR corrected %d BUSTED p-values", n)

    except Exception as exc:
        log.warning("FDR correction failed: %s", exc)


def run_step6_all_collect(args):
    """Collect all HyPhy results (BUSTED-PH + FEL + BUSTED + RELAX) from merged output."""
    from pipeline.layer2_evolution.selection import (
        load_selection_scores, load_fel_busted_scores, load_relax_scores, build_gene_og_map,
    )

    results_dir = Path(args.input_dir)
    gene_by_og = build_gene_og_map()
    all_results = _find_result_jsons(results_dir)

    meme_batch: dict[str, dict] = {}
    fb_batch: dict[str, dict] = {}
    relax_batch: dict[str, dict] = {}

    for rfile in all_results:
        with open(rfile) as f:
            r = json.load(f)
        og_id = r.get("og_id") or rfile.parent.name
        if r.get("status") in ("no_codon_alignment", "meme_failed", "skipped", "not_in_aligned"):
            continue

        # Map busted_ph_pvalue → dnds_pvalue so it scores via existing DB field.
        # selection_model="busted_ph" lets scoring.py distinguish from MEME results.
        if "busted_ph_pvalue" in r:
            r["dnds_pvalue"] = r["busted_ph_pvalue"]
            r["dnds_ratio"] = 0.0  # not applicable for BUSTED-PH

        meme_batch[og_id] = r

        if "fel_sites" in r or "busted_pvalue" in r:
            fb_batch[og_id] = {
                "og_id": og_id,
                "fel_sites": r.get("fel_sites", 0),
                "busted_pvalue": r.get("busted_pvalue", 1.0),
            }
        if "relax_k" in r or "relax_pvalue" in r:
            relax_batch[og_id] = {
                "og_id": og_id,
                "relax_k": r.get("relax_k"),
                "relax_pvalue": r.get("relax_pvalue"),
            }

    log.info("All-HyPhy collect: %d BUSTED-PH, %d FEL+BUSTED, %d RELAX results",
             len(meme_batch), len(fb_batch), len(relax_batch))

    # BH-FDR correction across all OGs before loading to DB
    _apply_fdr_correction(meme_batch, fb_batch)

    if meme_batch:
        load_selection_scores(meme_batch, gene_by_og)
        log.info("Loaded %d BUSTED-PH/proxy results to DB", len(meme_batch))
    if fb_batch:
        load_fel_busted_scores(fb_batch, gene_by_og)
        log.info("Loaded %d FEL+BUSTED results to DB", len(fb_batch))
    if relax_batch:
        load_relax_scores(relax_batch, gene_by_og)
        log.info("Loaded %d RELAX results to DB", len(relax_batch))


def run_step6_batch(args):
    """Batch HyPhy — codon alignment + BUSTED-PH / FEL / BUSTED / RELAX.

    Per-OG workflow:
      1. Build codon alignment + QC (requires ≥4 unique species, ≥60% CDS coverage)
      2. Run BUSTED-PH + FEL + BUSTED + RELAX in parallel (ThreadPoolExecutor)
    Falls back to protein-based aBSREL proxy if codon alignment unavailable.
    """
    import gc
    import json
    import shutil
    from concurrent.futures import ThreadPoolExecutor
    from pipeline.layer2_evolution.meme_selection import (
        load_cds_cache_pkl,
        fetch_cds_for_protein,
        protein_to_codon_alignment,
        run_busted_ph, parse_busted_ph_results,
        run_fel, parse_fel_results,
        run_busted, parse_busted_results,
        run_relax, parse_relax_results,
    )
    from pipeline.layer2_evolution.selection import write_hyphy_input, run_absrel, parse_absrel_results
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species
    from pipeline.config import get_local_storage_root

    og_id_file = getattr(args, 'og_id_file', None)
    if not og_id_file:
        raise ValueError("--og-id-file required for step6_batch")
    with open(og_id_file) as f:
        og_ids = [line.strip() for line in f if line.strip()]
    log.info("HyPhy batch: %d OGs from %s", len(og_ids), og_id_file)

    cds_cache_file = getattr(args, 'cds_cache', None)
    if cds_cache_file:
        log.info("Loading pre-fetched CDS cache from %s", cds_cache_file)
        load_cds_cache_pkl(cds_cache_file)
    else:
        log.warning("No CDS cache provided — batch tasks will skip NCBI fetching")

    data = _load_aligned(args.input_pkl)
    aligned = data.get("aligned", data)
    treefile = Path(args.treefile)
    out_dir = Path(args.output_dir or ".")
    out_dir.mkdir(parents=True, exist_ok=True)

    for og_id in og_ids:
        og_out = out_dir / og_id
        og_out.mkdir(parents=True, exist_ok=True)

        if og_id not in aligned:
            log.warning("OG %s not in aligned data, skipping", og_id)
            with open(og_out / "result.json", "w") as f:
                json.dump({"og_id": og_id, "status": "not_in_aligned"}, f)
            continue

        seqs = dict(aligned[og_id])
        species_ids = [label.split("|")[0] for label in seqs]

        cds_seqs: dict[str, str] = {}
        for label in seqs:
            parts = label.split("|")
            accession = parts[-1] if parts else label
            cds = fetch_cds_for_protein(accession)
            if cds:
                cds_seqs[label] = cds

        codon_aln = None
        cds_coverage = len(cds_seqs) / len(seqs) if seqs else 0
        log.info("OG %s: CDS coverage %d/%d (%.0f%%)", og_id, len(cds_seqs), len(seqs), cds_coverage * 100)

        if cds_coverage >= 0.60 and len(cds_seqs) >= 3:
            seqs_with_cds = {label: seq for label, seq in seqs.items() if label in cds_seqs}
            codon_aln = protein_to_codon_alignment(seqs_with_cds, cds_seqs)

        result = {"og_id": og_id, "cds_coverage": round(cds_coverage, 3)}
        if codon_aln:
            seen_sp = set()
            deduped_aln = {}
            for label, seq in codon_aln.items():
                sp = label.split("|")[0]
                if sp not in seen_sp:
                    seen_sp.add(sp)
                    deduped_aln[sp] = seq
            codon_aln = deduped_aln
            codon_species = list(codon_aln.keys())

            if len(codon_species) < 4:
                log.info("OG %s: only %d species after dedup, too few for HyPhy", og_id, len(codon_species))
                codon_aln = None

        if codon_aln:
            codon_species = list(codon_aln.keys())

            pruned_tree = prune_tree_to_species(treefile, codon_species)

            with ThreadPoolExecutor(max_workers=4) as executor:
                fut_busted_ph = executor.submit(run_busted_ph, codon_aln, pruned_tree, og_id)
                fut_fel = executor.submit(run_fel, codon_aln, pruned_tree, og_id)
                fut_busted = executor.submit(run_busted, codon_aln, pruned_tree, og_id)
                fut_relax = executor.submit(run_relax, codon_aln, pruned_tree, og_id)

                busted_ph_json = fut_busted_ph.result()
                fel_json = fut_fel.result()
                busted_json = fut_busted.result()
                relax_json = fut_relax.result()

            busted_ph_parsed = parse_busted_ph_results(busted_ph_json) if busted_ph_json else {"busted_ph_pvalue": 1.0}
            result["busted_ph_pvalue"] = busted_ph_parsed.get("busted_ph_pvalue", 1.0)
            result["selection_model"] = "busted_ph" if busted_ph_json is not None else "no_busted_ph"

            result["fel_sites"] = parse_fel_results(fel_json).get("fel_sites", 0) if fel_json else 0
            result["busted_pvalue"] = parse_busted_results(busted_json).get("busted_pvalue", 1.0) if busted_json else 1.0

            relax_parsed = parse_relax_results(relax_json) if relax_json else {}
            result["relax_k"] = relax_parsed.get("relax_k")
            result["relax_pvalue"] = relax_parsed.get("relax_pvalue", 1.0)
        else:
            pruned_tree = prune_tree_to_species(treefile, species_ids)
            aln_path, tree_path = write_hyphy_input(og_id, seqs, pruned_tree)
            raw = run_absrel(aln_path, tree_path, og_id)
            if raw is not None:
                result.update(parse_absrel_results(raw))
                result["selection_model"] = "proxy"
                # Filter branches to test (target-phenotype) species for scoring specificity
                from pipeline.layer2_evolution.meme_selection import _get_relax_branch_sets
                try:
                    test_set, _ref, _out = _get_relax_branch_sets()
                    if test_set and result.get("branches_under_selection"):
                        result["branches_under_selection"] = [
                            b for b in result["branches_under_selection"] if b in test_set
                        ]
                except Exception:
                    pass
            else:
                result["status"] = "no_codon_alignment"

        with open(og_out / "result.json", "w") as f:
            json.dump(result, f)

        gc.collect()

    log.info("HyPhy batch complete: %d OGs processed", len(og_ids))


def run_step6b_batch(args):
    """Batch FEL+BUSTED — reads codon alignments from MEME output dir, runs per-OG."""
    import json
    from Bio import SeqIO
    from pipeline.layer2_evolution.meme_selection import (
        run_fel, run_busted,
        parse_fel_results, parse_busted_results,
    )
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species

    meme_dir = Path(args.input_dir)
    treefile = Path(args.treefile)
    out_dir = Path(args.output_dir or ".")
    out_dir.mkdir(parents=True, exist_ok=True)

    og_dirs = [d for d in meme_dir.iterdir() if d.is_dir()]
    log.info("FEL+BUSTED batch: %d OG dirs in %s", len(og_dirs), meme_dir)

    for og_dir in og_dirs:
        og_id = og_dir.name
        og_out = out_dir / og_id
        og_out.mkdir(parents=True, exist_ok=True)

        codon_aln_path = og_dir / "codon_aln.fna"
        if not codon_aln_path.exists():
            log.info("OG %s: no codon alignment, skipping FEL+BUSTED", og_id)
            with open(og_out / "result.json", "w") as f:
                json.dump({"og_id": og_id, "status": "skipped"}, f)
            continue

        codon_aln = {r.id: str(r.seq) for r in SeqIO.parse(str(codon_aln_path), "fasta")}

        gene_treefile = og_dir / "gene.treefile"
        if gene_treefile.exists():
            tree = gene_treefile.read_text().strip()
        else:
            species_ids = [label.split("|")[0] for label in codon_aln]
            tree = prune_tree_to_species(treefile, species_ids)

        fel_json = run_fel(codon_aln, tree, og_id)
        busted_json = run_busted(codon_aln, tree, og_id)
        fel_result = parse_fel_results(fel_json) if fel_json else {"fel_sites": 0}
        busted_result = parse_busted_results(busted_json) if busted_json else {"busted_pvalue": 1.0}

        with open(og_out / "result.json", "w") as f:
            json.dump({"og_id": og_id, **fel_result, **busted_result}, f)

    log.info("FEL+BUSTED batch complete: %d OGs processed", len(og_dirs))


def run_step6c_batch(args):
    """Batch RELAX — reads codon alignments from MEME output dir, runs per-OG."""
    import json
    from Bio import SeqIO
    from pipeline.layer2_evolution.meme_selection import run_relax, parse_relax_results
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species

    meme_dir = Path(args.input_dir)
    treefile = Path(args.treefile)
    out_dir = Path(args.output_dir or ".")
    out_dir.mkdir(parents=True, exist_ok=True)

    og_dirs = [d for d in meme_dir.iterdir() if d.is_dir()]
    log.info("RELAX batch: %d OG dirs in %s", len(og_dirs), meme_dir)

    for og_dir in og_dirs:
        og_id = og_dir.name
        og_out = out_dir / og_id
        og_out.mkdir(parents=True, exist_ok=True)

        codon_aln_path = og_dir / "codon_aln.fna"
        if not codon_aln_path.exists():
            log.info("OG %s: no codon alignment, skipping RELAX", og_id)
            with open(og_out / "result.json", "w") as f:
                json.dump({"og_id": og_id, "status": "skipped"}, f)
            continue

        codon_aln = {r.id: str(r.seq) for r in SeqIO.parse(str(codon_aln_path), "fasta")}

        gene_treefile = og_dir / "gene.treefile"
        if gene_treefile.exists():
            tree = gene_treefile.read_text().strip()
        else:
            species_ids = [label.split("|")[0] for label in codon_aln]
            tree = prune_tree_to_species(treefile, species_ids)

        relax_json = run_relax(codon_aln, tree, og_id)
        relax_result = parse_relax_results(relax_json) if relax_json else {"relax_k": None, "relax_pvalue": 1.0}

        with open(og_out / "result.json", "w") as f:
            json.dump({"og_id": og_id, **relax_result}, f)

    log.info("RELAX batch complete: %d OGs processed", len(og_dirs))


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

    gene_treefile = codon_aln_path.parent / "gene.treefile"
    if gene_treefile.exists():
        tree = gene_treefile.read_text().strip()
    else:
        species_ids = [label.split("|")[0] for label in codon_aln]
        tree = prune_tree_to_species(treefile, species_ids)

    fel_json = run_fel(codon_aln, tree, og_id)
    busted_json = run_busted(codon_aln, tree, og_id)
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

    gene_treefile = codon_aln_path.parent / "gene.treefile"
    if gene_treefile.exists():
        tree = gene_treefile.read_text().strip()
    else:
        species_ids = [label.split("|")[0] for label in codon_aln]
        tree = prune_tree_to_species(treefile, species_ids)

    relax_json = run_relax(codon_aln, tree, og_id)
    relax_result = parse_relax_results(relax_json) if relax_json else {"relax_k": None, "relax_pvalue": 1.0}

    with open(out_dir / "result.json", "w") as f:
        json.dump({"og_id": og_id, **relax_result}, f)


def run_step6b_collect(args):
    """Collect per-OG FEL+BUSTED results and write to DB. Supports both flat and batch-nested dirs."""
    from pipeline.layer2_evolution.selection import load_fel_busted_scores, build_gene_og_map

    results_dir = Path(args.input_dir)
    gene_by_og = build_gene_og_map()
    batch: dict[str, dict] = {}
    for rfile in _find_result_jsons(results_dir):
        with open(rfile) as f:
            r = json.load(f)
        og_id = r.get("og_id") or rfile.parent.name
        if r and r.get("status") != "skipped":
            batch[og_id] = r

    log.info("FEL+BUSTED collect: %d valid results", len(batch))
    if batch:
        load_fel_busted_scores(batch, gene_by_og)
        log.info("Loaded %d FEL+BUSTED results to DB", len(batch))
    else:
        log.warning("No valid FEL+BUSTED results found to load")


def run_step6c_collect(args):
    """Collect per-OG RELAX results and write to DB. Supports both flat and batch-nested dirs."""
    from pipeline.layer2_evolution.selection import load_relax_scores, build_gene_og_map

    results_dir = Path(args.input_dir)
    gene_by_og = build_gene_og_map()
    batch: dict[str, dict] = {}
    for rfile in _find_result_jsons(results_dir):
        with open(rfile) as f:
            r = json.load(f)
        og_id = r.get("og_id") or rfile.parent.name
        if r and r.get("status") != "skipped":
            batch[og_id] = r

    log.info("RELAX collect: %d valid results", len(batch))
    if batch:
        load_relax_scores(batch, gene_by_og)
        log.info("Loaded %d RELAX results to DB", len(batch))
    else:
        log.warning("No valid RELAX results found to load")


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


def run_step6_prefetch_cds(args):
    """Pre-fetch all CDS sequences and output a single pickle for distribution."""
    from pipeline.layer2_evolution.meme_selection import export_cds_cache_pkl

    data = _load_aligned(args.input_pkl)
    aligned = data.get("aligned", data)
    out_path = getattr(args, 'output_dir', None) or "."
    pkl_out = Path(out_path) / "cds_cache.pkl"
    pkl_out.parent.mkdir(parents=True, exist_ok=True)
    export_cds_cache_pkl(aligned, str(pkl_out))
    log.info("CDS cache pickle: %s", pkl_out)


# ---------------------------------------------------------------------------
# Per-step "already done?" completion checks.
# Each function returns True if the step has sufficient data to skip.
# Used by --skip-if-done flag so a full step1→step9 run auto-skips completed steps.
# ---------------------------------------------------------------------------

def _db_count(sql: str, threshold: int = 1) -> bool:
    """Return True if the SQL scalar COUNT query returns ≥ threshold."""
    import psycopg2
    db_url = os.environ.get("DATABASE_URL", "")
    if not db_url:
        return False
    try:
        conn = psycopg2.connect(db_url, connect_timeout=10)
        cur = conn.cursor()
        cur.execute(sql)
        n = cur.fetchone()[0]
        conn.close()
        return (n or 0) >= threshold
    except Exception as exc:
        log.warning("_db_count check failed: %s", exc)
        return False


STEP_DONE_CHECKS: dict[str, callable] = {
    # step1: environment validated = DB reachable
    "step1": lambda: _db_count("SELECT 1", 1),

    # step2: ≥ 15 species downloaded and in Species table
    "step2": lambda: _db_count("SELECT COUNT(*) FROM species", 15),

    # step3: ≥ 50,000 orthologs (OrthoFinder done and loaded)
    "step3": lambda: _db_count("SELECT COUNT(*) FROM ortholog", 50_000),

    # step3b: same — orthologs loaded into DB
    "step3b": lambda: _db_count("SELECT COUNT(*) FROM ortholog", 50_000),

    # step3c: ≥ 10,000 nucleotide conservation scores in DB
    # (stored in gene.nucleotide_conservation_score or similar column)
    "step3c": lambda: _db_count(
        "SELECT COUNT(*) FROM gene WHERE nucleotide_conservation_score IS NOT NULL", 500
    ),

    # step4: ≥ 1,000,000 divergent motifs in DB (current run has 1.9M)
    "step4": lambda: _db_count("SELECT COUNT(*) FROM divergent_motif", 1_000_000),

    # step5: species tree built — check that treefile exists in storage via boto3
    # Stored at cache/species.treefile (same prefix as aligned_orthogroups.pkl)
    "step5": lambda: _s3_exists(
        os.environ.get("BIORESILIENT_STORAGE_ROOT", "."), "cache/species.treefile"
    ),

    # step3d: ≥ 500 genes with phylo_conservation_score set (step3d runs after step5)
    "step3d": lambda: _db_count(
        "SELECT COUNT(*) FROM gene WHERE phylo_conservation_score IS NOT NULL", 500
    ),

    # step7: ≥ 1,000 evolution_score rows with convergence_count > 0
    "step7": lambda: _db_count(
        "SELECT COUNT(*) FROM evolution_score WHERE convergence_count > 0", 1_000
    ),

    # step7b: ≥ 5,000 divergent_motif rows with convergent_aa_count > 0
    "step7b": lambda: _db_count(
        "SELECT COUNT(*) FROM divergent_motif WHERE convergent_aa_count > 0", 5_000
    ),

    # Steps intentionally NOT skipped — they need resume/re-run logic:
    # step4b: resume-safe but new motifs always need annotation  → do NOT skip
    # step4c: resume-safe (GPU) → do NOT skip; Python skips scored genes
    # step4d: re-run required (classifier was fixed)            → do NOT skip
    # step6: HyPhy scatter; OG-level skip via extract_og_ids    → do NOT skip
    # step8: GEO species list was fixed                         → do NOT skip
    # step8b: Bgee, best-effort                                 → do NOT skip
    # step9: always re-score after any upstream change          → do NOT skip
}


def _provide_step4_outputs(args) -> None:
    """Download aligned_orthogroups.pkl from S3 into the NF work dir when step4 is skipped."""
    dst = Path(args.output_dir or ".") / "aligned_orthogroups.pkl"
    if not dst.exists():
        storage = os.environ.get("BIORESILIENT_STORAGE_ROOT", ".")
        log.info("step4 skip: downloading aligned_orthogroups.pkl from %s/cache/", storage)
        _s3_download(storage, "cache/aligned_orthogroups.pkl", dst)
        log.info("aligned_orthogroups.pkl ready (%.1f MB)", dst.stat().st_size / 1e6)


def _provide_step5_outputs(args) -> None:
    """Download species.treefile from S3 into the NF work dir when step5 is skipped."""
    dst = Path("species.treefile")
    if not dst.exists():
        storage = os.environ.get("BIORESILIENT_STORAGE_ROOT", ".")
        log.info("step5 skip: downloading species.treefile from %s/cache/", storage)
        _s3_download(storage, "cache/species.treefile", dst)


# Registry of output-file providers called when a step is skipped.
# When --skip-if-done fires, the step function never runs, so Nextflow
# would not find declared output files. Each entry here copies the
# already-computed file from storage into the process work dir.
# To support a new step: add an entry here — no shell fallback needed in .nf files.
STEP_OUTPUT_PROVIDERS: dict[str, callable] = {
    "step4": _provide_step4_outputs,
    "step5": _provide_step5_outputs,
}


def check_done_and_exit(step: str, args=None) -> None:
    """If step already has sufficient data in RDS, provide NF outputs then exit 0."""
    check_fn = STEP_DONE_CHECKS.get(step)
    if check_fn is None:
        return  # no check registered → always run
    try:
        if check_fn():
            log.info("SKIP: step %s already complete — sufficient data found in DB.", step)
            provider = STEP_OUTPUT_PROVIDERS.get(step)
            if provider and args is not None:
                try:
                    provider(args)
                except Exception as exc:
                    log.warning("Output provider for %s failed (%s) — outputs may be missing.", step, exc)
            sys.exit(0)
        else:
            log.info("PROCEED: step %s data insufficient or absent — running step.", step)
    except SystemExit:
        raise
    except Exception as exc:
        log.warning("Done-check for step %s raised an error (%s) — proceeding.", step, exc)


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
    "step4c_chunk": run_step4c_chunk,
    "step4d": run_step4d,
    "step4d_scatter": run_step4d_scatter,
    "step4d_chunk": run_step4d_chunk,
    "step5": run_step5,
    "step6_single_og": run_step6_single_og,
    "step6_prefetch_cds": run_step6_prefetch_cds,
    "step6_batch": run_step6_batch,
    "step6_collect": run_step6_collect,
    "step6_all_collect": run_step6_all_collect,
    "step6b_single_og": run_step6b_single_og,
    "step6b_batch": run_step6b_batch,
    "step6b_collect": run_step6b_collect,
    "step6c_single_og": run_step6c_single_og,
    "step6c_batch": run_step6c_batch,
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
    parser.add_argument("--og-id-file", default=None)
    parser.add_argument("--input-pkl", default=None)
    parser.add_argument("--cds-cache", default=None)
    parser.add_argument("--input-dir", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--treefile", default=None)
    parser.add_argument("--codon-aln", default=None)
    parser.add_argument("--gene-id-file", default=None)
    parser.add_argument(
        "--skip-if-done",
        action="store_true",
        default=False,
        help=(
            "Exit 0 immediately if this step already has sufficient data in RDS. "
            "Enables seamless full-pipeline runs (step1→step9) that auto-skip completed steps."
        ),
    )

    args = parser.parse_args()
    _setup_env(args)

    if args.skip_if_done:
        check_done_and_exit(args.step, args)

    log.info("Running step: %s", args.step)
    STEP_MAP[args.step](args)
    log.info("Step %s complete", args.step)


if __name__ == "__main__":
    main()
