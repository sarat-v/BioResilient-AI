"""Phase 1 Pipeline Orchestrator.

Runs all steps in the correct dependency order. Each step is gated:
the pipeline will not advance past a step that fails.

Usage:
    python pipeline/orchestrator.py
    python pipeline/orchestrator.py --resume-from step5
    python pipeline/orchestrator.py --steps step1,step2
    python pipeline/orchestrator.py --phenotype cancer_resistance
    python pipeline/orchestrator.py --dry-run
"""

import argparse
import json
import logging
import os
import pickle
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Pipeline state file — read by the API for the frontend status panel
# ---------------------------------------------------------------------------

_REPO_ROOT = Path(__file__).resolve().parents[1]
_STATE_FILE = _REPO_ROOT / "pipeline_state.json"
_LOG_FILE = _REPO_ROOT / "pipeline.log"

STEP_LABELS = {
    "step1":   "Validate environment",
    "step2":   "Download proteomes",
    "step3":   "Run OrthoFinder",
    "step3b":  "Load orthologs into DB",
    "step3c":  "Nucleotide region extraction + alignment",
    "step3d":  "Phylogenetic conservation scoring (phyloP/PhastCons) — requires tree from step5",
    "step4":   "Align sequences & find divergent motifs",
    "step4b":  "Domain annotation (Pfam) + functional consequence (AlphaMissense)",
    "step4c":  "ESM-1v structural variant effect scoring",
    "step4d":  "Variant direction inference (GoF / LoF / neutral)",
    "step5":   "Build phylogenetic tree",
    "step6":   "Evolutionary selection (HyPhy MEME)",
    "step6b":  "FEL + BUSTED supplementary selection tests",
    "step6c":  "RELAX branch-specific rate acceleration",
    "step7":   "Convergence & conservation scoring",
    "step7b":  "True convergent amino acid substitution detection",
    "step8":   "Functional evidence (Open Targets + GTEx + DepMap)",
    "step8b":  "Functional evidence supplement (pass-through)",
    "step9":   "Compute composite scores",
    "step10":  "API ready",
    "step10b": "Regulatory divergence (AlphaGenome)",
    "step11":  "Disease annotation",
    "step11b": "Rare protective variant mapping (gnomAD)",
    "step11c": "Literature validation sanity check (PubMed)",
    "step11d": "Pathway-level convergence enrichment",
    "step12":  "Druggability assessment",
    "step12b": "P2Rank ML pocket prediction",
    "step13":  "Gene therapy feasibility",
    "step14":  "Safety pre-screen",
    "step14b": "DepMap essentiality + GTEx expression breadth",
    "step15":  "Final rescore",
    "step16":  "Pipeline complete",
}


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_state(state: dict) -> None:
    try:
        _STATE_FILE.write_text(json.dumps(state, indent=2))
    except Exception as exc:
        log.warning("Could not write pipeline state file %s: %s", _STATE_FILE, exc)


def _read_state() -> dict:
    try:
        if _STATE_FILE.exists():
            return json.loads(_STATE_FILE.read_text())
    except Exception as exc:
        log.warning("Could not read pipeline state file %s: %s — starting fresh", _STATE_FILE, exc)
    return {"status": "idle", "steps": {}}


def _mark_step(state: dict, step_name: str, status: str, elapsed: float | None = None) -> None:
    if "steps" not in state:
        state["steps"] = {}
    state["steps"][step_name] = {
        "status": status,
        "label": STEP_LABELS.get(step_name, step_name),
        "updated_at": _now_iso(),
        **({"elapsed_s": round(elapsed, 1)} if elapsed is not None else {}),
    }
    state["current_step"] = step_name
    state["updated_at"] = _now_iso()
    _write_state(state)
    run_id = state.get("run_id")
    if run_id:
        try:
            from db.models import PipelineRun
            from db.session import get_session
            with get_session() as session:
                run = session.get(PipelineRun, run_id)
                if run:
                    run.step_statuses = state.get("steps", {})
                    run.status = "running"
                    session.commit()
        except Exception as exc:
            log.warning("Could not update PipelineRun %s in DB: %s", run_id, exc)

# Add repo root to sys.path so imports work from any working directory
sys.path.insert(0, str(_REPO_ROOT))


def _load_species_registry(phenotype: str | None = None) -> list[dict]:
    registry_path = _REPO_ROOT / "config" / "species_registry.json"
    if not registry_path.exists():
        raise FileNotFoundError(
            f"Species registry not found at {registry_path}. "
            "Ensure config/species_registry.json is present."
        )
    with open(registry_path) as f:
        all_species = json.load(f)
    if not phenotype:
        return all_species
    # Keep species whose phenotype list includes the requested phenotype OR baseline/control
    filtered = []
    for sp in all_species:
        phenotypes = sp.get("phenotype", [])
        is_control = sp.get("is_control", False)
        if phenotype in phenotypes or "baseline" in phenotypes or is_control:
            filtered.append(sp)
    log.info(
        "Species filter (phenotype=%s): %d / %d species selected: %s",
        phenotype,
        len(filtered),
        len(all_species),
        [s["id"] for s in filtered],
    )
    return filtered


def step1_validate_environment() -> None:
    """Validate that all required tools are installed and the DB is reachable.

    In containerized (Nextflow) mode, bioinformatics tools live in separate
    Docker images so missing tools are logged as warnings instead of errors.
    """
    log.info("Step 1: Validating environment...")

    containerized = (
        os.environ.get("NXF_TASK_WORKDIR")
        or os.environ.get("FUSION_ENABLED")
        or os.environ.get("AWS_BATCH_JOB_ID")
        or os.environ.get("FUSION_GPU_USED") is not None
    )

    tools = [
        ("orthofinder", ["orthofinder"]),
        ("mafft",       ["mafft"]),
        ("iqtree2",     ["iqtree2", "iqtree"]),
        ("hyphy",       ["hyphy", "HYPHYMP"]),
        ("diamond",     ["diamond"]),
    ]
    missing = []
    for name, binaries in tools:
        found = any(
            subprocess.run(["which", b], capture_output=True).returncode == 0
            for b in binaries
        )
        if not found:
            missing.append(name)
            log.warning("  ✗ %s not found", name)
        else:
            log.info("  ✓ %s", name)

    if missing and not containerized:
        raise RuntimeError(f"Missing tools: {missing}. Run setup_local.sh first.")
    elif missing:
        log.info("  Containerized mode — %d tools not in this image (OK, each has its own container)", len(missing))

    # Validate DB connection
    from sqlalchemy import text
    from db.session import get_engine
    engine = get_engine()
    with engine.connect() as conn:
        conn.execute(text("SELECT 1"))
    log.info("  ✓ Database connection OK")

    # GPU check (optional — Phase 1 runs on CPU)
    try:
        import torch
        has_gpu = torch.cuda.is_available()
        log.info("  GPU available: %s", has_gpu)
    except ImportError:
        log.info("  torch not importable — GPU not checked.")


def step2_download_proteomes(species_list: list[dict], dry_run: bool = False) -> dict:
    """Download reference proteomes for all species."""
    log.info("Step 2: Downloading proteomes for %d species...", len(species_list))
    if dry_run:
        log.info("  [dry-run] skipping downloads")
        return {}

    from pipeline.layer1_sequence.download import run_downloads
    paths = run_downloads(species_list)
    log.info("  Downloaded %d proteomes.", len(paths))
    return paths


def step3_run_orthofinder(dry_run: bool = False) -> Path:
    """Run OrthoFinder on all proteomes."""
    log.info("Step 3: Running OrthoFinder (DIAMOND backend)...")
    if dry_run:
        log.info("  [dry-run] skipping OrthoFinder")
        return Path("/tmp/mock_orthofinder_results")

    from pipeline.config import get_local_storage_root, get_storage_root
    from pipeline.layer1_sequence.orthofinder import run_orthofinder

    proteomes_dir = Path(get_local_storage_root()) / "proteomes"

    if not any(proteomes_dir.glob("*.reheadered.faa")):
        _sync_proteomes_from_s3(proteomes_dir, get_storage_root())

    results_dir = run_orthofinder(proteomes_dir)
    log.info("  OrthoFinder results at: %s", results_dir)
    return results_dir


def _sync_proteomes_from_s3(proteomes_dir: Path, storage_root: str) -> None:
    """Pull all reheadered proteome FASTAs from S3 when running in a fresh container."""
    env_root = os.environ.get("BIORESILIENT_STORAGE_ROOT", "")
    if env_root.startswith("s3://"):
        storage_root = env_root
    if not storage_root.startswith("s3://"):
        return
    import boto3
    bucket = storage_root.replace("s3://", "").rstrip("/").split("/")[0]
    prefix = "proteomes/"
    proteomes_dir.mkdir(parents=True, exist_ok=True)

    s3 = boto3.client("s3")
    paginator = s3.get_paginator("list_objects_v2")
    count = 0
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            fname = key[len(prefix):]
            if not fname.endswith(".reheadered.faa") or fname.startswith("."):
                continue
            dest = proteomes_dir / fname
            if dest.exists() and dest.stat().st_size > 0:
                continue
            log.info("  S3 → %s (%s MB)", fname, round(obj["Size"] / 1e6, 1))
            s3.download_file(bucket, key, str(dest))
            count += 1
    log.info("  Synced %d proteome files from s3://%s/%s", count, bucket, prefix)


def step3b_load_orthologs(results_dir: Path, dry_run: bool = False) -> None:
    """Parse OrthoFinder output and load into DB."""
    log.info("Step 3b: Loading orthologs into database...")
    if dry_run:
        return

    from pipeline.config import get_local_storage_root, get_storage_root
    from pipeline.layer1_sequence.orthofinder import (
        flag_one_to_one_orthogroups,
        load_orthologs_to_db,
        load_sequence_map,
        parse_orthogroups,
    )

    proteomes_dir = Path(get_local_storage_root()) / "proteomes"
    if not any(proteomes_dir.glob("*.reheadered.faa")):
        _sync_proteomes_from_s3(proteomes_dir, get_storage_root())
    long_df = parse_orthogroups(results_dir, proteomes_dir)
    seq_map = load_sequence_map(proteomes_dir)

    # Build human gene map from the sequence headers
    human_gene_map = _build_human_gene_map(seq_map.get("human", {}))

    inserted = load_orthologs_to_db(long_df, seq_map, human_gene_map)
    log.info("  Loaded %d ortholog rows.", inserted)
    flagged = flag_one_to_one_orthogroups()
    log.info("  1:1 filter: %d many-to-many orthogroups excluded from downstream analysis.", flagged)


def step3c_nucleotide_conservation(species_list: list[dict], dry_run: bool = False) -> None:
    """Extract nucleotide regions (CDS/promoter/downstream) and align against human."""
    log.info("Step 3c: Nucleotide region extraction + alignment...")
    if dry_run:
        return

    from pipeline.layer1_sequence.nucleotide_scan import run_nucleotide_scan
    from pipeline.layer1_sequence.nucleotide_align import run_nucleotide_alignment

    n_regions = run_nucleotide_scan()
    log.info("  Nucleotide scan: %d regions extracted", n_regions)
    n_scored = run_nucleotide_alignment()
    log.info("  Nucleotide alignment: %d genes scored", n_scored)


def step3d_phylo_conservation(dry_run: bool = False) -> None:
    """Run phyloP / PhastCons conservation scoring on nucleotide alignments."""
    log.info("Step 3d: Phylogenetic conservation scoring (phyloP/PhastCons)...")
    if dry_run:
        return

    from pipeline.config import get_local_storage_root
    from pipeline.layer2_evolution.phylo_conservation import run_phylo_conservation
    from db.models import CandidateScore
    from db.session import get_session as _get_session

    treefile = Path(get_local_storage_root()) / "phylo" / "species.treefile"

    # Use Tier1+Tier2 gene IDs if step9 has already run; otherwise score all genes.
    tier_gene_ids: list[str] | None = None
    try:
        with _get_session() as s:
            rows = (
                s.query(CandidateScore.gene_id)
                .filter(CandidateScore.tier.in_(["Tier1", "Tier2", "Validated"]))
                .all()
            )
            if rows:
                tier_gene_ids = [r[0] for r in rows]
                log.info("  Restricting phylo scoring to %d Tier1/2 genes", len(tier_gene_ids))
    except Exception as exc:
        log.debug("step9 not yet run (no CandidateScore rows) — scoring all genes: %s", exc)

    n_scored = run_phylo_conservation(tier_gene_ids, treefile)
    log.info("  Phylogenetic conservation: %d genes scored", n_scored)


def _build_human_gene_map(human_sequences: dict[str, str]) -> dict[str, tuple[str, str]]:
    """Build {protein_id: (gene_uuid, gene_symbol)} for human proteins.

    In Phase 1 we use the UniProt accession as both gene_uuid and gene_symbol
    until we enrich with NCBI Gene IDs in Phase 2.
    """
    gene_map = {}
    for protein_id in human_sequences:
        # UniProt accession is like "human|P12345"
        acc = protein_id.split("|")[-1]
        # Use deterministic UUID from accession
        gene_uuid = str(uuid.uuid5(uuid.NAMESPACE_URL, f"human.{acc}"))
        gene_map[protein_id] = (gene_uuid, acc)
    return gene_map


def step4_alignment_and_divergence(
    species_list: list[dict],
    dry_run: bool = False,
) -> tuple[dict, dict]:
    """Align orthogroups and compute divergent motifs. Applies Gate 1 (before MAFFT) and Gate 2 (after divergence)."""
    log.info("Step 4: Gate 1 filter → MAFFT alignment → sliding window divergence → Gate 2 filter...")
    if dry_run:
        return {}, {}

    from pipeline.config import get_local_storage_root

    from pipeline.config import get_thresholds
    from pipeline.layer1_sequence.alignment import (
        align_all_orthogroups,
        filter_orthogroups_by_global_identity,
        load_orthogroup_sequences_from_db,
    )
    from pipeline.layer1_sequence.divergence import (
        apply_branch_length_weighting,
        branch_length_weights_from_tree,
        filter_by_independent_lineages,
        load_motifs_to_db,
        run_divergence_pipeline,
        update_sequence_identities,
    )

    orthogroups = load_orthogroup_sequences_from_db()
    log.info("  %d orthogroups loaded.", len(orthogroups))

    thresholds = get_thresholds()
    min_divergent = int(thresholds.get("divergence_min_species", 2))
    identity_max = float(thresholds.get("divergence_identity_max", 0.85))
    divergence_pct_min = 100.0 * (1.0 - identity_max)  # e.g. 0.85 → 15% divergence
    orthogroups_gate1 = filter_orthogroups_by_global_identity(
        orthogroups,
        min_divergent_species=min_divergent,
        divergence_pct_min=divergence_pct_min,
    )
    log.info("  %d orthogroups to align (after Gate 1).", len(orthogroups_gate1))

    aligned = align_all_orthogroups(orthogroups_gate1)
    update_sequence_identities(aligned)

    motifs_by_og_raw = run_divergence_pipeline(aligned)

    # Optional: scale divergence_score by inverse evolutionary distance.
    # Requires the IQ-TREE species tree from step 5 (chicken-and-egg on first run;
    # set divergence_apply_branch_norm: true in environment.yml for reruns).
    if thresholds.get("divergence_apply_branch_norm", False):
        from pipeline.config import get_local_storage_root as _glsr
        _treefile = Path(_glsr()) / "phylo" / "species.treefile"
        _weights = branch_length_weights_from_tree(_treefile)
        if _weights:
            motifs_by_og_raw = apply_branch_length_weighting(motifs_by_og_raw, _weights)
        else:
            log.info(
                "Branch-length weighting requested but tree unavailable — skipped. "
                "Re-run step4 after step5 to apply weighting."
            )

    species_to_lineage = {
        s["id"]: s.get("lineage_group", "Other")
        for s in species_list
        if not s.get("is_control", False) and "baseline" not in s.get("phenotype", [])
    }
    min_lineages = int(thresholds.get("convergence_min_lineages", 2))
    motifs_by_og = filter_by_independent_lineages(
        motifs_by_og_raw,
        species_to_lineage,
        min_lineages=min_lineages,
    )
    all_motifs_flat = [m for motifs in motifs_by_og.values() for m in motifs]
    load_motifs_to_db(all_motifs_flat, {})

    aligned_for_hyphy = {og: aligned[og] for og in motifs_by_og}
    log.info("  After Gate 2: %d orthogroups for HyPhy.", len(aligned_for_hyphy))

    # Save aligned orthogroups to disk for resume capability
    import pickle
    from pipeline.config import get_local_storage_root, sync_to_s3
    _cache_path = Path(get_local_storage_root()) / "aligned_orthogroups.pkl"
    with open(_cache_path, "wb") as _f:
        pickle.dump({"aligned": aligned, "motifs_by_og": motifs_by_og}, _f)
    log.info("  Saved aligned orthogroups cache → %s", _cache_path)
    sync_to_s3(_cache_path, "cache/aligned_orthogroups.pkl")

    return aligned, motifs_by_og


def step4b_domain_and_consequence(dry_run: bool = False) -> None:
    """Annotate divergent motifs with Pfam domain context and AlphaMissense functional scores."""
    log.info("Step 4b: Domain annotation (Pfam/InterPro) + functional consequence (AlphaMissense)...")
    if dry_run:
        return

    from pipeline.layer1_sequence.pfam import run_pfam_pipeline
    from pipeline.layer1_sequence.alphamissense import run_alphamissense_pipeline

    n_domains = run_pfam_pipeline()
    log.info("  Pfam: %d motifs in functional domains.", n_domains)

    n_scored = run_alphamissense_pipeline()
    log.info("  AlphaMissense: %d motifs received consequence scores.", n_scored)


def step4c_esm1v(dry_run: bool = False) -> None:
    """Score divergent motifs with ESM-1v variant effect predictions."""
    log.info("Step 4c: ESM-1v structural variant effect scoring...")
    if dry_run:
        return

    from pipeline.layer1_sequence.esm1v import run_esm1v_pipeline
    n = run_esm1v_pipeline()
    log.info("  ESM-1v: %d motifs scored.", n)


def step4d_variant_direction(dry_run: bool = False) -> None:
    """Classify each divergent motif as GoF / LoF / neutral using AM + ESM-1v + LOEUF."""
    log.info("Step 4d: Variant direction inference (GoF / LoF)...")
    if dry_run:
        return

    from pipeline.layer1_sequence.variant_direction import run_variant_direction_pipeline
    n = run_variant_direction_pipeline()
    log.info("  Variant direction: %d motifs classified.", n)


def step5_phylogenetic_tree(
    aligned_orthogroups: dict,
    dry_run: bool = False,
) -> Path:
    """Build species tree with IQ-TREE2."""
    log.info("Step 5: Building species phylogenetic tree with IQ-TREE2...")
    if dry_run:
        return Path("/tmp/mock.treefile")

    from pipeline.config import get_local_storage_root, sync_to_s3
    from pipeline.layer2_evolution.phylo_tree import build_concatenated_alignment, run_iqtree, _max_ogs

    concat = build_concatenated_alignment(aligned_orthogroups, max_ogs=_max_ogs())
    if concat is None:
        raise RuntimeError("Not enough single-copy orthogroups to build species tree.")

    treefile = run_iqtree(concat)
    sync_to_s3(treefile, "cache/species.treefile")
    return treefile


def step6_evolutionary_selection(
    aligned_orthogroups: dict,
    motifs_by_og: dict,
    treefile: Path,
    dry_run: bool = False,
) -> None:
    """Run MEME episodic selection analysis on candidate orthogroups.

    Tries codon-guided HyPhy MEME first (statistically rigorous episodic
    selection). Falls back to protein divergence proxy for any orthogroup
    where CDS fetching fails (e.g. UniProt accessions without NCBI CDS links).
    """
    log.info("Step 6: Running MEME episodic selection analysis...")
    if dry_run:
        return

    from pipeline.layer2_evolution.meme_selection import run_meme_pipeline
    from pipeline.layer2_evolution.selection import build_gene_og_map, load_selection_scores

    gene_by_og = build_gene_og_map()

    # Skip OGs already written to DB — safe to resume mid-run
    from db.models import EvolutionScore
    from db.session import get_session
    with get_session() as session:
        already_scored_genes = {
            r[0] for r in session.query(EvolutionScore.gene_id)
            .filter(EvolutionScore.dnds_ratio.isnot(None))
            .all()
        }
    already_scored_ogs = {og for og, gid in gene_by_og.items() if gid in already_scored_genes}
    if already_scored_ogs:
        log.info("Skipping %d already-scored orthogroups; resuming from checkpoint.",
                 len(already_scored_ogs))
        aligned_orthogroups = {k: v for k, v in aligned_orthogroups.items()
                                if k not in already_scored_ogs}

    # Flush callback: write to DB every 500 OGs so progress is saved incrementally
    def _flush(batch: dict[str, dict]) -> None:
        load_selection_scores(batch, gene_by_og)

    selection_results = run_meme_pipeline(
        aligned_orthogroups, motifs_by_og, treefile, flush_callback=_flush
    )
    if selection_results:
        load_selection_scores(selection_results, gene_by_og)


def step6b_fel_busted(dry_run: bool = False) -> None:
    """Run FEL (pervasive selection) and BUSTED (gene-level episodic) supplementary tests.

    Reuses codon alignments from step6. Skips gracefully if alignments are missing.
    Updates EvolutionScore.fel_sites and EvolutionScore.busted_pvalue.
    """
    log.info("Step 6b: FEL + BUSTED supplementary selection tests...")
    if dry_run:
        return

    from pipeline.config import get_local_storage_root
    from pipeline.layer2_evolution.meme_selection import run_fel_busted_pipeline
    from pipeline.layer2_evolution.selection import build_gene_og_map, load_fel_busted_scores

    import pickle
    _cache_path = Path(get_local_storage_root()) / "aligned_orthogroups.pkl"
    if not _cache_path.exists():
        log.warning("  Aligned orthogroups cache not found; skipping step 6b.")
        return

    with open(_cache_path, "rb") as _f:
        cached = pickle.load(_f)

    aligned = cached.get("aligned", {})
    motifs_by_og = cached.get("motifs_by_og", {})

    treefile = Path(get_local_storage_root()) / "phylo" / "species.treefile"
    if not treefile.exists():
        log.warning("  Species treefile not found; skipping step 6b.")
        return

    gene_by_og = build_gene_og_map()

    def _flush_fel_busted(batch: dict) -> None:
        load_fel_busted_scores(batch, gene_by_og)

    run_fel_busted_pipeline(aligned, motifs_by_og, treefile, flush_callback=_flush_fel_busted)
    log.info("  FEL/BUSTED: complete.")


def step6c_relax(dry_run: bool = False) -> None:
    """Run RELAX to detect branch-specific rate acceleration in resilient species.

    Reuses codon alignments from step6. Skips gracefully if alignments are missing.
    Updates EvolutionScore.relax_k and EvolutionScore.relax_pvalue.
    """
    log.info("Step 6c: RELAX branch acceleration tests...")
    if dry_run:
        return

    from pipeline.config import get_local_storage_root
    from pipeline.layer2_evolution.meme_selection import run_relax_pipeline
    from pipeline.layer2_evolution.selection import build_gene_og_map, load_relax_scores

    import pickle
    _cache_path = Path(get_local_storage_root()) / "aligned_orthogroups.pkl"
    if not _cache_path.exists():
        log.warning("  Aligned orthogroups cache not found; skipping step 6c.")
        return

    with open(_cache_path, "rb") as f:
        cached = pickle.load(f)

    aligned = cached.get("aligned", {})
    motifs_by_og = cached.get("motifs_by_og", {})

    treefile = Path(get_local_storage_root()) / "phylo" / "species.treefile"
    if not treefile.exists():
        log.warning("  Species treefile not found; skipping step 6c.")
        return

    gene_by_og = build_gene_og_map()

    def _flush_relax(batch: dict) -> None:
        load_relax_scores(batch, gene_by_og)

    run_relax_pipeline(aligned, motifs_by_og, treefile, flush_callback=_flush_relax)
    log.info("  RELAX: complete.")


def step7_convergence(dry_run: bool = False) -> None:
    log.info("Step 7: Computing convergence detection + PhyloP enrichment...")
    if dry_run:
        return

    from pipeline.layer2_evolution.convergence import (
        build_chrom_map,
        enrich_phylop_scores,
        run_convergence_pipeline,
        run_convergence_permutation_test,
    )
    from db.models import EvolutionScore, PhyloConservationScore
    from db.session import get_session
    from pipeline.config import get_tool_config

    run_convergence_pipeline()

    # Fix 2: run permutation-based null model to assign significance p-values.
    # Skip if already populated (idempotent) — check a sample of rows.
    cfg = get_tool_config()
    n_iter = int(cfg.get("convergence_permutation_iterations", 200))
    with get_session() as session:
        sample = (
            session.query(EvolutionScore.convergence_pval)
            .filter(EvolutionScore.convergence_pval.isnot(None))
            .limit(1)
            .first()
        )
    if sample is None:
        log.info("Running convergence permutation test (%d iterations)...", n_iter)
        run_convergence_permutation_test(n_iterations=n_iter)
    else:
        log.info(
            "convergence_pval already populated — skipping permutation test "
            "(delete rows or set convergence_pval=NULL to re-run)."
        )

    # Seed evolution_score.phylop_score from step 3d output (phylo_conservation_score)
    # before calling the slower UCSC API. This avoids re-querying the ~11K genes that
    # step 3d already resolved, reducing the sequential UCSC loop from ~12K iterations
    # to only the ~1-2K genes with no existing conservation data.
    with get_session() as session:
        pcs_rows = session.query(
            PhyloConservationScore.gene_id,
            PhyloConservationScore.cds_phylo_score,
        ).filter(PhyloConservationScore.cds_phylo_score.isnot(None)).all()

    if pcs_rows:
        log.info(
            "Seeding phylop_score from step 3d data for %d genes before UCSC enrichment...",
            len(pcs_rows),
        )
        with get_session() as session:
            ev_map = {ev.gene_id: ev for ev in session.query(EvolutionScore).all()}
            seeded = 0
            for gene_id, cds_score in pcs_rows:
                ev = ev_map.get(gene_id)
                if ev is None:
                    ev = EvolutionScore(gene_id=gene_id)
                    session.add(ev)
                    ev_map[gene_id] = ev
                # Only seed if not already set by run_convergence_pipeline
                if ev.phylop_score is None:
                    ev.phylop_score = cds_score
                    seeded += 1
            session.commit()
        log.info("Seeded phylop_score for %d genes from step 3d.", seeded)

    # Only UCSC-query genes that still have no phylop_score after the seed above
    with get_session() as session:
        missing_ids = [
            ev.gene_id
            for ev in session.query(EvolutionScore).all()
            if ev.phylop_score is None
        ]

    if missing_ids:
        log.info(
            "Building NCBI chrom map for %d genes missing phylop_score...", len(missing_ids)
        )
        chrom_map = build_chrom_map(gene_ids=missing_ids)
        if chrom_map:
            enrich_phylop_scores(chrom_map)
    else:
        log.info("All genes have phylop_score from step 3d — skipping UCSC enrichment.")


def step7b_convergent_aa(dry_run: bool = False) -> None:
    """Detect true convergent amino acid substitutions at each motif position.

    Goes beyond lineage counting: checks whether the *same* (or biochemically
    equivalent) amino acid change evolved independently in multiple lineages.
    Stores convergent_aa_count on DivergentMotif.
    """
    log.info("Step 7b: True convergent amino acid substitution detection...")
    if dry_run:
        return

    from pipeline.layer2_evolution.convergent_aa import run_convergent_aa_pipeline
    n = run_convergent_aa_pipeline()
    log.info("  Convergent AA: %d motifs with true convergence (≥2 lineages, same substitution).", n)


def step8_functional_evidence(phenotype: str = "cancer_resistance", dry_run: bool = False) -> None:
    """Gather functional evidence (Open Targets + GTEx + DepMap) for all genes.

    Replaces GEO/DESeq2 (old step 8) and Bgee (old step 8b) with three
    well-curated, phenotype-configurable data sources:

      1. Open Targets Platform  — gene-disease association scores
      2. GTEx tissue expression — median TPM in phenotype-relevant tissues
      3. DepMap essentiality    — CRISPR Chronos scores (cancer/dna_repair only)

    Config: config/functional_evidence_config.json
    Generalised for any phenotype via per-phenotype EFO/MONDO disease IDs,
    GTEx tissue lists, and optional DepMap enable flag.
    """
    log.info("Step 8: Functional evidence (OT + GTEx + DepMap), phenotype=%r...", phenotype)
    if dry_run:
        return

    from pipeline.layer1_sequence.functional_evidence import run_functional_evidence
    summary = run_functional_evidence(phenotype=phenotype)
    log.info("  Functional evidence complete: %s", summary)


# Keep old name as alias for backward compatibility with any direct callers.
def step8_expression(species_list: list[dict], dry_run: bool = False) -> None:
    """Deprecated: delegates to step8_functional_evidence."""
    phenotype = "cancer_resistance"
    step8_functional_evidence(phenotype=phenotype, dry_run=dry_run)


def step8b_bgee(dry_run: bool = False) -> None:
    """Step 8b: No-op pass-through (unified with step 8 functional evidence)."""
    log.info("Step 8b: Skipping — functional evidence is now unified in step 8.")


def step9_composite_score(dry_run: bool = False) -> None:
    """Assemble composite scores and assign tiers."""
    log.info("Step 9: Assembling composite scores...")
    if dry_run:
        return

    from pipeline.scoring import get_top_candidates, run_scoring
    run_scoring(phase="phase1")

    top = get_top_candidates(n=10)
    log.info("  Top 10 candidates:")
    for i, c in enumerate(top, 1):
        log.info("    %2d. %-12s  composite=%.3f  %s",
                 i, c["gene_symbol"], c["composite_score"], c["tier"])


def step10_start_api(dry_run: bool = False) -> None:
    """Start the FastAPI server (background process)."""
    log.info("Step 10: API ready. Start with: uvicorn api.main:app --host 0.0.0.0 --port 8000")
    if dry_run:
        return


def step10b_alphagenome(dry_run: bool = False) -> None:
    """Track B: AlphaGenome regulatory divergence scan.

    Runs on expression-positive genes not fully explained by Track A (protein divergence).
    Requires ALPHAGENOME_API_KEY environment variable. Silently skips if key is absent.
    """
    log.info("Step 10b: Track B — AlphaGenome regulatory divergence...")
    if dry_run:
        return

    from pipeline.layer_regulatory.alphagenome import run_alphagenome_track
    n = run_alphagenome_track()
    if n:
        log.info("  Track B: %d regulatory divergence rows written.", n)


def _get_tier12_gene_ids(trait_id: str = "") -> list[str]:
    """Return gene IDs for Tier1 and Tier2 candidates (for funnel gating of Phase 2 layers).

    Phase 1 stores CandidateScore rows with trait_id='' (empty string).
    If no rows match the requested trait_id, fall back to trait_id='' so that
    Phase 2 steps always see the scored genes regardless of how Phase 1 was run.
    """
    from db.models import CandidateScore
    from db.session import get_session
    with get_session() as session:
        q = session.query(CandidateScore.gene_id).filter(
            CandidateScore.tier.in_(["Tier1", "Tier2"])
        )
        if trait_id:
            q = q.filter(CandidateScore.trait_id == trait_id)
        rows = q.all()
        if not rows and trait_id:
            # Phase 1 may have stored scores with empty trait_id — fall back
            log.info(
                "_get_tier12_gene_ids: no rows for trait_id=%r, falling back to trait_id=''", trait_id
            )
            rows = (
                session.query(CandidateScore.gene_id)
                .filter(
                    CandidateScore.tier.in_(["Tier1", "Tier2"]),
                    CandidateScore.trait_id == "",
                )
                .all()
            )
    return [r[0] for r in rows]


def step11_disease_annotation(dry_run: bool = False, trait_id: str = "") -> None:
    """Layer 3: OpenTargets, GWAS, gnomAD, IMPC, Human Protein Atlas (Tier1+Tier2 only).

    All five external API calls are independent — they run concurrently via
    ThreadPoolExecutor to reduce wall-clock time from ~30 min to ~6-8 min
    (bounded by the slowest API, typically GWAS Catalog).
    """
    log.info("Step 11: Disease annotation (Layer 3)...")
    if dry_run:
        return
    gene_ids = _get_tier12_gene_ids(trait_id)
    if not gene_ids:
        log.warning("  No Tier1/Tier2 candidates; skipping disease annotation.")
        return

    from concurrent.futures import ThreadPoolExecutor, as_completed
    from pipeline.layer3_disease import opentargets, gwas, gnomad, impc, protein_atlas, pathways

    # Independent API sources — run in parallel (I/O bound, not CPU bound)
    tasks = {
        "opentargets": lambda: opentargets.annotate_genes_opentargets(gene_ids),
        "gwas":        lambda: gwas.annotate_genes_gwas(gene_ids),
        "gnomad":      lambda: gnomad.annotate_genes_gnomad(gene_ids),
        "impc":        lambda: impc.annotate_genes_impc(gene_ids),
        "protein_atlas": lambda: protein_atlas.annotate_genes_protein_atlas(gene_ids),
    }

    with ThreadPoolExecutor(max_workers=len(tasks)) as executor:
        future_to_name = {executor.submit(fn): name for name, fn in tasks.items()}
        for future in as_completed(future_to_name):
            name = future_to_name[future]
            try:
                future.result()
                log.info("  %s annotation done.", name)
            except Exception as exc:
                log.warning("  %s annotation failed (non-fatal, continuing): %s", name, exc)

    # Pathways enrichment runs after OpenTargets completes (uses disease IDs)
    try:
        n_pathways = pathways.annotate_genes_pathways(gene_ids)
        log.info("  Pathways/GO annotated for %d genes.", n_pathways)
    except Exception as exc:
        log.warning("  Pathways annotation failed (non-fatal): %s", exc)


def step11b_rare_variants(dry_run: bool = False, trait_id: str = "") -> None:
    """Map divergent motif positions to rare protective variants in gnomAD (PCSK9 paradigm)."""
    log.info("Step 11b: Rare protective variant mapping (gnomAD + GWAS Catalog)...")
    if dry_run:
        return
    gene_ids = _get_tier12_gene_ids(trait_id)
    if not gene_ids:
        log.warning("  No Tier1/Tier2 candidates; skipping rare variant mapping.")
        return
    from pipeline.layer3_disease.rare_variants import run_rare_variants_pipeline
    n = run_rare_variants_pipeline(gene_ids)
    log.info("  Rare variants: %d genes with protective variant matches.", n)


def step11c_literature(dry_run: bool = False, trait_id: str = "") -> None:
    """Check top candidates against known resilience/longevity literature (PubMed sanity check)."""
    log.info("Step 11c: Literature validation (PubMed sanity check)...")
    if dry_run:
        return
    gene_ids = _get_tier12_gene_ids(trait_id)
    if not gene_ids:
        log.warning("  No Tier1/Tier2 candidates; skipping literature check.")
        return
    from pipeline.layer3_disease.literature import run_literature_pipeline
    n = run_literature_pipeline(gene_ids)
    log.info("  Literature: %d genes with citations in resilience literature.", n)


def step11d_pathway_convergence(dry_run: bool = False) -> None:
    """Compute pathway-level convergence enrichment across all candidate genes."""
    log.info("Step 11d: Pathway-level convergence scoring...")
    if dry_run:
        return
    from pipeline.layer3_disease.pathway_convergence import run_pathway_convergence_pipeline
    run_pathway_convergence_pipeline()


def step12_druggability(dry_run: bool = False, trait_id: str = "") -> None:
    """Layer 4: AlphaFold structures, fpocket, ChEMBL, CanSAR, peptide (Tier1+Tier2 only)."""
    log.info("Step 12: Druggability (Layer 4)...")
    if dry_run:
        return
    gene_ids = _get_tier12_gene_ids(trait_id)
    if not gene_ids:
        return
    from pipeline.layer4_druggability import structure, pockets, chembl, cansar, peptide
    gene_to_pdb = structure.ensure_structures_for_genes(gene_ids)
    pockets.annotate_pockets(gene_to_pdb)
    chembl.annotate_genes_chembl(gene_ids)
    cansar.annotate_genes_cansar(gene_ids)
    peptide.annotate_motifs_peptide(None)


def step12b_p2rank(dry_run: bool = False, trait_id: str = "") -> None:
    """Run P2Rank ML pocket prediction on AlphaFold structures."""
    log.info("Step 12b: P2Rank ML pocket prediction...")
    if dry_run:
        return
    gene_ids = _get_tier12_gene_ids(trait_id)
    if not gene_ids:
        return
    from pipeline.layer4_druggability.structure import ensure_structures_for_genes
    from pipeline.layer4_druggability.p2rank import run_p2rank_pipeline
    gene_to_pdb = ensure_structures_for_genes(gene_ids)
    n = run_p2rank_pipeline(gene_to_pdb, gene_ids)
    log.info("  P2Rank: %d genes with pocket predictions.", n)


def step13_gene_therapy(dry_run: bool = False, trait_id: str = "") -> None:
    """Layer 5: AAV + CRISPR (Tier1 only)."""
    log.info("Step 13: Gene therapy (Layer 5)...")
    if dry_run:
        return
    from db.models import CandidateScore
    from db.session import get_session
    with get_session() as session:
        q = session.query(CandidateScore.gene_id).filter(CandidateScore.tier == "Tier1")
        if trait_id:
            q = q.filter(CandidateScore.trait_id == trait_id)
        tier1 = q.all()
    gene_ids = [r[0] for r in tier1]
    if not gene_ids:
        return
    from pipeline.layer5_gene_therapy import aav, crispr
    aav.annotate_genes_aav(gene_ids)
    crispr.annotate_genes_crispr(gene_ids)


def step14_safety(dry_run: bool = False, trait_id: str = "") -> None:
    """Layer 6: PheWAS, STRING, selectivity (Tier1+Tier2 only)."""
    log.info("Step 14: Safety (Layer 6)...")
    if dry_run:
        return
    gene_ids = _get_tier12_gene_ids(trait_id)
    if not gene_ids:
        return
    from pipeline.layer6_safety import phewas, network, selectivity
    phewas.annotate_genes_phewas(gene_ids)
    network.annotate_genes_network(gene_ids)
    selectivity.annotate_genes_selectivity(gene_ids)


def step14b_depmap_gtex(dry_run: bool = False, trait_id: str = "") -> None:
    """Annotate safety with DepMap CRISPR essentiality and GTEx expression breadth."""
    log.info("Step 14b: DepMap essentiality + GTEx expression breadth...")
    if dry_run:
        return
    gene_ids = _get_tier12_gene_ids(trait_id)
    if not gene_ids:
        return
    from pipeline.layer6_safety.depmap import run_depmap_pipeline
    from pipeline.layer6_safety.gtex import run_gtex_pipeline
    n_depmap = run_depmap_pipeline(gene_ids)
    log.info("  DepMap: %d genes annotated.", n_depmap)
    n_gtex = run_gtex_pipeline(gene_ids)
    log.info("  GTEx: %d genes annotated.", n_gtex)


def step15_rescore(dry_run: bool = False) -> None:
    """Re-run composite scoring with Phase 2 weights, then apply control species penalty."""
    log.info("Step 15: Rescore (phase2 weights) + control divergence penalty...")
    if dry_run:
        return
    from pipeline.scoring import get_top_candidates, run_scoring
    run_scoring(phase="phase2")

    from pipeline.layer2_evolution.convergence import (
        apply_control_divergence_penalty,
        compute_control_divergence_fractions,
    )
    fractions = compute_control_divergence_fractions()
    if fractions:
        apply_control_divergence_penalty(fractions)

    top = get_top_candidates(n=10)
    for i, c in enumerate(top, 1):
        log.info("    %2d. %-12s  composite=%.3f  %s",
                 i, c["gene_symbol"], c["composite_score"], c["tier"])


def step16_start_api(dry_run: bool = False) -> None:
    """Pipeline complete; API can be started separately."""
    log.info("Step 16: Pipeline complete. Start API: uvicorn api.main:app --host 0.0.0.0 --port 8000")
    if dry_run:
        return


STEPS = {
    "step1":   step1_validate_environment,
    "step2":   step2_download_proteomes,
    "step3":   step3_run_orthofinder,
    "step3b":  step3b_load_orthologs,
    "step3c":  step3c_nucleotide_conservation,
    "step4":   step4_alignment_and_divergence,
    "step4b":  step4b_domain_and_consequence,
    "step4c":  step4c_esm1v,
    "step4d":  step4d_variant_direction,
    "step5":   step5_phylogenetic_tree,
    "step3d":  step3d_phylo_conservation,
    "step6":   step6_evolutionary_selection,
    "step6b":  step6b_fel_busted,
    "step6c":  step6c_relax,
    "step7":   step7_convergence,
    "step7b":  step7b_convergent_aa,
    "step8":   step8_functional_evidence,
    "step8b":  step8b_bgee,
    "step9":   step9_composite_score,
    "step10":  step10_start_api,
    "step10b": step10b_alphagenome,
    "step11":  step11_disease_annotation,
    "step11b": step11b_rare_variants,
    "step11c": step11c_literature,
    "step11d": step11d_pathway_convergence,
    "step12":  step12_druggability,
    "step12b": step12b_p2rank,
    "step13":  step13_gene_therapy,
    "step14":  step14_safety,
    "step14b": step14b_depmap_gtex,
    "step15":  step15_rescore,
    "step16":  step16_start_api,
}


def run_pipeline(
    species_list: list[dict],
    resume_from: str = "step1",
    only_steps: list[str] | None = None,
    dry_run: bool = False,
    trait_id: str = "",
) -> None:
    step_order = list(STEPS.keys())
    if only_steps:
        step_order = [s for s in step_order if s in only_steps]
    elif resume_from:
        start_idx = step_order.index(resume_from) if resume_from in step_order else 0
        step_order = step_order[start_idx:]

    start_time = time.time()
    log.info("=" * 60)
    log.info("BioResilient AI — Phase 1 Pipeline")
    log.info("Steps: %s", step_order)
    log.info("Dry run: %s", dry_run)
    if trait_id:
        log.info("Phenotype / trait: %s", trait_id)
    log.info("=" * 60)

    # Shared state passed between steps
    aligned_orthogroups: dict = {}
    motifs_by_og: dict = {}
    results_dir: Path | None = None
    treefile: Path | None = None

    # If resuming past step3, recover results_dir from disk
    from pipeline.layer1_sequence.orthofinder import _orthofinder_output_dir
    from pipeline.config import get_local_storage_root
    _proteomes_dir = Path(get_local_storage_root()) / "proteomes"
    _existing = _orthofinder_output_dir(_proteomes_dir)
    if _existing is not None:
        results_dir = _existing

    # If resuming past step4, recover aligned_orthogroups from disk cache.
    # This triggers when:
    #   - resume_from is a step after step4 (e.g. --resume-from step5), OR
    #   - only_steps contains a step that needs aligned_orthogroups (step5, step6, step7...)
    _EARLY_STEPS = {"step1", "step2", "step3", "step3b", "step4", "step4b", "step4c", "step4d"}
    _needs_pkl = (
        resume_from not in _EARLY_STEPS
        or (only_steps and not set(only_steps).issubset(_EARLY_STEPS))
    )
    if _needs_pkl:
        _pkl_path = Path(get_local_storage_root()) / "aligned_orthogroups.pkl"

        # Auto-sync from S3 if missing locally — happens after instance stop/start
        if not _pkl_path.exists():
            _s3_key = "cache/aligned_orthogroups.pkl"
            log.info(
                "aligned_orthogroups.pkl not found locally — attempting S3 sync from %s ...",
                _s3_key,
            )
            try:
                from pipeline.config import sync_from_s3
                sync_from_s3(_s3_key, _pkl_path)
                log.info("Synced aligned_orthogroups.pkl from S3 (%s)", _pkl_path)
            except Exception as _sync_err:
                log.warning("Could not sync aligned_orthogroups.pkl from S3: %s", _sync_err)

        if _pkl_path.exists():
            try:
                with open(_pkl_path, "rb") as _f:
                    _data = pickle.load(_f)
                aligned_orthogroups = _data["aligned"]
                motifs_by_og = _data["motifs_by_og"]
                log.info("Recovered %d aligned orthogroups from disk cache", len(aligned_orthogroups))
            except Exception as _e:
                log.warning("Could not load aligned orthogroups cache: %s", _e)

        if not aligned_orthogroups:
            # Last resort: reconstruct from DB (sequences only, no alignment gaps)
            try:
                from db.models import Ortholog
                from db.session import get_session as _gs
                _aligned: dict = {}
                with _gs() as _s:
                    for row in _s.query(Ortholog).all():
                        og = row.orthofinder_og or row.gene_id
                        if og not in _aligned:
                            _aligned[og] = {}
                        if row.protein_seq:
                            _aligned[og][row.species_id] = row.protein_seq
                aligned_orthogroups = _aligned
                motifs_by_og = {og: [] for og in _aligned}
                log.info("Recovered %d orthogroups from DB for resume", len(aligned_orthogroups))
            except Exception as _e:
                log.warning("Could not recover aligned_orthogroups from DB: %s", _e)

    # If resuming past step5, recover treefile from disk or S3
    _tree_candidate = Path(get_local_storage_root()) / "phylo" / "species.treefile"
    if not _tree_candidate.exists():
        _tree_candidate.parent.mkdir(parents=True, exist_ok=True)
        try:
            from pipeline.config import sync_from_s3
            sync_from_s3("cache/species.treefile", _tree_candidate)
            log.info("Synced species.treefile from S3")
        except Exception:
            pass
    if _tree_candidate.exists():
        treefile = _tree_candidate

    # Initialise pipeline state file and PipelineRun row for API/DB tracking
    state = {
        "status": "running",
        "started_at": _now_iso(),
        "updated_at": _now_iso(),
        "dry_run": dry_run,
        "steps": {
            s: {"status": "pending", "label": STEP_LABELS.get(s, s)}
            for s in step_order
        },
    }
    run_id = os.environ.get("PIPELINE_RUN_ID")
    try:
        from db.models import PipelineRun
        from db.session import get_session
        with get_session() as session:
            if run_id:
                run = session.get(PipelineRun, run_id)
                if run:
                    run.species_ids = [s.get("id") for s in species_list] if species_list else []
                    run.step_statuses = state["steps"]
                    run.status = "running"
                    session.commit()
            else:
                run = PipelineRun(
                    status="running",
                    species_ids=[s.get("id") for s in species_list] if species_list else [],
                    step_statuses=state["steps"],
                    phase="phase1",
                    pid=os.getpid(),
                )
                session.add(run)
                session.commit()
                run_id = run.id
            state["run_id"] = run_id
    except Exception as exc:
        log.warning("Could not create/update PipelineRun row in DB: %s", exc)
    _write_state(state)

    # Also add a file handler so the log streams to pipeline.log for SSE tailing
    file_handler = logging.FileHandler(_LOG_FILE, mode="a", encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s  %(levelname)-7s  %(message)s", datefmt="%H:%M:%S"))
    logging.getLogger().addHandler(file_handler)

    step_name = step_order[0] if step_order else "unknown"
    try:
        for step_name in step_order:
            log.info("")
            step_start = time.time()
            _mark_step(state, step_name, "running")

            if step_name == "step1":
                step1_validate_environment()

            elif step_name == "step2":
                step2_download_proteomes(species_list, dry_run=dry_run)

            elif step_name == "step3":
                results_dir = step3_run_orthofinder(dry_run=dry_run)

            elif step_name == "step3b":
                if results_dir:
                    step3b_load_orthologs(results_dir, dry_run=dry_run)
                else:
                    raise RuntimeError(
                        "step3b: results_dir is None — OrthoFinder has not been run "
                        "or its output directory could not be found. "
                        "Run step3 first, or set results_dir manually."
                    )

            elif step_name == "step3c":
                step3c_nucleotide_conservation(species_list, dry_run=dry_run)

            elif step_name == "step3d":
                step3d_phylo_conservation(dry_run=dry_run)

            elif step_name == "step4":
                aligned_orthogroups, motifs_by_og = step4_alignment_and_divergence(
                    species_list, dry_run=dry_run
                )

            elif step_name == "step4b":
                step4b_domain_and_consequence(dry_run=dry_run)

            elif step_name == "step4c":
                step4c_esm1v(dry_run=dry_run)

            elif step_name == "step4d":
                step4d_variant_direction(dry_run=dry_run)

            elif step_name == "step5":
                treefile = step5_phylogenetic_tree(aligned_orthogroups, dry_run=dry_run)

            elif step_name == "step6":
                if not treefile:
                    raise RuntimeError(
                        "step6: treefile is None — phylogenetic tree has not been built. "
                        "Run step5 first."
                    )
                if not motifs_by_og:
                    raise RuntimeError(
                        "step6: motifs_by_og is empty — no divergent motifs found in step4. "
                        "Ensure step4 completed successfully with motifs in the database."
                    )
                aligned_for_hyphy = {
                    og: aligned_orthogroups[og]
                    for og in motifs_by_og
                    if og in aligned_orthogroups
                }
                step6_evolutionary_selection(
                    aligned_for_hyphy, motifs_by_og, treefile, dry_run=dry_run
                )

            elif step_name == "step6b":
                step6b_fel_busted(dry_run=dry_run)

            elif step_name == "step6c":
                step6c_relax(dry_run=dry_run)

            elif step_name == "step7":
                step7_convergence(dry_run=dry_run)

            elif step_name == "step7b":
                step7b_convergent_aa(dry_run=dry_run)

            elif step_name == "step8":
                step8_functional_evidence(phenotype=phenotype, dry_run=dry_run)

            elif step_name == "step8b":
                step8b_bgee(dry_run=dry_run)

            elif step_name == "step9":
                step9_composite_score(dry_run=dry_run)

            elif step_name == "step10":
                step10_start_api(dry_run=dry_run)

            elif step_name == "step10b":
                step10b_alphagenome(dry_run=dry_run)

            elif step_name == "step11":
                step11_disease_annotation(dry_run=dry_run, trait_id=trait_id)

            elif step_name == "step11b":
                step11b_rare_variants(dry_run=dry_run, trait_id=trait_id)

            elif step_name == "step11c":
                step11c_literature(dry_run=dry_run, trait_id=trait_id)

            elif step_name == "step11d":
                step11d_pathway_convergence(dry_run=dry_run)

            elif step_name == "step12":
                step12_druggability(dry_run=dry_run, trait_id=trait_id)

            elif step_name == "step12b":
                step12b_p2rank(dry_run=dry_run, trait_id=trait_id)

            elif step_name == "step13":
                step13_gene_therapy(dry_run=dry_run, trait_id=trait_id)

            elif step_name == "step14":
                step14_safety(dry_run=dry_run, trait_id=trait_id)

            elif step_name == "step14b":
                step14b_depmap_gtex(dry_run=dry_run, trait_id=trait_id)

            elif step_name == "step15":
                step15_rescore(dry_run=dry_run)

            elif step_name == "step16":
                step16_start_api(dry_run=dry_run)

            elapsed = time.time() - step_start
            log.info("  %s complete (%.1fs)", step_name, elapsed)
            _mark_step(state, step_name, "complete", elapsed)

            # Write step cache so reruns can skip completed steps
            _cache_path = _REPO_ROOT / "pipeline_cache.json"
            try:
                existing_cache = json.loads(_cache_path.read_text()) if _cache_path.exists() else {}
                completed_list = existing_cache.get("completed", [])
                if step_name not in completed_list:
                    completed_list.append(step_name)
                existing_cache["completed"] = completed_list
                _cache_path.write_text(json.dumps(existing_cache, indent=2))
            except Exception as exc:
                log.warning("Could not write step cache file: %s", exc)

    except Exception as exc:
        log.error("Pipeline FAILED at %s: %s", step_name, exc, exc_info=True)
        _mark_step(state, step_name, "failed")
        state["status"] = "failed"
        state["failed_at"] = _now_iso()
        state["error"] = str(exc)
        _write_state(state)
        run_id = state.get("run_id")
        if run_id:
            try:
                from db.models import PipelineRun
                from db.session import get_session
                with get_session() as session:
                    run = session.get(PipelineRun, run_id)
                    if run:
                        run.status = "failed"
                        run.finished_at = datetime.now(timezone.utc)
                        run.step_statuses = state.get("steps", {})
                        session.commit()
            except Exception as exc:
                log.warning("Could not update PipelineRun %s status to 'failed' in DB: %s", run_id, exc)
        sys.exit(1)

    total = time.time() - start_time
    log.info("")
    log.info("=" * 60)
    log.info("Phase 1 pipeline complete in %.1f minutes.", total / 60)
    log.info("=" * 60)
    state["status"] = "complete"
    state["completed_at"] = _now_iso()
    state["elapsed_total_s"] = round(total, 1)
    _write_state(state)
    run_id = state.get("run_id")
    if run_id:
        try:
            from db.models import PipelineRun
            from db.session import get_session
            with get_session() as session:
                run = session.get(PipelineRun, run_id)
                if run:
                    run.status = "completed"
                    run.finished_at = datetime.now(timezone.utc)
                    run.step_statuses = state.get("steps", {})
                    session.commit()
        except Exception as exc:
            log.warning("Could not update PipelineRun %s status to 'completed' in DB: %s", run_id, exc)


def main() -> None:
    parser = argparse.ArgumentParser(description="BioResilient AI Phase 1 Pipeline")
    parser.add_argument(
        "--resume-from",
        default="step1",
        help="Resume from a specific step (default: step1)",
    )
    parser.add_argument(
        "--steps",
        help="Comma-separated list of steps to run (overrides --resume-from)",
    )
    parser.add_argument(
        "--phenotype",
        default="",
        help=(
            "Filter species registry to this phenotype (e.g. cancer_resistance). "
            "Human baseline and control species are always included. "
            "Also used as trait_id for CandidateScore records."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run through step logic without executing tools or DB writes",
    )
    args = parser.parse_args()

    known_steps = set(STEPS.keys())

    # Validate --resume-from before starting the pipeline
    if args.resume_from and args.resume_from not in known_steps:
        parser.error(
            f"Unknown step '{args.resume_from}' for --resume-from. "
            f"Valid steps: {', '.join(STEPS.keys())}"
        )

    # Validate --steps entries
    only_steps: list[str] | None = None
    if args.steps:
        only_steps = [s.strip() for s in args.steps.split(",")]
        invalid = [s for s in only_steps if s not in known_steps]
        if invalid:
            parser.error(
                f"Unknown step(s) in --steps: {', '.join(invalid)}. "
                f"Valid steps: {', '.join(STEPS.keys())}"
            )

    trait_id = args.phenotype or ""
    species_list = _load_species_registry(phenotype=trait_id or None)

    run_pipeline(
        species_list=species_list,
        resume_from=args.resume_from,
        only_steps=only_steps,
        dry_run=args.dry_run,
        trait_id=trait_id,
    )


if __name__ == "__main__":
    main()
