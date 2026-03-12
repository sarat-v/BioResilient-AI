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
import sys
import time
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
    "step3d":  "Phylogenetic conservation scoring (phyloP/PhastCons)",
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
    "step8":   "Expression analysis (GEO/DESeq2)",
    "step8b":  "Bgee cross-species expression supplement",
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
    except Exception:
        pass


def _read_state() -> dict:
    try:
        if _STATE_FILE.exists():
            return json.loads(_STATE_FILE.read_text())
    except Exception:
        pass
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
        except Exception:
            pass

# Add repo root to sys.path so imports work from any working directory
sys.path.insert(0, str(_REPO_ROOT))


def _load_species_registry(phenotype: str | None = None) -> list[dict]:
    registry_path = _REPO_ROOT / "config" / "species_registry.json"
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
    """Validate that all required tools are installed and the DB is reachable."""
    import subprocess

    log.info("Step 1: Validating environment...")
    # Each entry is (canonical_name, [accepted_binary_names])
    tools = [
        ("orthofinder", ["orthofinder"]),
        ("mafft",       ["mafft"]),
        ("iqtree2",     ["iqtree2", "iqtree"]),   # bioconda installs as either name
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
            log.error("  ✗ %s not found", name)
        else:
            log.info("  ✓ %s", name)

    if missing:
        raise RuntimeError(f"Missing tools: {missing}. Run setup_local.sh first.")

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

    from pipeline.config import get_storage_root
    from pipeline.layer1_sequence.orthofinder import run_orthofinder

    proteomes_dir = Path(get_storage_root()) / "proteomes"
    results_dir = run_orthofinder(proteomes_dir)
    log.info("  OrthoFinder results at: %s", results_dir)
    return results_dir


def step3b_load_orthologs(results_dir: Path, dry_run: bool = False) -> None:
    """Parse OrthoFinder output and load into DB."""
    log.info("Step 3b: Loading orthologs into database...")
    if dry_run:
        return

    from pipeline.config import get_storage_root
    from pipeline.layer1_sequence.orthofinder import (
        flag_one_to_one_orthogroups,
        load_orthologs_to_db,
        load_sequence_map,
        parse_orthogroups,
    )

    proteomes_dir = Path(get_storage_root()) / "proteomes"
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

    from pipeline.config import get_storage_root
    from pipeline.layer2_evolution.phylo_conservation import run_phylo_conservation
    from db.models import CandidateScore
    from db.session import get_session as _get_session

    treefile = Path(get_storage_root()) / "phylo" / "species.treefile"

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
    except Exception:
        pass  # step9 not yet run — score all genes

    n_scored = run_phylo_conservation(tier_gene_ids, treefile)
    log.info("  Phylogenetic conservation: %d genes scored", n_scored)


def _build_human_gene_map(human_sequences: dict[str, str]) -> dict[str, tuple[str, str]]:
    """Build {protein_id: (gene_uuid, gene_symbol)} for human proteins.

    In Phase 1 we use the UniProt accession as both gene_uuid and gene_symbol
    until we enrich with NCBI Gene IDs in Phase 2.
    """
    import uuid
    gene_map = {}
    for protein_id in human_sequences:
        # UniProt accession is like "human|P12345"
        acc = protein_id.split("|")[-1]
        # Use deterministic UUID from accession
        gene_uuid = str(uuid.uuid5(uuid.NAMESPACE_DNS, f"human.{acc}"))
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

    from pipeline.config import get_storage_root

    from pipeline.config import get_thresholds
    from pipeline.layer1_sequence.alignment import (
        align_all_orthogroups,
        filter_orthogroups_by_global_identity,
        load_orthogroup_sequences_from_db,
    )
    from pipeline.layer1_sequence.divergence import (
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
    species_to_lineage = {s["id"]: s.get("lineage_group", "Other") for s in species_list}
    min_lineages = int(thresholds.get("convergence_min_lineages", 2))
    motifs_by_og = filter_by_independent_lineages(
        motifs_by_og_raw,
        species_to_lineage,
        min_lineages=min_lineages,
    )
    for og_id, motifs in motifs_by_og.items():
        load_motifs_to_db(motifs, {})

    aligned_for_hyphy = {og: aligned[og] for og in motifs_by_og}
    log.info("  After Gate 2: %d orthogroups for HyPhy.", len(aligned_for_hyphy))

    # Save aligned orthogroups to disk for resume capability
    import pickle
    from pipeline.config import sync_to_s3
    _cache_path = Path(get_storage_root()) / "aligned_orthogroups.pkl"
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

    from pipeline.config import get_storage_root, sync_to_s3
    from pipeline.layer2_evolution.phylo_tree import build_concatenated_alignment, run_iqtree

    concat = build_concatenated_alignment(aligned_orthogroups)
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

    selection_results = run_meme_pipeline(aligned_orthogroups, motifs_by_og, treefile)
    gene_by_og = build_gene_og_map()
    load_selection_scores(selection_results, gene_by_og)


def step6b_fel_busted(dry_run: bool = False) -> None:
    """Run FEL (pervasive selection) and BUSTED (gene-level episodic) supplementary tests.

    Reuses codon alignments from step6. Skips gracefully if alignments are missing.
    Updates EvolutionScore.fel_sites and EvolutionScore.busted_pvalue.
    """
    log.info("Step 6b: FEL + BUSTED supplementary selection tests...")
    if dry_run:
        return

    from pipeline.config import get_storage_root
    from pipeline.layer2_evolution.meme_selection import run_fel_busted_pipeline
    from pipeline.layer2_evolution.selection import build_gene_og_map, load_fel_busted_scores

    import pickle
    _cache_path = Path(get_storage_root()) / "aligned_orthogroups.pkl"
    if not _cache_path.exists():
        log.warning("  Aligned orthogroups cache not found; skipping step 6b.")
        return

    with open(_cache_path, "rb") as _f:
        cached = pickle.load(_f)

    aligned = cached.get("aligned", {})
    motifs_by_og = cached.get("motifs_by_og", {})

    treefile = Path(get_storage_root()) / "phylo" / "species.treefile"
    if not treefile.exists():
        log.warning("  Species treefile not found; skipping step 6b.")
        return

    results = run_fel_busted_pipeline(aligned, motifs_by_og, treefile)
    gene_by_og = build_gene_og_map()
    updated = load_fel_busted_scores(results, gene_by_og)
    log.info("  FEL/BUSTED: updated %d genes.", updated)


def step6c_relax(dry_run: bool = False) -> None:
    """Run RELAX to detect branch-specific rate acceleration in resilient species.

    Reuses codon alignments from step6. Skips gracefully if alignments are missing.
    Updates EvolutionScore.relax_k and EvolutionScore.relax_pvalue.
    """
    log.info("Step 6c: RELAX branch acceleration tests...")
    if dry_run:
        return

    from pipeline.config import get_storage_root
    from pipeline.layer2_evolution.meme_selection import run_relax_pipeline
    from pipeline.layer2_evolution.selection import build_gene_og_map, load_relax_scores

    import pickle
    _cache_path = Path(get_storage_root()) / "aligned_orthogroups.pkl"
    if not _cache_path.exists():
        log.warning("  Aligned orthogroups cache not found; skipping step 6c.")
        return

    with open(_cache_path, "rb") as f:
        cached = pickle.load(f)

    aligned = cached.get("aligned", {})
    motifs_by_og = cached.get("motifs_by_og", {})

    treefile = Path(get_storage_root()) / "phylo" / "species.treefile"
    if not treefile.exists():
        log.warning("  Species treefile not found; skipping step 6c.")
        return

    results = run_relax_pipeline(aligned, motifs_by_og, treefile)
    gene_by_og = build_gene_og_map()
    updated = load_relax_scores(results, gene_by_og)
    log.info("  RELAX: updated %d genes.", updated)


def step7_convergence(dry_run: bool = False) -> None:
    log.info("Step 7: Computing convergence detection + PhyloP enrichment...")
    if dry_run:
        return

    from pipeline.layer2_evolution.convergence import (
        build_chrom_map,
        enrich_phylop_scores,
        run_convergence_pipeline,
    )

    run_convergence_pipeline()
    chrom_map = build_chrom_map()
    if chrom_map:
        enrich_phylop_scores(chrom_map)


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


def step8_expression(species_list: list[dict], dry_run: bool = False) -> None:
    """Fetch GEO datasets and run DESeq2 expression analysis."""
    log.info("Step 8: Running expression annotation (GEO + DESeq2)...")
    if dry_run:
        return

    from pipeline.layer1_sequence.expression import run_expression_pipeline, save_expression_scores
    scores_by_species = run_expression_pipeline(species_list)
    save_expression_scores(scores_by_species)


def step8b_bgee(dry_run: bool = False) -> None:
    """Supplement expression data with Bgee curated cross-species calls."""
    log.info("Step 8b: Bgee cross-species expression supplement...")
    if dry_run:
        return

    from pipeline.layer1_sequence.bgee import run_bgee_pipeline
    n = run_bgee_pipeline()
    log.info("  Bgee: %d expression rows added.", n)


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
    """Return gene IDs for Tier1 and Tier2 candidates (for funnel gating of Phase 2 layers)."""
    from db.models import CandidateScore
    from db.session import get_session
    with get_session() as session:
        q = session.query(CandidateScore.gene_id).filter(
            CandidateScore.tier.in_(["Tier1", "Tier2"])
        )
        if trait_id:
            q = q.filter(CandidateScore.trait_id == trait_id)
        rows = q.all()
    return [r[0] for r in rows]


def step11_disease_annotation(dry_run: bool = False, trait_id: str = "") -> None:
    """Layer 3: OpenTargets, GWAS, gnomAD, IMPC, Human Protein Atlas (Tier1+Tier2 only)."""
    log.info("Step 11: Disease annotation (Layer 3)...")
    if dry_run:
        return
    gene_ids = _get_tier12_gene_ids(trait_id)
    if not gene_ids:
        log.warning("  No Tier1/Tier2 candidates; skipping disease annotation.")
        return
    from pipeline.layer3_disease import opentargets, gwas, gnomad, impc, protein_atlas, pathways
    opentargets.annotate_genes_opentargets(gene_ids)
    gwas.annotate_genes_gwas(gene_ids)
    gnomad.annotate_genes_gnomad(gene_ids)
    impc.annotate_genes_impc(gene_ids)
    protein_atlas.annotate_genes_protein_atlas(gene_ids)
    n_pathways = pathways.annotate_genes_pathways(gene_ids)
    log.info("  Pathways/GO annotated for %d genes.", n_pathways)


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
    "step4":   step4_alignment_and_divergence,
    "step4b":  step4b_domain_and_consequence,
    "step4c":  step4c_esm1v,
    "step4d":  step4d_variant_direction,
    "step5":   step5_phylogenetic_tree,
    "step6":   step6_evolutionary_selection,
    "step6b":  step6b_fel_busted,
    "step6c":  step6c_relax,
    "step7":   step7_convergence,
    "step7b":  step7b_convergent_aa,
    "step8":   step8_expression,
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
    from pipeline.config import get_storage_root
    _proteomes_dir = Path(get_storage_root()) / "proteomes"
    _existing = _orthofinder_output_dir(_proteomes_dir)
    if _existing is not None:
        results_dir = _existing

    # If resuming past step4, recover aligned_orthogroups from disk cache
    if resume_from not in ("step1","step2","step3","step3b","step4"):
        _pkl_path = Path(get_storage_root()) / "aligned_orthogroups.pkl"
        if _pkl_path.exists():
            try:
                import pickle
                _data = pickle.load(open(_pkl_path, "rb"))
                aligned_orthogroups = _data["aligned"]
                motifs_by_og = _data["motifs_by_og"]
                log.info("Recovered %d aligned orthogroups from disk cache", len(aligned_orthogroups))
            except Exception as _e:
                log.warning("Could not load aligned orthogroups cache: %s", _e)
        else:
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

    # If resuming past step5, recover treefile from disk
    _tree_candidate = Path(get_storage_root()) / "phylo" / "species.treefile"
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
    import os
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
    except Exception:
        pass
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
                if treefile and motifs_by_og:
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
                step8_expression(species_list, dry_run=dry_run)

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
            except Exception:
                pass

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
            except Exception:
                pass
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
        except Exception:
            pass


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

    trait_id = args.phenotype or ""
    species_list = _load_species_registry(phenotype=trait_id or None)
    only_steps = [s.strip() for s in args.steps.split(",")] if args.steps else None

    run_pipeline(
        species_list=species_list,
        resume_from=args.resume_from,
        only_steps=only_steps,
        dry_run=args.dry_run,
        trait_id=trait_id,
    )


if __name__ == "__main__":
    main()
