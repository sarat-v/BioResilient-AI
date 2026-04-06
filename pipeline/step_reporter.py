"""BioResilient AI — Per-step reporter and scientific plausibility validator.

Called by run_cancer_resistance_stepwise.sh after every pipeline step.

Produces two files per step in step_cache/:
  <step>.json  — machine-readable structured data (all counts, tables, distributions)
  <step>.md    — human-readable Markdown report with explanation, numbers, and flags

Also provides validate(step) → ValidationResult with PASS / WARN / FAIL status
and per-check detail, used by the stepwise runner to drive the go/stop prompt.

Usage:
    python pipeline/step_reporter.py --step step4
    python pipeline/step_reporter.py --step step4 --show      # pretty-print .md to terminal
    python pipeline/step_reporter.py --step step9 --validate  # print validation result only
"""

import argparse
import json
import logging
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parents[1]
_CACHE_DIR = _REPO_ROOT / "step_cache"

# ── Vertebrate species that must have ≥10k proteins to be acceptable ─────────
_VERTEBRATE_SPECIES = {
    "naked_mole_rat", "blind_mole_rat", "damaraland_mole_rat", "beaver",
    "bowhead_whale", "sperm_whale", "african_elephant", "asian_elephant",
    "greenland_shark", "little_skate", "painted_turtle", "little_brown_bat",
    "elephant_shark", "human", "rat", "macaque",
}

# ── Benchmark cancer-resistance genes used in step9 gate ─────────────────────
_CANCER_BENCHMARKS = {
    "P53_HUMAN", "ATM_HUMAN", "BRCA1_HUMAN", "BRCA2_HUMAN",
    "CDN2A_HUMAN", "MDM2_HUMAN", "PTEN_HUMAN", "PCNA_HUMAN", "ERCC1_HUMAN",
    "CASP8_HUMAN", "CASP3_HUMAN", "BCL2_HUMAN", "BAX_HUMAN", "CHK2_HUMAN",
    "TP53", "LIF6", "ATM", "BRCA1", "BRCA2", "CDKN2A", "MDM2", "PTEN", "PCNA", "ERCC1",
}

# ── Step explanation text (mirrors the bash explain() blocks) ────────────────
_EXPLANATIONS: dict[str, str] = {
    "step1": (
        "Validates that every required bioinformatics tool is installed (OrthoFinder, MAFFT, "
        "IQ-TREE2, HyPhy, DIAMOND), the AWS RDS PostgreSQL database is reachable, and optional "
        "GPU acceleration is detected. Nothing analytical — purely a safety check."
    ),
    "step2": (
        "Downloads reference proteomes (all annotated proteins) from NCBI RefSeq for all 18 species: "
        "15 cancer-resistant foreground species, 1 human baseline, 2 controls (rat, macaque). "
        "Headers are rewritten to a standardised format OrthoFinder requires. "
        "Proteome quality (protein count, file size) determines the ceiling for all downstream steps."
    ),
    "step3": (
        "OrthoFinder clusters all proteins across 18 species using all-vs-all DIAMOND BLAST + MCL graph "
        "clustering into orthogroups (OGs) — groups of proteins sharing a common ancestor. "
        "Only OGs with a human member are kept. Single-copy OGs (1:1 orthologs) are flagged for the "
        "species tree in step 5."
    ),
    "step3b": (
        "Parses OrthoFinder output and loads every ortholog relationship into the database. "
        "Gene table gets one row per human protein; Ortholog table gets one row per (gene × species)."
    ),
    "step3c": (
        "Nucleotide conservation layer (DNA-level, supplementary to protein analysis). "
        "For each gene × species pair, downloads the genome GFF3 annotation and genomic FASTA from NCBI, "
        "extracts three sequence windows: CDS (coding sequence), promoter (5 kb upstream of TSS), "
        "and downstream (2 kb after gene end). Aligns each resilient species against human using "
        "minimap2 -x asm5 (asm-to-asm preset). Computes conservation_score = "
        "(pct_identity × alignment_length) / region_length. "
        "For promoter regions, counts regulatory_divergence_count (mutations shared by ≥3 resilient species, "
        "absent in human + controls) and regulatory_convergence_count (same mutation in ≥2 distinct "
        "evolutionary lineage clusters — e.g. whale + bat but not controls). "
        "This captures regulatory evolution that protein alignments cannot detect."
    ),
    "step3d": (
        "Phylogenetic conservation scoring with phyloP or PhastCons (PHAST package). "
        "Builds a multi-species FASTA alignment for each gene region from NucleotideRegion rows, "
        "then runs phyloP / PhastCons using the IQ-TREE2 species tree from step 5. "
        "phyloP scores: positive = strongly conserved site; negative = accelerated evolution "
        "(common signal of adaptive change). Restricted to Tier1+Tier2 genes after step 9 "
        "for performance. NULL scores are written when the tools are not installed (acceptable; "
        "the nucleotide_divergence_score falls back to minimap2-derived signal alone)."
    ),
    "step4": (
        "Gate 1: removes OGs where all resilient species are >85% identical to human (no divergence to find). "
        "MAFFT alignment: aligns surviving OGs column-by-column. "
        "Sliding-window divergence: flags 10–20 aa windows that differ significantly between resilient and human. "
        "Gate 2: keeps only OGs where divergence appears independently in ≥3 distinct evolutionary lineages "
        "(random drift would not cause the same region to change in a shark, bat, AND elephant)."
    ),
    "step4b": (
        "Pfam/InterPro: annotates which protein domain each divergent motif falls in. "
        "AlphaMissense: predicts whether each substitution is likely pathogenic (disrupts function) "
        "or benign. High AM score in a resilient species = possible gain-of-function or domain repurposing."
    ),
    "step4c": (
        "ESM-1v (Meta's protein language model, 650M params) computes a log-likelihood ratio (LLR) "
        "for each substitution. Negative LLR = unusual change given evolutionary context = potentially adaptive. "
        "Complements AlphaMissense with an evolutionary fitness lens."
    ),
    "step4d": (
        "Classifies each divergent motif as GoF (gain-of-function), LoF (loss-of-function), or neutral "
        "using a voting rule across AlphaMissense score, ESM-1v LLR, and gnomAD LOEUF constraint. "
        "For cancer resistance: GoF in tumour suppressors and LoF in growth signalling are the expected patterns."
    ),
    "step5": (
        "IQ-TREE2 builds a maximum-likelihood species tree from a concatenated alignment of ~100–300 "
        "single-copy orthogroups. Uses the best-fit substitution model (ModelFinder) + 1000 ultrafast "
        "bootstrap replicates. The tree is essential for all convergence and selection tests — a wrong "
        "topology would conflate inherited changes with convergent ones."
    ),
    "step6": (
        "HyPhy MEME tests each codon site in every candidate orthogroup for episodic positive selection "
        "(dN/dS > 1 in a subset of branches). Returns a p-value per site. Genes where divergent motifs "
        "overlap MEME-positive sites have statistical evidence the change was actively selected for."
    ),
    "step6b": (
        "FEL detects pervasive (site-wide) selection across all branches — complementary to MEME's episodic test. "
        "BUSTED is a gene-level test asking whether any sites anywhere show positive selection. "
        "Together with MEME they give three statistically orthogonal lines of selection evidence."
    ),
    "step6c": (
        "RELAX asks whether selective constraints are intensified (k>1) or relaxed (k<1) in resilient species "
        "vs background. k>1 in DNA repair genes = enhanced purifying selection (consistent with elephant's "
        "extra TP53). k<1 in growth signalling = relaxed constraint on the pathway."
    ),
    "step7": (
        "Counts how many independent evolutionary lineage groups show the same divergent motif. "
        "convergence_count is the core signal: 3+ independent lineages achieving the same change is "
        "statistically very unlikely by chance. PhyloP enrichment rewards sites that are highly conserved "
        "in most vertebrates but specifically changed in resilient species."
    ),
    "step7b": (
        "Checks whether independent lineages evolved the SAME amino acid substitution (or biochemically "
        "equivalent one, per BLOSUM62). Same-substitution convergence is much stronger evidence than "
        "merely 'both regions diverged'."
    ),
    "step8": (
        "Searches NCBI GEO for RNA-seq datasets per species, downloads count matrices, runs DESeq2 to "
        "find differentially expressed genes. Expression evidence validates that a protein with an "
        "adaptive sequence change is also functionally active in the relevant tissue."
    ),
    "step8b": (
        "Supplements GEO data with Bgee pre-curated cross-species expression calls for species with "
        "limited GEO coverage (painted turtle, ocean quahog clam, Hydra)."
    ),
    "step9": (
        "Phase 1 composite scoring: convergence×0.40 + selection×0.35 + expression×0.25. "
        "Tier1 ≥0.70, Tier2 ≥0.40. This is the PRIMARY CHECKPOINT — the first biological preview "
        "of which genes the evolutionary evidence puts at the top. If known benchmark genes "
        "(TP53, ATM, ERCC1) are not in Tier1/2, the pipeline needs investigation before Phase 2."
    ),
    "step10": "API ready. No computation — FastAPI server can now serve Phase 1 results.",
    "step10b": (
        "AlphaGenome (Track B): even if the protein sequence is similar, the gene's regulatory region "
        "may be divergent. AlphaGenome predicts chromatin/TF-binding from DNA sequence to detect "
        "promoter-level adaptation (especially relevant for NMR p16/p27 upregulation)."
    ),
    "step11": (
        "Phase 2 disease annotation (Tier1+Tier2 only). OpenTargets: cancer association scores. "
        "GWAS Catalog: cancer predisposition hits. gnomAD: pLI constraint. IMPC: mouse KO phenotypes. "
        "Human Protein Atlas: tumour tissue expression."
    ),
    "step11b": (
        "PCSK9 paradigm search: maps divergent motif positions to rare human gnomAD variants. "
        "If humans with rare variants at the same positions the resilient species diverged at are "
        "protected from cancer, that is the most directly translatable signal in the pipeline."
    ),
    "step11c": "PubMed sanity check: scores genes by co-mentions with cancer resistance / longevity.",
    "step11d": (
        "Hypergeometric pathway enrichment over Reactome + GO. Identifies which biological pathways "
        "are over-represented among the top candidates — e.g. DNA repair, cell cycle, apoptosis."
    ),
    "step12": (
        "Druggability: AlphaFold structures + fpocket pocket detection + ChEMBL existing drugs + "
        "CanSAR tractability tiers + peptide surface accessibility."
    ),
    "step12b": (
        "P2Rank ML pocket prediction (random forest on point clouds). More accurate than fpocket on "
        "AlphaFold models. If a divergent motif residue lines the predicted pocket, the adaptive "
        "change is directly in the binding site."
    ),
    "step13": (
        "Gene therapy feasibility (Tier1 only). AAV compatibility (gene size, tissue tropism). "
        "CRISPR guide efficiency and off-target risk for each candidate."
    ),
    "step14": (
        "Safety pre-screen: PheWAS disease associations, STRING network hub risk (degree > 100), "
        "protein family off-target selectivity."
    ),
    "step14b": (
        "DepMap CRISPR essentiality (broadly essential genes are risky to inhibit). "
        "GTEx expression breadth (narrow expression = safer therapeutic window)."
    ),
    "step15": (
        "Final rescore with Phase 2 weights: convergence×0.22 + selection×0.18 + disease×0.20 + "
        "druggability×0.15 + expression×0.10 + safety×0.10 + regulatory×0.05. "
        "Control divergence penalty removes signals also seen in rat/macaque. "
        "'Validated' tier = Tier1/2 + human genetics evidence."
    ),
    "step16": "Pipeline complete.",
}


# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class CheckResult:
    name: str
    result: str          # "PASS" | "WARN" | "FAIL" | "INFO"
    value: Any
    threshold: str
    message: str = ""


@dataclass
class ValidationResult:
    step: str
    status: str          # "PASS" | "WARN" | "FAIL"
    checks: list[CheckResult] = field(default_factory=list)
    recommendation: str = ""

    def add(self, name: str, result: str, value: Any, threshold: str, message: str = "") -> None:
        self.checks.append(CheckResult(name, result, value, threshold, message))
        if result == "FAIL" and self.status != "FAIL":
            self.status = "FAIL"
        elif result == "WARN" and self.status == "PASS":
            self.status = "WARN"

    def to_dict(self) -> dict:
        d = asdict(self)
        return d


# ─────────────────────────────────────────────────────────────────────────────
# DB helpers
# ─────────────────────────────────────────────────────────────────────────────

def _db_session():
    sys.path.insert(0, str(_REPO_ROOT))
    from db.session import get_session
    return get_session


def _scalar(query_fn):
    """Run a callable that returns a scalar and swallow errors."""
    try:
        return query_fn()
    except Exception:
        return None


# ─────────────────────────────────────────────────────────────────────────────
# Per-step data collection
# ─────────────────────────────────────────────────────────────────────────────

def _collect_step1() -> dict:
    import subprocess
    import time

    data: dict[str, Any] = {"step": "step1", "timestamp": _now()}
    # Singleton tools — must be present
    single_tools = ["orthofinder", "mafft", "diamond", "fpocket"]
    # Aliased tools — at least one name in each pair must be present
    aliased_tools = [("iqtree2", "iqtree"), ("hyphy", "HYPHYMP")]

    def _which(t: str) -> bool:
        return subprocess.run(["which", t], capture_output=True).returncode == 0

    found, missing = [], []
    for t in single_tools:
        (found if _which(t) else missing).append(t)
    for aliases in aliased_tools:
        hit = next((a for a in aliases if _which(a)), None)
        if hit:
            found.append(hit)
        else:
            missing.append(aliases[0])  # report canonical name only

    data["tools_found"] = found
    data["tools_missing"] = missing

    # DB ping
    try:
        from db.session import get_engine
        from sqlalchemy import text
        t0 = time.time()
        with get_engine().connect() as c:
            c.execute(text("SELECT 1"))
        data["db_ping_ms"] = round((time.time() - t0) * 1000, 1)
        data["db_ok"] = True
    except Exception as e:
        data["db_ok"] = False
        data["db_error"] = str(e)

    # GPU
    try:
        import torch
        data["gpu_available"] = torch.cuda.is_available()
        data["gpu_name"] = torch.cuda.get_device_name(0) if torch.cuda.is_available() else None
    except ImportError:
        data["gpu_available"] = False
        data["gpu_name"] = None

    return data


def _collect_step2() -> dict:
    from db.session import get_session
    from db.models import Species
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step2", "timestamp": _now()}
    with get_session() as s:
        species = s.query(Species).all()
        data["species_count"] = len(species)
        data["species"] = [
            {
                "id": sp.id,
                "proteome_path": sp.proteome_path,
                "is_control": sp.is_control,
                "lineage_group": sp.lineage_group,
            }
            for sp in species
        ]

    # Count proteins per FASTA
    from pipeline.config import get_local_storage_root, get_deployment
    storage = get_local_storage_root()
    if storage.startswith("s3://"):
        proteomes_dir = Path("/tmp/bioresilient/proteomes")
    else:
        proteomes_dir = Path(storage) / "proteomes"

    protein_counts: dict[str, int] = {}
    if proteomes_dir.exists():
        for faa in proteomes_dir.glob("*.reheadered.faa"):
            sid = faa.name.replace(".reheadered.faa", "")
            try:
                from Bio import SeqIO
                protein_counts[sid] = sum(1 for _ in SeqIO.parse(str(faa), "fasta"))
            except Exception:
                protein_counts[sid] = -1
    data["protein_counts"] = protein_counts
    return data


def _collect_step3() -> dict:
    from db.session import get_session
    from db.models import Gene, Ortholog
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step3", "timestamp": _now()}
    with get_session() as s:
        data["gene_count"] = s.query(func.count(Gene.id)).scalar() or 0
        data["ortholog_count"] = s.query(func.count(Ortholog.id)).scalar() or 0
        data["one_to_one_count"] = (
            s.query(func.count(Ortholog.id)).filter(Ortholog.is_one_to_one == True).scalar() or 0
        )
        # Species coverage in orthogroups
        coverage = dict(
            s.query(Ortholog.species_id, func.count(Ortholog.gene_id))
            .group_by(Ortholog.species_id)
            .all()
        )
        data["species_ortholog_coverage"] = coverage

        # Top 10 largest orthogroups
        top_ogs = (
            s.query(Ortholog.orthofinder_og, func.count(Ortholog.id))
            .group_by(Ortholog.orthofinder_og)
            .order_by(func.count(Ortholog.id).desc())
            .limit(10)
            .all()
        )
        data["top10_orthogroups"] = [{"og": og, "members": n} for og, n in top_ogs]
    return data


def _collect_step3c() -> dict:
    from db.session import get_session
    from db.models import NucleotideRegion, NucleotideScore, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step3c", "timestamp": _now()}
    with get_session() as s:
        data["nucleotide_regions_total"] = s.query(func.count(NucleotideRegion.id)).scalar() or 0
        data["genes_with_regions"] = (
            s.query(func.count(func.distinct(NucleotideRegion.gene_id))).scalar() or 0
        )
        data["species_with_regions"] = (
            s.query(func.count(func.distinct(NucleotideRegion.species_id))).scalar() or 0
        )

        # Breakdown by region_type
        type_counts = dict(
            s.query(NucleotideRegion.region_type, func.count(NucleotideRegion.id))
            .group_by(NucleotideRegion.region_type)
            .all()
        )
        data["regions_by_type"] = type_counts

        # NucleotideScore stats
        data["genes_scored"] = (
            s.query(func.count(func.distinct(NucleotideScore.gene_id))).scalar() or 0
        )
        # Promoter divergence
        total_div = s.query(func.sum(NucleotideScore.regulatory_divergence_count)).scalar() or 0
        total_conv = s.query(func.sum(NucleotideScore.regulatory_convergence_count)).scalar() or 0
        data["total_regulatory_divergence_events"] = int(total_div)
        data["total_regulatory_convergence_events"] = int(total_conv)

        # Mean conservation by region type
        for rtype in ("cds", "promoter", "downstream"):
            avg = (
                s.query(func.avg(NucleotideScore.conservation_score))
                .filter(NucleotideScore.region_type == rtype,
                        NucleotideScore.conservation_score != None)
                .scalar()
            )
            data[f"mean_{rtype}_conservation"] = round(float(avg), 3) if avg is not None else None

        # Top 10 genes by promoter divergence_count
        top_div = (
            s.query(NucleotideScore.gene_id, NucleotideScore.regulatory_divergence_count)
            .filter(NucleotideScore.region_type == "promoter",
                    NucleotideScore.regulatory_divergence_count != None)
            .order_by(NucleotideScore.regulatory_divergence_count.desc())
            .limit(10)
            .all()
        )
        gene_ids = [r[0] for r in top_div]
        sym = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())
        data["top10_promoter_divergence"] = [
            {"gene": sym.get(gid, gid), "divergence_count": n} for gid, n in top_div
        ]
    return data


def _collect_step3d() -> dict:
    from db.session import get_session
    from db.models import PhyloConservationScore, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step3d", "timestamp": _now()}
    with get_session() as s:
        data["genes_with_phylo_score"] = s.query(func.count(PhyloConservationScore.gene_id)).scalar() or 0

        for field in ("cds_phylo_score", "promoter_phylo_score", "downstream_phylo_score"):
            attr = getattr(PhyloConservationScore, field)
            avg = s.query(func.avg(attr)).filter(attr != None).scalar()
            data[f"mean_{field}"] = round(float(avg), 4) if avg is not None else None

        # Genes with accelerated promoter evolution (negative phyloP)
        data["genes_accelerated_promoter"] = (
            s.query(func.count(PhyloConservationScore.gene_id))
            .filter(PhyloConservationScore.promoter_phylo_score < 0)
            .scalar() or 0
        )
        # Top 10 by most accelerated promoter (most negative score)
        top_acc = (
            s.query(PhyloConservationScore.gene_id, PhyloConservationScore.promoter_phylo_score)
            .filter(PhyloConservationScore.promoter_phylo_score != None)
            .order_by(PhyloConservationScore.promoter_phylo_score.asc())
            .limit(10)
            .all()
        )
        gene_ids = [r[0] for r in top_acc]
        sym = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())
        data["top10_accelerated_promoter"] = [
            {"gene": sym.get(gid, gid), "promoter_phylo_score": round(float(sc), 4) if sc else None}
            for gid, sc in top_acc
        ]
    return data


def _collect_step4() -> dict:
    from db.session import get_session
    from db.models import DivergentMotif, Ortholog, Gene, Species
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step4", "timestamp": _now()}
    with get_session() as s:
        data["motif_count"] = s.query(func.count(DivergentMotif.id)).scalar() or 0

        # genes_with_motifs: count distinct genes via ortholog join
        data["genes_with_motifs"] = (
            s.query(func.count(func.distinct(Ortholog.gene_id)))
            .join(DivergentMotif, DivergentMotif.ortholog_id == Ortholog.id)
            .scalar() or 0
        )

        # Per-species motif coverage via ortholog join
        coverage_q = (
            s.query(Ortholog.species_id, func.count(DivergentMotif.id))
            .join(DivergentMotif, DivergentMotif.ortholog_id == Ortholog.id)
            .group_by(Ortholog.species_id)
            .all()
        )
        data["motifs_per_species"] = dict(coverage_q)

        # Lineage coverage — distinct genes per lineage
        lineage_q = (
            s.query(Species.lineage_group, func.count(func.distinct(Ortholog.gene_id)))
            .join(Ortholog, Ortholog.species_id == Species.id)
            .join(DivergentMotif, DivergentMotif.ortholog_id == Ortholog.id)
            .group_by(Species.lineage_group)
            .all()
        )
        data["lineage_motif_coverage"] = dict(lineage_q)

        # Top 20 genes ranked by distinct species showing divergence (most broadly diverged genes)
        top_genes = (
            s.query(
                Ortholog.gene_id,
                func.count(DivergentMotif.id).label("motif_count"),
                func.count(func.distinct(Ortholog.species_id)).label("species_count"),
                func.avg(DivergentMotif.divergence_score).label("avg_div"),
            )
            .join(DivergentMotif, DivergentMotif.ortholog_id == Ortholog.id)
            .group_by(Ortholog.gene_id)
            .order_by(
                func.count(func.distinct(Ortholog.species_id)).desc(),
                func.count(DivergentMotif.id).desc(),
                func.avg(DivergentMotif.divergence_score).desc(),
            )
            .limit(20)
            .all()
        )
        gene_ids = [r[0] for r in top_genes]
        gene_symbols = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())
        data["top20_genes_by_motifs"] = [
            {
                "gene_id": gid,
                "symbol": gene_symbols.get(gid, gid),
                "motif_count": n,
                "species_count": sp,
                "avg_divergence": round(float(avg_div), 4) if avg_div else None,
            }
            for gid, n, sp, avg_div in top_genes
        ]
    return data


def _collect_step4b() -> dict:
    from db.session import get_session
    from db.models import DivergentMotif, Ortholog, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step4b", "timestamp": _now()}
    with get_session() as s:
        total = s.query(func.count(DivergentMotif.id)).scalar() or 0
        data["motif_count"] = total
        in_domain = (
            s.query(func.count(DivergentMotif.id))
            .filter(DivergentMotif.in_functional_domain == True).scalar() or 0
        )
        data["motifs_in_functional_domain"] = in_domain
        data["domain_pct"] = round(in_domain / total * 100, 1) if total else 0

        # Top Pfam domains
        top_domains = (
            s.query(DivergentMotif.domain_name, func.count(DivergentMotif.id))
            .filter(DivergentMotif.domain_name != None)
            .group_by(DivergentMotif.domain_name)
            .order_by(func.count(DivergentMotif.id).desc())
            .limit(10)
            .all()
        )
        data["top10_pfam_domains"] = [{"domain": d, "count": n} for d, n in top_domains]

        # AlphaMissense distribution
        am_scored = (
            s.query(func.count(DivergentMotif.id))
            .filter(DivergentMotif.consequence_score != None).scalar() or 0
        )
        data["am_scored_count"] = am_scored
        data["am_pct_coverage"] = round(am_scored / total * 100, 1) if total else 0

        # Pathogenic (AM > 0.564)
        am_pathogenic = (
            s.query(func.count(DivergentMotif.id))
            .filter(DivergentMotif.consequence_score > 0.564).scalar() or 0
        )
        data["am_likely_pathogenic"] = am_pathogenic

        # Top 10 highest AM score motifs with gene symbol
        top_am = (
            s.query(Ortholog.gene_id, DivergentMotif.consequence_score,
                    DivergentMotif.domain_name)
            .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
            .filter(DivergentMotif.consequence_score != None)
            .order_by(DivergentMotif.consequence_score.desc())
            .limit(10)
            .all()
        )
        gene_ids = [r[0] for r in top_am]
        gene_symbols = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())
        data["top10_am_motifs"] = [
            {"gene": gene_symbols.get(gid, gid), "am_score": round(float(am), 3),
             "domain": dom}
            for gid, am, dom in top_am
        ]
    return data


def _collect_step4c() -> dict:
    from db.session import get_session
    from db.models import DivergentMotif, Ortholog, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step4c", "timestamp": _now()}
    with get_session() as s:
        total = s.query(func.count(DivergentMotif.id)).scalar() or 0
        scored = (
            s.query(func.count(DivergentMotif.id))
            .filter(DivergentMotif.esm1v_score != None).scalar() or 0
        )
        data["motif_count"] = total
        data["esm1v_scored"] = scored
        data["esm1v_pct_coverage"] = round(scored / total * 100, 1) if total else 0

        # Mean LLR
        avg = s.query(func.avg(DivergentMotif.esm1v_score)).filter(
            DivergentMotif.esm1v_score != None
        ).scalar()
        data["mean_llr"] = round(float(avg), 3) if avg is not None else None

        # Histogram bins: < -4, -4 to -2, -2 to 0, > 0
        bins = {}
        for label, lo, hi in [("< -4", None, -4), ("-4 to -2", -4, -2),
                                ("-2 to 0", -2, 0), ("> 0", 0, None)]:
            q = s.query(func.count(DivergentMotif.id)).filter(
                DivergentMotif.esm1v_score != None
            )
            if lo is not None:
                q = q.filter(DivergentMotif.esm1v_score >= lo)
            if hi is not None:
                q = q.filter(DivergentMotif.esm1v_score < hi)
            bins[label] = q.scalar() or 0
        data["llr_distribution"] = bins

        # Top 10 most negative LLR (most unusual substitutions)
        top_unusual = (
            s.query(Ortholog.gene_id, DivergentMotif.esm1v_score,
                    DivergentMotif.start_pos, DivergentMotif.human_seq, DivergentMotif.animal_seq)
            .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
            .filter(DivergentMotif.esm1v_score != None)
            .order_by(DivergentMotif.esm1v_score.asc())
            .limit(10)
            .all()
        )
        gene_ids = [r[0] for r in top_unusual]
        gene_symbols = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())
        data["top10_unusual_motifs"] = [
            {"gene": gene_symbols.get(gid, gid), "llr": round(float(llr), 3),
             "pos": pos, "human": hseq, "animal": aseq}
            for gid, llr, pos, hseq, aseq in top_unusual
        ]
    return data


def _collect_step4d() -> dict:
    from db.session import get_session
    from db.models import DivergentMotif, Ortholog, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step4d", "timestamp": _now()}
    with get_session() as s:
        direction_counts = dict(
            s.query(DivergentMotif.motif_direction, func.count(DivergentMotif.id))
            .group_by(DivergentMotif.motif_direction)
            .all()
        )
        data["direction_counts"] = {k or "unclassified": v for k, v in direction_counts.items()}

        for direction in ("gain_of_function", "loss_of_function"):
            top = (
                s.query(Ortholog.gene_id, func.count(DivergentMotif.id))
                .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
                .filter(DivergentMotif.motif_direction == direction)
                .group_by(Ortholog.gene_id)
                .order_by(func.count(DivergentMotif.id).desc())
                .limit(10)
                .all()
            )
            gene_ids = [r[0] for r in top]
            sym = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())
            data[f"top10_{direction}"] = [
                {"gene": sym.get(gid, gid), "motif_count": n} for gid, n in top
            ]
        # functional_shift count (previously 'likely_pathogenic' — renamed to clarify
        # these are hypothesis labels from human-centric ML models, not hard filters)
        data["functional_shift_count"] = direction_counts.get("functional_shift", 0)
    return data


def _collect_step5() -> dict:
    from pipeline.config import get_local_storage_root
    import subprocess

    data: dict[str, Any] = {"step": "step5", "timestamp": _now()}
    storage = get_local_storage_root()

    # Check all candidate treefile locations (local or S3-backed)
    candidates = [
        Path(storage) / "phylo" / "species.treefile",
        Path("/tmp/bioresilient/phylo/species.treefile"),
        Path("/tmp/species.treefile"),
        Path("step_cache/species.treefile"),
    ]
    treefile = next((p for p in candidates if p.exists()), None)

    # Fallback: pull from S3 cache if none found locally
    if treefile is None:
        s3_paths = [
            "s3://bioresilient-data/cache/species.treefile",
            "s3://bioresilient-data/step_cache/cancer_resistance/species.treefile",
        ]
        local_tmp = Path("/tmp/species_reporter.treefile")
        for s3p in s3_paths:
            try:
                r = subprocess.run(
                    ["aws", "s3", "cp", s3p, str(local_tmp)],
                    capture_output=True, timeout=30
                )
                if r.returncode == 0 and local_tmp.exists():
                    treefile = local_tmp
                    break
            except Exception:
                pass

    data["treefile_exists"] = treefile is not None and treefile.exists()
    data["treefile_path"] = str(treefile) if treefile else str(Path(storage) / "phylo" / "species.treefile")

    if treefile is not None and treefile.exists():
        newick = treefile.read_text().strip()
        data["newick"] = newick

        # Bootstrap support stats via ete3
        try:
            from ete3 import Tree
            t = Tree(newick)
            supports = [float(n.support) for n in t.traverse() if n.support and n.support > 0]
            if supports:
                data["bootstrap_min"] = round(min(supports), 1)
                data["bootstrap_mean"] = round(sum(supports) / len(supports), 1)
                data["bootstrap_pct_above_90"] = round(
                    sum(1 for s in supports if s >= 90) / len(supports) * 100, 1
                )
            else:
                data["bootstrap_mean"] = None
        except ImportError:
            data["bootstrap_mean"] = None
            data["bootstrap_note"] = "ete3 not installed — bootstrap stats unavailable"
        except Exception as e:
            data["bootstrap_mean"] = None
            data["bootstrap_note"] = str(e)
    return data


def _collect_step6() -> dict:
    from db.session import get_session
    from db.models import EvolutionScore, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step6", "timestamp": _now()}
    with get_session() as s:
        data["genes_with_selection"] = (
            s.query(func.count(func.distinct(EvolutionScore.gene_id))).scalar() or 0
        )

        # PAML branch-site model: dnds_pvalue is the primary signal (LRT p-value)
        # meme_qvalue is the legacy HyPhy field — may be NULL for PAML runs
        paml_count = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.dnds_pvalue != None).scalar() or 0
        )
        data["paml_scored_genes"] = paml_count

        # Significant positive selection (LRT p < 0.05)
        paml_positive = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.dnds_pvalue != None, EvolutionScore.dnds_pvalue < 0.05)
            .scalar() or 0
        )
        data["paml_significant_genes"] = paml_positive

        # Strong positive selection (foreground omega > 1)
        omega_positive = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.dnds_pvalue < 0.05, EvolutionScore.dnds_ratio > 1)
            .scalar() or 0
        )
        data["omega_positive_genes"] = omega_positive

        # Mean dN/dS for significant genes
        avg_omega = (
            s.query(func.avg(EvolutionScore.dnds_ratio))
            .filter(EvolutionScore.dnds_pvalue < 0.05)
            .scalar()
        )
        data["mean_omega_significant"] = round(float(avg_omega), 4) if avg_omega is not None else None

        # Model distribution
        model_counts = dict(
            s.query(EvolutionScore.selection_model, func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.selection_model != None)
            .group_by(EvolutionScore.selection_model)
            .all()
        )
        data["model_counts"] = model_counts

        # Also check legacy HyPhy MEME field (0 for PAML runs)
        meme_positive = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.meme_qvalue != None, EvolutionScore.meme_qvalue < 0.1)
            .scalar() or 0
        )
        data["meme_positive_genes"] = meme_positive

        # Top 20 genes by PAML significance (lowest p-value)
        top = (
            s.query(EvolutionScore.gene_id, EvolutionScore.dnds_pvalue, EvolutionScore.dnds_ratio)
            .filter(EvolutionScore.dnds_pvalue != None)
            .order_by(EvolutionScore.dnds_pvalue.asc())
            .limit(20)
            .all()
        )
        gene_ids = [r[0] for r in top]
        sym = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())

        # Check for known benchmarks in top results
        top_symbols = {sym.get(gid, "") for gid, _, _ in top}
        data["benchmark_genes_in_top20"] = list(
            _CANCER_BENCHMARKS & {gs.upper() for gs in top_symbols if gs}
        )
        # Keep legacy field name for backward compat
        data["benchmark_genes_in_top20_meme"] = data["benchmark_genes_in_top20"]

        data["top20_paml_genes"] = [
            {"gene": sym.get(gid, gid), "dnds_pvalue": round(float(p), 6) if p else None,
             "dnds_ratio": round(float(d), 3) if d else None}
            for gid, p, d in top
        ]
        # Keep legacy field for backward compat
        data["top20_meme_genes"] = [
            {"gene": r["gene"], "meme_qvalue": r["dnds_pvalue"], "dnds": r["dnds_ratio"]}
            for r in data["top20_paml_genes"]
        ]
    return data


def _collect_step6b() -> dict:
    from db.session import get_session
    from db.models import EvolutionScore, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step6b", "timestamp": _now()}
    with get_session() as s:
        data["fel_positive_genes"] = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.fel_sites != None, EvolutionScore.fel_sites > 0)
            .scalar() or 0
        )
        data["busted_positive_genes"] = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.busted_pvalue != None, EvolutionScore.busted_pvalue < 0.05)
            .scalar() or 0
        )
        # Overlap: MEME + BUSTED both positive
        both = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(
                EvolutionScore.meme_qvalue < 0.1,
                EvolutionScore.busted_pvalue < 0.05,
            )
            .scalar() or 0
        )
        data["meme_and_busted_both_positive"] = both
    return data


def _collect_step6c() -> dict:
    from db.session import get_session
    from db.models import EvolutionScore, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step6c", "timestamp": _now()}
    with get_session() as s:
        data["relax_intensified"] = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.relax_k != None, EvolutionScore.relax_k > 1,
                    EvolutionScore.relax_pvalue < 0.05)
            .scalar() or 0
        )
        data["relax_relaxed"] = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.relax_k != None, EvolutionScore.relax_k < 1,
                    EvolutionScore.relax_pvalue < 0.05)
            .scalar() or 0
        )

        for label, filt in [("intensified", EvolutionScore.relax_k > 1),
                             ("relaxed", EvolutionScore.relax_k < 1)]:
            top = (
                s.query(EvolutionScore.gene_id, EvolutionScore.relax_k, EvolutionScore.relax_pvalue)
                .filter(EvolutionScore.relax_k != None, EvolutionScore.relax_pvalue < 0.05, filt)
                .order_by(EvolutionScore.relax_pvalue.asc())
                .limit(10)
                .all()
            )
            gene_ids = [r[0] for r in top]
            sym = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())
            data[f"top10_{label}"] = [
                {"gene": sym.get(gid, gid), "k": round(float(k), 3) if k else None,
                 "pvalue": round(float(p), 5) if p else None}
                for gid, k, p in top
            ]
    return data


def _collect_step7() -> dict:
    from db.session import get_session
    from db.models import EvolutionScore, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step7", "timestamp": _now()}
    with get_session() as s:
        # Distribution table
        dist: dict[str, int] = {}
        for threshold in range(1, 9):
            n = (
                s.query(func.count(EvolutionScore.gene_id))
                .filter(EvolutionScore.convergence_count == threshold)
                .scalar() or 0
            )
            if n > 0:
                dist[str(threshold)] = n
        n_plus = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.convergence_count >= 8)
            .scalar() or 0
        )
        if n_plus:
            dist["8+"] = n_plus
        data["convergence_count_distribution"] = dist

        data["genes_conv_ge3"] = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.convergence_count >= 3).scalar() or 0
        )
        data["genes_conv_ge4"] = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.convergence_count >= 4).scalar() or 0
        )
        data["genes_conv_ge5"] = (
            s.query(func.count(EvolutionScore.gene_id))
            .filter(EvolutionScore.convergence_count >= 5).scalar() or 0
        )

        # Top 20 by convergence_count
        top = (
            s.query(EvolutionScore.gene_id, EvolutionScore.convergence_count,
                    EvolutionScore.phylop_score)
            .filter(EvolutionScore.convergence_count != None)
            .order_by(EvolutionScore.convergence_count.desc(), EvolutionScore.phylop_score.desc())
            .limit(20)
            .all()
        )
        gene_ids = [r[0] for r in top]
        sym = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())
        data["top20_convergent_genes"] = [
            {"gene": sym.get(gid, gid), "convergence_count": cc,
             "phylop": round(float(pp), 3) if pp else None}
            for gid, cc, pp in top
        ]
    return data


def _collect_step7b() -> dict:
    from db.session import get_session
    from db.models import DivergentMotif, Ortholog, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step7b", "timestamp": _now()}
    with get_session() as s:
        data["motifs_with_convergent_aa"] = (
            s.query(func.count(DivergentMotif.id))
            .filter(DivergentMotif.convergent_aa_count >= 2).scalar() or 0
        )
        top = (
            s.query(Ortholog.gene_id, DivergentMotif.convergent_aa_count,
                    DivergentMotif.start_pos, DivergentMotif.human_seq, DivergentMotif.animal_seq)
            .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
            .filter(DivergentMotif.convergent_aa_count >= 2)
            .order_by(DivergentMotif.convergent_aa_count.desc())
            .limit(20)
            .all()
        )
        gene_ids = [r[0] for r in top]
        sym = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())
        data["top20_convergent_aa_motifs"] = [
            {"gene": sym.get(gid, gid), "convergent_aa_count": n, "pos": pos,
             "human_seq": hseq, "animal_seq": aseq}
            for gid, n, pos, hseq, aseq in top
        ]
    return data


def _collect_step8() -> dict:
    from db.session import get_session
    from db.models import ExpressionResult, Gene
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step8", "timestamp": _now()}
    with get_session() as s:
        data["expression_rows"] = s.query(func.count(ExpressionResult.id)).scalar() or 0
        data["species_with_expression"] = (
            s.query(func.count(func.distinct(ExpressionResult.comparison))).scalar() or 0
        )
        # DE genes: log2FC >= 1 and padj < 0.05
        de = (
            s.query(func.count(func.distinct(ExpressionResult.gene_id)))
            .filter(
                ExpressionResult.log2fc != None,
                func.abs(ExpressionResult.log2fc) >= 1,
                ExpressionResult.padj != None,
                ExpressionResult.padj < 0.05,
            )
            .scalar() or 0
        )
        data["de_genes"] = de

        # Top 20 DE by abs log2FC
        top = (
            s.query(ExpressionResult.gene_id, ExpressionResult.log2fc,
                    ExpressionResult.padj, ExpressionResult.comparison)
            .filter(
                ExpressionResult.log2fc != None,
                ExpressionResult.padj != None,
                ExpressionResult.padj < 0.05,
            )
            .order_by(func.abs(ExpressionResult.log2fc).desc())
            .limit(20)
            .all()
        )
        gene_ids = [r[0] for r in top]
        sym = dict(s.query(Gene.id, Gene.gene_symbol).filter(Gene.id.in_(gene_ids)).all())
        data["top20_de_genes"] = [
            {"gene": sym.get(gid, gid), "log2fc": round(float(l), 3) if l else None,
             "padj": round(float(p), 5) if p else None, "species": sp}
            for gid, l, p, sp in top
        ]
    return data


def _collect_step9() -> dict:
    from db.session import get_session
    from db.models import CandidateScore, Gene, EvolutionScore
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step9", "timestamp": _now()}
    with get_session() as s:
        tier_counts = dict(
            s.query(CandidateScore.tier, func.count(CandidateScore.gene_id))
            .group_by(CandidateScore.tier).all()
        )
        data["tier_counts"] = tier_counts

        # Full top-20 with all sub-scores
        top = (
            s.query(CandidateScore, Gene)
            .join(Gene, CandidateScore.gene_id == Gene.id)
            .order_by(CandidateScore.composite_score.desc())
            .limit(20)
            .all()
        )
        data["top20_candidates"] = [
            {
                "rank": i + 1,
                "gene": g.gene_symbol,
                "tier": cs.tier,
                "composite": round(float(cs.composite_score or 0), 3),
                "convergence": round(float(cs.convergence_score or 0), 3),
                "selection": round(float(cs.selection_score or 0), 3),
                "expression": round(float(cs.expression_score or 0), 3),
            }
            for i, (cs, g) in enumerate(top)
        ]

        # Score distribution histogram
        all_scores = [float(r[0]) for r in s.query(CandidateScore.composite_score).all() if r[0]]
        if all_scores:
            bins = {"0.0-0.2": 0, "0.2-0.4": 0, "0.4-0.6": 0, "0.6-0.8": 0, "0.8-1.0": 0}
            for sc in all_scores:
                if sc < 0.2: bins["0.0-0.2"] += 1
                elif sc < 0.4: bins["0.2-0.4"] += 1
                elif sc < 0.6: bins["0.4-0.6"] += 1
                elif sc < 0.8: bins["0.6-0.8"] += 1
                else: bins["0.8-1.0"] += 1
            data["score_distribution"] = bins
    return data


def _collect_step11() -> dict:
    from db.session import get_session
    from db.models import DiseaseAnnotation
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step11", "timestamp": _now()}
    with get_session() as s:
        data["disease_annotations"] = s.query(func.count(DiseaseAnnotation.id)).scalar() or 0
        data["gwas_hits"] = (
            s.query(func.count(DiseaseAnnotation.id))
            .filter(DiseaseAnnotation.gwas_hits != None,
                    DiseaseAnnotation.gwas_hits.cast("text") != "[]")
            .scalar() or 0
        )
        data["high_pli_genes"] = (
            s.query(func.count(DiseaseAnnotation.id))
            .filter(DiseaseAnnotation.gnomad_pli != None, DiseaseAnnotation.gnomad_pli >= 0.9)
            .scalar() or 0
        )
    return data


def _collect_step12() -> dict:
    from db.session import get_session
    from db.models import DrugTarget
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step12", "timestamp": _now()}
    with get_session() as s:
        data["drug_targets_annotated"] = s.query(func.count(DrugTarget.id)).scalar() or 0
        data["with_druggable_pocket"] = (
            s.query(func.count(DrugTarget.id))
            .filter(DrugTarget.pocket_count != None, DrugTarget.pocket_count > 0)
            .scalar() or 0
        )
        data["with_chembl_drug"] = (
            s.query(func.count(DrugTarget.id))
            .filter(DrugTarget.chembl_id != None)
            .scalar() or 0
        )
    return data


def _collect_step14() -> dict:
    from db.session import get_session
    from db.models import SafetyFlag
    from sqlalchemy import func

    data: dict[str, Any] = {"step": "step14", "timestamp": _now()}
    with get_session() as s:
        data["safety_flags_annotated"] = s.query(func.count(SafetyFlag.id)).scalar() or 0
        data["hub_genes"] = (
            s.query(func.count(SafetyFlag.id))
            .filter(SafetyFlag.string_degree != None, SafetyFlag.string_degree > 100)
            .scalar() or 0
        )
        data["depmap_essential"] = (
            s.query(func.count(SafetyFlag.id))
            .filter(SafetyFlag.depmap_score != None, SafetyFlag.depmap_score < -0.5)
            .scalar() or 0
        )
    return data


# ─────────────────────────────────────────────────────────────────────────────
# Validation gates
# ─────────────────────────────────────────────────────────────────────────────

def validate(step: str, data: dict) -> ValidationResult:
    vr = ValidationResult(step=step, status="PASS")

    if step == "step1":
        missing = data.get("tools_missing", [])
        critical_missing = [t for t in missing if t in ("orthofinder", "mafft", "diamond")]
        if critical_missing:
            vr.add("critical_tools", "FAIL", critical_missing, "all present",
                   f"Missing critical tools: {critical_missing}. Run scripts/setup_cloud.sh first.")
        elif missing:
            vr.add("optional_tools", "WARN", missing, "all present",
                   f"Optional tools missing: {missing}")
        else:
            vr.add("tools", "PASS", "all present", "all present")

        if not data.get("db_ok"):
            vr.add("database", "FAIL", "unreachable", "reachable",
                   f"DB error: {data.get('db_error')}. Check RDS_HOST env var.")
        else:
            vr.add("database", "PASS", f"{data.get('db_ping_ms')}ms", "reachable")

        if not data.get("gpu_available"):
            vr.add("gpu", "INFO", "not detected",
                   "optional",
                   "No GPU — step 4c (ESM embeddings) will run on CPU. Slower but functional.")
        else:
            vr.add("gpu", "PASS", data.get("gpu_name"), "optional")

    elif step == "step2":
        counts = data.get("protein_counts", {})
        total_downloaded = len(counts)
        vr.add("species_downloaded", "PASS" if total_downloaded >= 14 else "FAIL",
               total_downloaded, "≥14",
               "" if total_downloaded >= 14 else "Too few species downloaded — convergence signal will be inadequate.")

        for sid, cnt in counts.items():
            if cnt < 0:
                vr.add(f"proteome_{sid}", "WARN", "parse error", ">0 proteins",
                       f"Could not count proteins in {sid} FASTA — file may be corrupt.")
            elif sid in _VERTEBRATE_SPECIES and cnt < 10000:
                vr.add(f"proteome_{sid}", "FAIL" if cnt < 1000 else "WARN",
                       cnt, "≥10000 (vertebrate)",
                       f"{sid} has only {cnt} proteins — check NCBI assembly or retry download.")
            elif cnt < 3000:
                vr.add(f"proteome_{sid}", "WARN", cnt, "≥3000",
                       f"{sid} has {cnt} proteins (invertebrate; lower counts expected but may be incomplete).")

    elif step in ("step3", "step3b"):
        og_count = data.get("gene_count", 0)
        vr.add("gene_count", "PASS" if og_count >= 5000 else "FAIL",
               og_count, "≥5000",
               "" if og_count >= 5000 else "Very few genes — OrthoFinder may not have completed.")

        orth = data.get("ortholog_count", 0)
        vr.add("ortholog_count", "PASS" if orth >= 10000 else "WARN",
               orth, "≥10000")

        one2one = data.get("one_to_one_count", 0)
        vr.add("one_to_one", "PASS" if one2one >= 500 else "WARN",
               one2one, "≥500",
               "" if one2one >= 500 else "Fewer single-copy OGs than expected — tree quality may suffer.")

        # Check any species with zero coverage
        coverage = data.get("species_ortholog_coverage", {})
        zero_species = [s for s, n in coverage.items() if n == 0]
        if zero_species:
            vr.add("species_coverage", "WARN", zero_species, "all species present",
                   f"Species with zero ortholog membership: {zero_species}. Check reheader step.")

    elif step == "step3c":
        genes_scored = data.get("genes_scored", 0)
        genes_with_regions = data.get("genes_with_regions", 0)
        gene_coverage_pct = round(genes_with_regions / max(genes_scored, 1) * 100, 1) if genes_scored else 0

        vr.add("regions_extracted", "PASS" if genes_with_regions > 0 else "WARN",
               genes_with_regions, ">0 genes",
               "" if genes_with_regions > 0 else "No nucleotide regions extracted — check genome assembly accessions in species_registry.json.")

        # WARN if fewer than 50% of genes have promoter regions (invertebrates with poor annotations)
        promoter_count = data.get("regions_by_type", {}).get("promoter", 0)
        if genes_with_regions > 0:
            prom_pct = round(promoter_count / genes_with_regions * 100, 1)
            vr.add("promoter_coverage", "PASS" if prom_pct >= 50 else "WARN",
                   f"{prom_pct}%", "≥50%",
                   "" if prom_pct >= 50 else "Less than half of genes have promoter regions — GFF3 annotation may be incomplete for some species (expected for invertebrates).")

        div_events = data.get("total_regulatory_divergence_events", 0)
        conv_events = data.get("total_regulatory_convergence_events", 0)
        vr.add("regulatory_divergence_events", "INFO", div_events, "any",
               f"{div_events} divergence events, {conv_events} convergent events detected across all promoters.")

    elif step == "step3d":
        scored = data.get("genes_with_phylo_score", 0)
        if scored == 0:
            vr.add("phylo_tool", "INFO", "no scores computed", "optional",
                   "phyloP/PhastCons not installed or treefile not available — NULL scores are acceptable. "
                   "Install PHAST: conda install -c bioconda phast")
        else:
            vr.add("phylo_scores", "PASS", f"{scored} genes scored", f">{0}")
            acc = data.get("genes_accelerated_promoter", 0)
            vr.add("accelerated_promoter_genes", "INFO", acc, "any",
                   f"{acc} genes show accelerated promoter evolution (negative phyloP score) — candidates for regulatory adaptation.")

    elif step == "step4":
        motifs = data.get("motif_count", 0)
        genes = data.get("genes_with_motifs", 0)
        vr.add("gate2_output", "PASS" if genes >= 50 else "FAIL",
               genes, "≥50 genes",
               "" if genes >= 50 else "Gate 2 passed too few orthogroups. Lower convergence_min_lineages or divergence_identity_max.")

        lineage_cov = data.get("lineage_motif_coverage", {})
        zero_lineages = [lg for lg, n in lineage_cov.items() if n == 0]
        if zero_lineages:
            vr.add("lineage_coverage", "WARN", zero_lineages, "all lineages",
                   f"Lineages with zero motifs: {zero_lineages}. These species will not contribute convergence signal.")
        else:
            vr.add("lineage_coverage", "PASS", f"{len(lineage_cov)} lineages covered", "≥6 lineages")

        if motifs > 0:
            vr.add("motif_count", "PASS" if motifs >= 100 else "WARN", motifs, "≥100")

    elif step == "step4b":
        domain_pct = data.get("domain_pct", 0)
        vr.add("domain_coverage", "PASS" if domain_pct >= 20 else "WARN",
               f"{domain_pct}%", "≥20%",
               "" if domain_pct >= 20 else "Fewer motifs in domains than expected — may indicate generic surface noise.")

        am_pct = data.get("am_pct_coverage", 0)
        vr.add("alphamissense_coverage", "PASS" if am_pct >= 50 else "WARN",
               f"{am_pct}%", "≥50%",
               "" if am_pct >= 50 else "AlphaMissense covered fewer motifs than expected — possible API rate limit.")

    elif step == "step4c":
        pct = data.get("esm1v_pct_coverage", 0)
        vr.add("esm1v_coverage", "PASS" if pct >= 30 else "WARN",
               f"{pct}%", "≥30%",
               "" if pct >= 30 else "ESM-1v scored fewer motifs than expected — check fair-esm is installed.")

        mean_llr = data.get("mean_llr")
        if mean_llr is not None:
            vr.add("mean_llr", "PASS" if mean_llr < -0.5 else "WARN",
                   mean_llr, "< -0.5",
                   "" if mean_llr < -0.5 else "Mean LLR close to zero — motifs may be mostly neutral drift.")

    elif step == "step4d":
        dc = data.get("direction_counts", {})
        total = sum(dc.values())
        neutral = dc.get("neutral", 0) + dc.get("unclassified", 0)
        neutral_pct = round(neutral / total * 100, 1) if total else 100
        vr.add("direction_signal", "WARN" if neutral_pct > 80 else "PASS",
               f"{neutral_pct}% neutral/unclassified", "<80%",
               "" if neutral_pct <= 80 else "Most motifs unclassified — AM and ESM-1v scores may be incomplete.")

    elif step == "step5":
        if not data.get("treefile_exists"):
            vr.add("treefile", "FAIL", "missing", "exists",
                   "Treefile not produced. Check IQ-TREE2 log for errors.")
        else:
            vr.add("treefile", "PASS", "present", "exists")

        bs = data.get("bootstrap_mean")
        if bs is None:
            vr.add("bootstrap", "INFO", "ete3 not installed", "≥70 mean",
                   "Install ete3 for bootstrap stats: pip install ete3")
        elif bs < 50:
            vr.add("bootstrap", "FAIL", bs, "≥50",
                   "Very low bootstrap support — tree is unreliable. This will affect convergence attribution.")
        elif bs < 70:
            vr.add("bootstrap", "WARN", bs, "≥70",
                   "Moderate bootstrap support. Acceptable but review tree topology.")
        else:
            vr.add("bootstrap", "PASS", bs, "≥70")

    elif step == "step6":
        paml_scored = data.get("paml_scored_genes", 0)
        paml_sig = data.get("paml_significant_genes", 0)
        omega_pos = data.get("omega_positive_genes", 0)

        if paml_scored > 0:
            # PAML branch-site run — validate on LRT p-values
            vr.add("paml_scored_genes", "PASS" if paml_scored >= 1000 else "WARN",
                   paml_scored, "≥1000",
                   "" if paml_scored >= 1000 else "Fewer orthogroups scored than expected — some PAML jobs may have failed.")
            vr.add("paml_significant_genes", "PASS" if paml_sig >= 100 else "WARN",
                   paml_sig, "≥100",
                   f"{paml_sig} genes with LRT p<0.05 positive selection signal. {omega_pos} have foreground ω>1.")
        else:
            # Legacy HyPhy MEME path
            meme_pos = data.get("meme_positive_genes", 0)
            vr.add("meme_positive", "PASS" if meme_pos >= 100 else "WARN",
                   meme_pos, "≥100",
                   "" if meme_pos >= 100 else "Few MEME-positive genes — check HyPhy ran on all orthogroups.")

        benchmarks = data.get("benchmark_genes_in_top20", data.get("benchmark_genes_in_top20_meme", []))
        if benchmarks:
            vr.add("benchmark_recall", "PASS", benchmarks, "any known gene",
                   f"Known cancer genes in top selection results: {benchmarks}")
        else:
            vr.add("benchmark_recall", "INFO", "none in top 20",
                   "any known gene",
                   "Known benchmark genes (TP53, ATM, BRCA1) not in top 20 by p-value — "
                   "this is expected since PAML tests all OGs including those without divergent motifs.")

    elif step in ("step6b", "step6c"):
        pass  # informational only

    elif step in ("step7", "step7b"):
        if step == "step7":
            conv3 = data.get("genes_conv_ge3", 0)
            if conv3 == 0:
                vr.add("convergence_signal", "FAIL", 0, ">0 genes at ≥3 lineages",
                       "No genes with convergence across ≥3 lineages. Gate 2 may have been too permissive or species lineage groups are wrong.")
                vr.recommendation = "Check lineage_group values in species_registry.json and rerun step4."
            elif conv3 < 10:
                vr.add("convergence_signal", "WARN", conv3, "≥10 genes",
                       f"Only {conv3} genes at ≥3 lineages. Signal is present but weak.")
            else:
                vr.add("convergence_signal", "PASS", conv3, "≥10 genes")

            conv4 = data.get("genes_conv_ge4", 0)
            vr.add("strong_convergence", "PASS" if conv4 >= 3 else "INFO",
                   conv4, "≥3 genes at ≥4 lineages",
                   "" if conv4 >= 3 else "Few genes at ≥4 lineages — strong multi-lineage signal is sparse.")
        else:  # step7b
            motifs_with_conv = data.get("motifs_with_convergent_aa", 0)
            vr.add("convergent_aa", "PASS" if motifs_with_conv > 0 else "WARN",
                   motifs_with_conv, ">0 motifs",
                   "" if motifs_with_conv > 0 else "No motifs with convergent amino acids found.")

    elif step in ("step8", "step8b"):
        de = data.get("de_genes", 0)
        vr.add("de_genes", "INFO", de, "any >0",
               "Expression data is supplementary. Low counts expected for invertebrate species.")

    elif step == "step9":
        tier_counts = data.get("tier_counts", {})
        tier1 = tier_counts.get("Tier1", 0)
        tier2 = tier_counts.get("Tier2", 0)

        if tier1 == 0:
            vr.add("tier1_candidates", "FAIL", 0, ">0",
                   "No Tier1 candidates. Phase 1 signal is insufficient — investigate convergence and selection sub-scores.")
            vr.recommendation = "Check step7 (convergence_count distribution) and step6 (MEME signal). Do not proceed to Phase 2."
        elif tier1 < 5:
            vr.add("tier1_candidates", "WARN", tier1, "≥5",
                   f"Only {tier1} Tier1 candidates. Thresholds may be too strict for this species panel.")
        else:
            vr.add("tier1_candidates", "PASS", tier1, "≥5")

        vr.add("tier2_candidates", "PASS" if tier2 >= 10 else "WARN", tier2, "≥10")

        # Benchmark recall
        top20 = data.get("top20_candidates", [])
        top20_symbols = {c["gene"].upper() for c in top20 if c.get("gene")}
        found = _CANCER_BENCHMARKS & top20_symbols
        vr.add("benchmark_recall_top20", "PASS" if found else "WARN",
               list(found) or "none", "≥1 known gene",
               f"Known genes in top 20: {list(found)}" if found else
               "No known cancer resistance genes in top 20. Review sub-score breakdown.")
        if not found:
            vr.recommendation = (
                "Investigate: run `python scripts/benchmark_recall.py --trait cancer_resistance` "
                "for detailed breakdown. Common causes: low motif count, no MEME signal for known genes, "
                "or lineage group mis-assignment."
            )

    elif step == "step12":
        pockets = data.get("with_druggable_pocket", 0)
        total = data.get("drug_targets_annotated", 0)
        if total > 0:
            vr.add("druggable_pockets", "PASS" if pockets >= 5 else "WARN",
                   f"{pockets}/{total}", "≥5 genes",
                   "" if pockets >= 5 else "Few druggable pockets found — check fpocket installed and AlphaFold structures downloaded.")

    elif step == "step15":
        tier_counts = data.get("tier_counts", {})
        validated = tier_counts.get("Validated", 0)
        tier1 = tier_counts.get("Tier1", 0)
        vr.add("validated_candidates", "PASS" if validated >= 1 else "WARN",
               validated, "≥1",
               "" if validated >= 1 else "No 'Validated' tier genes (Tier1 + human genetics evidence). Check disease annotation coverage.")
        vr.add("final_tier1", "PASS" if tier1 >= 3 else "WARN", tier1, "≥3")

    return vr


# ─────────────────────────────────────────────────────────────────────────────
# Markdown report generator
# ─────────────────────────────────────────────────────────────────────────────

def _render_md(step: str, data: dict, vr: ValidationResult) -> str:
    from pipeline.config import get_config, get_deployment
    cfg = get_config()
    thresholds = cfg.get("thresholds", {})
    phenotype = cfg.get("target_phenotype", "")

    lines = [
        f"# Step Report: {step}",
        f"",
        f"**Generated:** {_now()}  ",
        f"**Deployment:** {get_deployment()}  ",
        f"**Phenotype:** {phenotype}  ",
        f"",
        f"---",
        f"",
        f"## What this step did",
        f"",
        _EXPLANATIONS.get(step, "No explanation available."),
        f"",
        f"---",
        f"",
        f"## Active thresholds",
        f"",
    ]
    for k, v in thresholds.items():
        lines.append(f"- `{k}`: {v}")
    lines += ["", "---", "", "## Validation", ""]

    status_icon = {"PASS": "✅", "WARN": "⚠️", "FAIL": "❌"}.get(vr.status, "ℹ️")
    lines.append(f"**Overall status: {status_icon} {vr.status}**")
    lines.append("")

    for chk in vr.checks:
        icon = {"PASS": "✅", "WARN": "⚠️", "FAIL": "❌", "INFO": "ℹ️"}.get(chk.result, "ℹ️")
        lines.append(f"- {icon} **{chk.name}**: `{chk.value}` (threshold: {chk.threshold})")
        if chk.message:
            lines.append(f"  - {chk.message}")

    if vr.recommendation:
        lines += ["", f"> **Recommendation:** {vr.recommendation}"]

    lines += ["", "---", "", "## Results", ""]

    # Render key data per step
    if step == "step1":
        lines.append(f"- Tools found: {', '.join(data.get('tools_found', []))}")
        if data.get("tools_missing"):
            lines.append(f"- Tools missing: {', '.join(data['tools_missing'])}")
        lines.append(f"- DB ping: {data.get('db_ping_ms')}ms  OK: {data.get('db_ok')}")
        lines.append(f"- GPU: {data.get('gpu_name') or 'not detected'}")

    elif step == "step2":
        counts = data.get("protein_counts", {})
        lines.append("| Species | Protein count |")
        lines.append("|---------|--------------|")
        for sid, n in sorted(counts.items()):
            lines.append(f"| {sid} | {n:,} |")

    elif step in ("step3", "step3b"):
        lines.append(f"- Genes (human proteins): {data.get('gene_count', 0):,}")
        lines.append(f"- Ortholog rows: {data.get('ortholog_count', 0):,}")
        lines.append(f"- 1:1 orthogroups: {data.get('one_to_one_count', 0):,}")
        lines.append("")
        cov = data.get("species_ortholog_coverage", {})
        if cov:
            lines.append("**Species ortholog coverage:**")
            lines.append("| Species | Orthologs |")
            lines.append("|---------|----------|")
            for sid, n in sorted(cov.items(), key=lambda x: -x[1]):
                lines.append(f"| {sid} | {n:,} |")

    elif step == "step3c":
        lines.append(f"- Nucleotide regions extracted: {data.get('nucleotide_regions_total', 0):,}")
        lines.append(f"- Genes with regions: {data.get('genes_with_regions', 0):,}")
        lines.append(f"- Species with regions: {data.get('species_with_regions', 0)}")
        lines.append(f"- Genes with alignment scores: {data.get('genes_scored', 0):,}")
        lines.append(f"- Total regulatory divergence events: {data.get('total_regulatory_divergence_events', 0):,}")
        lines.append(f"- Total regulatory convergence events: {data.get('total_regulatory_convergence_events', 0):,}")
        rbt = data.get("regions_by_type", {})
        if rbt:
            lines += ["", "**Regions by type:**"]
            for rt, n in sorted(rbt.items()):
                lines.append(f"- {rt}: {n:,}")
        lines += ["", "**Mean conservation by region:**"]
        for rtype in ("cds", "promoter", "downstream"):
            val = data.get(f"mean_{rtype}_conservation")
            lines.append(f"- {rtype}: {val if val is not None else '—'}")
        top = data.get("top10_promoter_divergence", [])
        if top:
            lines += ["", "**Top 10 genes by promoter divergence count:**", "",
                      "| Gene | Divergence count |", "|------|-----------------|"]
            for r in top:
                lines.append(f"| {r['gene']} | {r['divergence_count']} |")

    elif step == "step3d":
        lines.append(f"- Genes with phylo scores: {data.get('genes_with_phylo_score', 0):,}")
        for field in ("cds_phylo_score", "promoter_phylo_score", "downstream_phylo_score"):
            val = data.get(f"mean_{field}")
            lines.append(f"- Mean {field}: {val if val is not None else '— (tool not installed)'}")
        lines.append(f"- Genes with accelerated promoter evolution: {data.get('genes_accelerated_promoter', 0):,}")
        top = data.get("top10_accelerated_promoter", [])
        if top:
            lines += ["", "**Top 10 genes by accelerated promoter evolution:**", "",
                      "| Gene | Promoter phyloP score |", "|------|----------------------|"]
            for r in top:
                lines.append(f"| {r['gene']} | {r['promoter_phylo_score']} |")

    elif step in ("step4", "step4b", "step4c", "step4d"):
        lines.append(f"- Divergent motifs: {data.get('motif_count', 0):,}")
        lines.append(f"- Genes with motifs: {data.get('genes_with_motifs', 0):,}")
        lc = data.get("lineage_motif_coverage", {})
        if lc:
            lines.append("")
            lines.append("**Motifs per lineage:**")
            for lg, n in sorted(lc.items(), key=lambda x: -x[1]):
                lines.append(f"- {lg}: {n:,}")
        top = data.get("top20_genes_by_motifs", [])
        if top:
            lines += ["", "**Top 20 genes by species breadth (most independently diverged):**", ""]
            lines.append("| Gene | Species diverged | Motif count | Avg divergence |")
            lines.append("|------|-----------------|-------------|----------------|")
            for r in top:
                avg = r.get("avg_divergence")
                avg_str = f"{avg:.4f}" if avg is not None else "—"
                sp = r.get("species_count", "—")
                lines.append(f"| {r['symbol']} | {sp} | {r['motif_count']} | {avg_str} |")
        if step == "step4b":
            lines.append(f"\n- Motifs in functional domains: {data.get('motifs_in_functional_domain', 0):,} ({data.get('domain_pct', 0)}%)")
            lines.append(f"- AlphaMissense scored: {data.get('am_scored_count', 0):,} ({data.get('am_pct_coverage', 0)}%)")
            lines.append(f"- Likely pathogenic (AM > 0.564): {data.get('am_likely_pathogenic', 0):,}")
        if step == "step4c":
            lines.append(f"\n- ESM-1v scored: {data.get('esm1v_scored', 0):,} ({data.get('esm1v_pct_coverage', 0)}%)")
            lines.append(f"- Mean LLR: {data.get('mean_llr')}")
            dist = data.get("llr_distribution", {})
            if dist:
                lines.append("\n**LLR distribution:**")
                for b, n in dist.items():
                    lines.append(f"- {b}: {n:,}")
        if step == "step4d":
            dc = data.get("direction_counts", {})
            lines.append("\n**Direction counts:**")
            for d, n in dc.items():
                lines.append(f"- {d}: {n:,}")

    elif step == "step5":
        lines.append(f"- Treefile: {data.get('treefile_path')}")
        lines.append(f"- Bootstrap mean: {data.get('bootstrap_mean')}")
        lines.append(f"- Bootstrap pct above 90: {data.get('bootstrap_pct_above_90')}%")
        if data.get("newick"):
            lines += ["", "**Newick topology:**", "", "```", data["newick"][:300] + ("..." if len(data["newick"]) > 300 else ""), "```"]

    elif step in ("step6", "step6b", "step6c"):
        if step == "step6":
            lines.append(f"- Genes with selection scores: {data.get('genes_with_selection', 0):,}")
            if data.get("paml_scored_genes", 0) > 0:
                lines.append(f"- PAML branch-site scored: {data.get('paml_scored_genes', 0):,}")
                lines.append(f"- PAML significant (LRT p<0.05): **{data.get('paml_significant_genes', 0):,}**")
                lines.append(f"- Foreground ω>1 (positive selection): {data.get('omega_positive_genes', 0):,}")
                lines.append(f"- Mean ω for significant genes: {data.get('mean_omega_significant')}")
                mc = data.get("model_counts", {})
                if mc:
                    lines.append(f"- Models: {', '.join(f'{k}: {v}' for k,v in mc.items())}")
            else:
                lines.append(f"- MEME-positive genes (p<0.1): {data.get('meme_positive_genes', 0):,}")
                lines.append(f"- Mean MEME q-value: {data.get('avg_meme_pvalue')}")
            bm = data.get("benchmark_genes_in_top20", data.get("benchmark_genes_in_top20_meme", []))
            if bm:
                lines.append(f"- Known benchmark genes in top 20: **{', '.join(bm)}**")
            top = data.get("top20_paml_genes", data.get("top20_meme_genes", []))
            if top:
                if data.get("paml_scored_genes", 0) > 0:
                    lines += ["", "**Top 20 by PAML branch-site LRT p-value:**", "", "| Gene | LRT p-value | dN/dS (ω) |", "|------|------------|-----------|"]
                    for r in top:
                        lines.append(f"| {r['gene']} | {r.get('dnds_pvalue', r.get('meme_qvalue'))} | {r.get('dnds_ratio', r.get('dnds'))} |")
                else:
                    lines += ["", "**Top 20 by MEME evidence:**", "", "| Gene | MEME q-value | dN/dS |", "|------|-------------|-------|"]
                    for r in top:
                        lines.append(f"| {r['gene']} | {r.get('meme_qvalue')} | {r.get('dnds')} |")
        elif step == "step6b":
            lines.append(f"- FEL-positive genes: {data.get('fel_positive_genes', 0):,}")
            lines.append(f"- BUSTED-positive genes: {data.get('busted_positive_genes', 0):,}")
            lines.append(f"- MEME+BUSTED both positive: {data.get('meme_and_busted_both_positive', 0):,}")
        elif step == "step6c":
            lines.append(f"- RELAX intensified (k>1, p<0.05): {data.get('relax_intensified', 0):,}")
            lines.append(f"- RELAX relaxed (k<1, p<0.05): {data.get('relax_relaxed', 0):,}")

    elif step in ("step7", "step7b"):
        lines.append(f"- Genes at convergence ≥3: **{data.get('genes_conv_ge3', 0):,}**")
        lines.append(f"- Genes at convergence ≥4: {data.get('genes_conv_ge4', 0):,}")
        lines.append(f"- Genes at convergence ≥5: {data.get('genes_conv_ge5', 0):,}")
        dist = data.get("convergence_count_distribution", {})
        if dist:
            lines += ["", "**Convergence count distribution:**", ""]
            for k, v in sorted(dist.items()):
                lines.append(f"- {k} lineages: {v:,} genes")
        top = data.get("top20_convergent_genes", [])
        if top:
            lines += ["", "**Top 20 convergent genes:**", "", "| Gene | Convergence count | PhyloP |", "|------|------------------|--------|"]
            for r in top:
                lines.append(f"| {r['gene']} | {r['convergence_count']} | {r['phylop']} |")

    elif step in ("step8", "step8b"):
        lines.append(f"- Expression rows: {data.get('expression_rows', 0):,}")
        lines.append(f"- Species with expression data: {data.get('species_with_expression', 0)}")
        lines.append(f"- DE genes (|log2FC| ≥1, padj < 0.05): {data.get('de_genes', 0):,}")
        top = data.get("top20_de_genes", [])
        if top:
            lines += ["", "**Top 20 DE genes:**", "", "| Gene | log2FC | padj | Species |", "|------|--------|------|---------|"]
            for r in top:
                lines.append(f"| {r['gene']} | {r['log2fc']} | {r['padj']} | {r['species']} |")

    elif step in ("step9", "step15"):
        tc = data.get("tier_counts", {})
        lines.append(f"- Tier1: **{tc.get('Tier1', 0)}**")
        lines.append(f"- Tier2: {tc.get('Tier2', 0)}")
        lines.append(f"- Tier3: {tc.get('Tier3', 0)}")
        lines.append(f"- Validated: {tc.get('Validated', 0)}")
        top = data.get("top20_candidates", [])
        if top:
            lines += ["", "**Top 20 candidates:**", "",
                      "| Rank | Gene | Tier | Composite | Convergence | Selection | Expression |",
                      "|------|------|------|-----------|-------------|-----------|------------|"]
            for r in top:
                lines.append(
                    f"| {r['rank']} | **{r['gene']}** | {r['tier']} | {r['composite']} "
                    f"| {r['convergence']} | {r['selection']} | {r['expression']} |"
                )

    elif step in ("step11", "step11b", "step11c", "step11d"):
        lines.append(f"- Disease annotations: {data.get('disease_annotations', 0):,}")
        lines.append(f"- Genes with GWAS hits: {data.get('gwas_hits', 0):,}")
        lines.append(f"- High pLI genes (≥0.9): {data.get('high_pli_genes', 0):,}")

    elif step in ("step12", "step12b"):
        lines.append(f"- Drug targets annotated: {data.get('drug_targets_annotated', 0):,}")
        lines.append(f"- With druggable pocket: {data.get('with_druggable_pocket', 0):,}")
        lines.append(f"- With ChEMBL drug: {data.get('with_chembl_drug', 0):,}")

    elif step in ("step14", "step14b"):
        lines.append(f"- Safety flags annotated: {data.get('safety_flags_annotated', 0):,}")
        lines.append(f"- Hub genes (>100 STRING edges): {data.get('hub_genes', 0):,}")
        lines.append(f"- DepMap essential (Chronos < -0.5): {data.get('depmap_essential', 0):,}")

    lines += ["", "---", f"*Report generated by pipeline/step_reporter.py*"]
    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Step data router
# ─────────────────────────────────────────────────────────────────────────────

def collect(step: str) -> dict:
    sys.path.insert(0, str(_REPO_ROOT))
    collectors = {
        "step1": _collect_step1,
        "step2": _collect_step2,
        "step3": _collect_step3,
        "step3b": _collect_step3,
        "step3c": _collect_step3c,
        "step3d": _collect_step3d,
        "step4": _collect_step4,
        "step4b": _collect_step4b,
        "step4c": _collect_step4c,
        "step4d": _collect_step4d,
        "step5": _collect_step5,
        "step6": _collect_step6,
        "step6b": _collect_step6b,
        "step6c": _collect_step6c,
        "step7": _collect_step7,
        "step7b": _collect_step7b,
        "step8": _collect_step8,
        "step8b": _collect_step8,
        "step9": _collect_step9,
        "step10": lambda: {"step": "step10", "timestamp": _now(), "note": "API ready — no DB data collected"},
        "step10b": lambda: {"step": "step10b", "timestamp": _now(), "note": "AlphaGenome track — see regulatory_divergence table"},
        "step11": _collect_step11,
        "step11b": _collect_step11,
        "step11c": _collect_step11,
        "step11d": _collect_step11,
        "step12": _collect_step12,
        "step12b": _collect_step12,
        "step13": lambda: {"step": "step13", "timestamp": _now(), "note": "See gene_therapy_score table"},
        "step14": _collect_step14,
        "step14b": _collect_step14,
        "step15": _collect_step9,  # same scoring table
        "step16": lambda: {"step": "step16", "timestamp": _now(), "note": "Pipeline complete"},
    }
    fn = collectors.get(step)
    if fn is None:
        return {"step": step, "timestamp": _now(), "note": "No collector defined for this step"}
    try:
        return fn()
    except Exception as e:
        log.warning("Collector failed for %s: %s", step, e)
        return {"step": step, "timestamp": _now(), "error": str(e)}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


# ─────────────────────────────────────────────────────────────────────────────
# Report writing
# ─────────────────────────────────────────────────────────────────────────────

def write_report(step: str) -> tuple[Path, Path, ValidationResult]:
    """Collect data, validate, write JSON + MD. Returns (json_path, md_path, vr)."""
    _CACHE_DIR.mkdir(exist_ok=True)
    data = collect(step)
    vr = validate(step, data)
    data["_validation"] = vr.to_dict()

    json_path = _CACHE_DIR / f"{step}.json"
    json_path.write_text(json.dumps(data, indent=2, default=str))

    md_path = _CACHE_DIR / f"{step}.md"
    md_path.write_text(_render_md(step, data, vr))

    # Update run metadata
    meta_path = _CACHE_DIR / "run_metadata.json"
    try:
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    except Exception:
        meta = {}
    meta[step] = {"completed_at": _now(), "validation_status": vr.status}
    if "started_at" not in meta:
        meta["started_at"] = _now()
    meta_path.write_text(json.dumps(meta, indent=2))

    return json_path, md_path, vr


def show_report(step: str) -> None:
    """Pretty-print the Markdown report to stdout."""
    md_path = _CACHE_DIR / f"{step}.md"
    if md_path.exists():
        print(md_path.read_text())
    else:
        print(f"No report found at {md_path}. Run without --show first.")


def print_validation(vr: ValidationResult) -> None:
    """Print validation summary to terminal (used by stepwise runner)."""
    status_icon = {"PASS": "✅", "WARN": "⚠ ", "FAIL": "❌"}.get(vr.status, "ℹ ")
    print(f"\n  Validation: {status_icon} {vr.status}")
    for chk in vr.checks:
        icon = {"PASS": "✅", "WARN": "⚠ ", "FAIL": "❌", "INFO": "ℹ "}.get(chk.result, "ℹ ")
        print(f"  {icon} {chk.name}: {chk.value}  (threshold: {chk.threshold})")
        if chk.message:
            print(f"     {chk.message}")
    if vr.recommendation:
        print(f"\n  Recommendation: {vr.recommendation}")
    print()


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry point
# ─────────────────────────────────────────────────────────────────────────────

def main() -> int:
    import os
    os.chdir(_REPO_ROOT)
    sys.path.insert(0, str(_REPO_ROOT))

    parser = argparse.ArgumentParser(description="BioResilient per-step reporter")
    parser.add_argument("--step", required=True, help="Step name e.g. step4")
    parser.add_argument("--show", action="store_true", help="Print Markdown report to terminal")
    parser.add_argument("--validate", action="store_true", help="Print validation result only")
    args = parser.parse_args()

    if args.show:
        show_report(args.step)
        return 0

    json_path, md_path, vr = write_report(args.step)
    print(f"  Report written: {json_path.name}  {md_path.name}")
    print_validation(vr)

    if args.validate:
        return 0 if vr.status != "FAIL" else 1

    return 0 if vr.status != "FAIL" else 1


if __name__ == "__main__":
    sys.exit(main())
