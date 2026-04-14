"""Step 9 — Composite score assembly.

Phase 1 uses a rank-product method instead of arbitrary weighted composites:

Evidence layers (Phase 1):
  1. Selection:           PAML branch-site LRT p-value (lower = stronger positive selection)
  2. Convergence:         PhyloP score (higher = more convergent evolution)
  3. Convergent AAs:      max convergent_aa_count per gene (higher = more convergent)
  4. Functional evidence: Open Targets + GTEx + DepMap combined score (step 8, optional)

Scoring procedure:
  1. Rank each gene within each evidence layer (best = rank 1)
  2. Compute rank product: RP = product(rank_i) / n^k (normalised)
  3. Compute rank-product p-value via permutation (1000 shuffles)
  4. Apply Benjamini-Hochberg FDR correction
  5. Assign tiers: Tier1 = FDR < 0.05, Tier2 = FDR < 0.20, Tier3 = rest

Reference: Breitling et al. (2004) "Rank products" FEBS Letters 573:83-92
"""

import logging
import math
from collections import defaultdict
from datetime import datetime, timezone
from typing import Optional

from sqlalchemy import func

from db.models import (
    CandidateScore,
    ConvergentPositionAnnotation,
    DiseaseAnnotation,
    DivergentMotif,
    DrugTarget,
    EvolutionScore,
    Gene,
    NucleotideScore,
    Ortholog,
    PhyloConservationScore,
    RegulatoryDivergence,
    SafetyFlag,
    Species,
)
from db.session import get_session
from pipeline.config import get_scoring_weights, get_tier_thresholds, get_thresholds

_test_species_cache: Optional[frozenset[str]] = None


def _get_test_species() -> frozenset[str]:
    """Return the set of test (target-phenotype) species IDs from the DB.

    Excludes baseline, control, and outgroup species — the same classification
    used by RELAX branch labeling.  Cached for the process lifetime.
    """
    global _test_species_cache
    if _test_species_cache is not None:
        return _test_species_cache
    try:
        with get_session() as session:
            rows = session.query(Species.id, Species.phenotypes, Species.is_control).all()
        _test_species_cache = frozenset(
            r.id for r in rows
            if (r.phenotypes or [])
            and "outgroup" not in (r.phenotypes or [])
            and not r.is_control
            and "baseline" not in (r.phenotypes or [])
        )
    except Exception:
        _test_species_cache = frozenset()
    return _test_species_cache
from pipeline.stats import apply_bh_correction, apply_bh_to_evolution_scores

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sub-score computation
# ---------------------------------------------------------------------------


def convergence_score(convergence_count: int, max_lineages: int = 8, phylo_weight: Optional[float] = None) -> float:
    """Score in [0, 1] based on convergent evolution signal.

    When phylo_weight is available (from phylogenetically-weighted convergence),
    it is used directly after normalisation — this rewards convergence across
    more phylogenetically distant lineages (e.g. Rodents + Birds scores higher
    than Rodents + other Rodents).

    When only convergence_count is available (test runs, legacy data), falls
    back to a sigmoid on count alone.

    Max expected phylo_weight:
      - 8 lineages at 310 MY avg spacing ≈ 8 × 2.06 = 16.5
      - Normalise by capping at 16.0 for a [0,1] range.
    """
    if phylo_weight is not None and phylo_weight > 0:
        return round(min(phylo_weight / 16.0, 1.0), 4)
    if convergence_count is None or convergence_count <= 0:
        return 0.0
    return round(1.0 / (1.0 + math.exp(-1.5 * (convergence_count - 2))), 4)


def selection_score(
    dnds_ratio: Optional[float],
    dnds_pvalue: Optional[float],
    fel_sites: Optional[int] = None,
    busted_pvalue: Optional[float] = None,
    relax_k: Optional[float] = None,
    relax_pvalue: Optional[float] = None,
    branches_under_selection: Optional[list[str]] = None,
    selection_model: Optional[str] = None,
) -> float:
    """Score in [0, 1] combining selection evidence.

    PAML branch-site model A (selection_model="paml_branch_site"):
        70% p-value component  — LRT p-value (H0 vs H1, chi-squared df=2)
        30% omega component    — foreground dN/dS in positive selection class
        This is the primary pathway going forward.

    HyPhy suite (selection_model="busted_ph", legacy):
        FEL       35%  — pervasive positive selection (site-level)
        BUSTED    30%  — gene-wide episodic selection
        BUSTED-PH 25%  — phenotype-specific episodic selection (foreground)
        RELAX     10%  — selection intensity shift in test vs reference branches

    Legacy MEME results (selection_model="meme") use backward-compatible scoring.
    """
    # --- PAML branch-site model A (primary pathway) ---
    if selection_model == "paml_branch_site":
        if dnds_pvalue is None or dnds_pvalue >= 1.0:
            return 0.0
        # Guard: ω ≤ 1 means no positive selection regardless of p-value.
        # Stored p=0 with ω=1 is a PAML numerical artifact (LRT=0 → p should be 1.0).
        if dnds_ratio is None or dnds_ratio <= 1.0:
            return 0.0
        p_score = round(min(-math.log10(max(dnds_pvalue, 1e-10)) / 10.0, 1.0), 4)
        omega_score = round(min((dnds_ratio or 0.0) / 5.0, 1.0), 4)
        return round(0.70 * p_score + 0.30 * omega_score, 4)

    # --- BUSTED-PH / MEME component (phenotype-specific, stored in dnds_pvalue) ---
    # New HyPhy runs: selection_model="busted_ph", dnds_pvalue = BUSTED-PH p-value
    # Legacy:         selection_model="meme",      dnds_pvalue = MEME pseudo p-value
    pheno_score = 0.0
    if dnds_pvalue is not None and 0 < dnds_pvalue <= 1:
        pheno_score = round(min(-math.log10(dnds_pvalue) / 10.0, 1.0), 4)
    elif selection_model == "meme" and dnds_ratio is not None and dnds_pvalue is not None:
        dnds_norm = min(dnds_ratio / 5.0, 1.0)
        if dnds_pvalue <= 0:
            pval_weight = 1.0
        elif dnds_pvalue >= 1:
            pval_weight = 0.0
        else:
            pval_weight = min(-math.log10(dnds_pvalue) / 10.0, 1.0)
        pheno_score = round(0.5 * dnds_norm + 0.5 * pval_weight, 4)
        if pheno_score > 0 and branches_under_selection:
            test_set = _get_test_species()
            if test_set:
                resistant_count = sum(1 for b in branches_under_selection if b in test_set)
                resistant_frac = resistant_count / len(branches_under_selection)
                pheno_score = round(pheno_score * (0.5 + 0.5 * resistant_frac), 4)

    # If no supplementary tests, return phenotype score scaled down
    if fel_sites is None and busted_pvalue is None and relax_k is None:
        return round(pheno_score * 0.65, 4)

    # --- FEL component (35%) ---
    fel_score = round(min(fel_sites / 10.0, 1.0), 4) if fel_sites and fel_sites > 0 else 0.0

    # --- BUSTED gene-wide component (30%) ---
    if busted_pvalue is not None and 0 < busted_pvalue <= 1:
        busted_score = round(min(-math.log10(busted_pvalue) / 10.0, 1.0), 4)
    else:
        busted_score = 0.0

    # --- RELAX component (10%): intensification (k>1) with p-value significance ---
    if relax_k is not None and relax_k > 1.0 and relax_pvalue is not None and 0 < relax_pvalue <= 1:
        k_component = min((relax_k - 1.0) / 2.0, 1.0)
        p_component = min(-math.log10(relax_pvalue) / 10.0, 1.0)
        relax_score = round(k_component * p_component, 4)
    else:
        relax_score = 0.0

    return round(0.35 * fel_score + 0.30 * busted_score + 0.25 * pheno_score + 0.10 * relax_score, 4)


def expression_score_from_db(gene_id: str, session) -> float:
    """Retrieve expression_score already stored in CandidateScore (from Step 8).

    Note: This function is for external callers and tests only. The main
    run_scoring() path uses a pre-fetched bulk map (expr_map) to avoid N+1 queries.
    """
    cs = session.query(CandidateScore).filter_by(gene_id=gene_id, trait_id="").first()
    if cs is None:
        return 0.0
    return cs.expression_score or 0.0


def disease_score(
    ann: Optional[DiseaseAnnotation],
    phenotype_disease_ids: Optional[set[str]] = None,
) -> float:
    """Score in [0, 1] from OpenTargets, GWAS, gnomAD, known drugs, protective variants.

    Scoring breakdown (max 1.0):
      0.30  OpenTargets association score — disease relevance breadth
      0.20  GWAS p-value                 — human genetic evidence
      0.15  gnomAD pLI or LOEUF          — LoF intolerance (functional importance)
      0.20  known drug phase             — clinical validation (phenotype-aware):
              Full credit (0.20) if the drug's indication overlaps with the target phenotype
              Partial credit (0.08) if the drug demonstrates tractability but for a different disease
              No credit (0.00) if phase is unknown
      0.15  protective variant           — PCSK9-paradigm human validation

    phenotype_disease_ids: EFO/MONDO IDs for the current phenotype from
        functional_evidence_config.json. When provided, drug phase credit is
        conditional on indication overlap to avoid rewarding pharmaceutical attention
        for off-phenotype drugs (e.g. a Phase 4 statin should not boost a
        cancer_resistance run).
    """
    if ann is None:
        return 0.0
    s = 0.0

    # OpenTargets association breadth
    if ann.opentargets_score is not None:
        s += min(ann.opentargets_score, 1.0) * 0.30

    # Human genetic evidence (GWAS)
    if ann.gwas_pvalue is not None and ann.gwas_pvalue > 0:
        s += min(-math.log10(ann.gwas_pvalue) / 15.0, 1.0) * 0.20

    # LoF constraint — prefer LOEUF (v4, continuous) over pLI (binary-ish).
    # LOEUF < 0.6 = intolerant; invert so low LOEUF → high score.
    loeuf = getattr(ann, "gnomad_loeuf", None)
    pli = ann.gnomad_pli
    if loeuf is not None:
        # 1 − LOEUF/1.0 maps [0, 1] LOEUF to [1, 0] score; clamp at 0
        s += max(0.0, min(1.0 - loeuf, 1.0)) * 0.15
    elif pli is not None:
        s += min(pli, 1.0) * 0.15

    # Known drug phase — phenotype-aware credit
    # Full credit: drug is approved/studied for an indication that matches the
    #   current phenotype's disease ontology IDs → confirms clinical relevance.
    # Partial credit: drug exists (target is tractable) but indication is off-phenotype
    #   → confirms the gene can be modulated, but not for this phenotype.
    known_phase = getattr(ann, "known_drug_phase", None)
    if known_phase is not None and known_phase > 0:
        indication_ids: list = getattr(ann, "known_drug_indication_ids", None) or []
        phase_fraction = min(known_phase / 4.0, 1.0)
        if phenotype_disease_ids and indication_ids:
            # Normalise OT ID format: "EFO_0000311" ↔ "EFO:0000311"
            normalised_indications = {
                i.replace("_", ":").upper() for i in indication_ids
            }
            normalised_phenotype = {
                i.replace("_", ":").upper() for i in phenotype_disease_ids
            }
            if normalised_indications & normalised_phenotype:
                # Drug is relevant to this phenotype → full clinical validation credit
                s += phase_fraction * 0.20
            else:
                # Drug is tractability evidence only — gene can be modulated, but
                # indication is off-phenotype. Reduced credit (0.08 max).
                s += phase_fraction * 0.08
        else:
            # No indication data stored yet (old pipeline runs) or no phenotype context
            # provided — fall back to full credit to avoid penalising existing results.
            s += phase_fraction * 0.20

    # Rare protective variant (PCSK9 paradigm) — binary bonus
    pv_count = ann.protective_variant_count
    pv_pval = ann.protective_variant_pvalue
    if pv_count is not None and pv_count >= 1 and pv_pval is not None and pv_pval < 5e-8:
        s += 0.15

    return round(min(s, 1.0), 4)


def druggability_score(dt: Optional[DrugTarget]) -> float:
    """Score in [0, 1] from pockets, ChEMBL, CanSAR tier, P2Rank, tractability, and convergent proximity.

    Scoring breakdown (max 1.0):
      0.20  fpocket pocket count   — proxy for surface pocketability
      0.20  fpocket top pocket score — direct druggability estimate
      0.10  P2Rank ML score        — independent ML pocket confidence
      0.15  ChEMBL target / existing drugs — chemical matter exists
      0.10  CanSAR tier            — curated druggability annotation
      0.10  OT tractability (SM/AB/PROTAC) — any modality tractable
      0.15  Convergent-pocket proximal — convergent residue near top pocket
    """
    if dt is None:
        return 0.0
    s = 0.0

    # Pocket existence and quality
    if dt.pocket_count is not None and dt.pocket_count > 0:
        s += min(dt.pocket_count / 5.0, 1.0) * 0.20
    if dt.top_pocket_score is not None:
        s += min(dt.top_pocket_score, 1.0) * 0.20

    # P2Rank ML prediction
    p2rank = getattr(dt, "p2rank_score", None)
    if p2rank is not None and p2rank > 0:
        s += min(p2rank, 1.0) * 0.10

    # Chemical matter evidence
    if dt.chembl_target_id:
        s += 0.10
    if dt.existing_drugs and len(dt.existing_drugs) > 0:
        s += 0.05

    # CanSAR tier
    tier = (dt.druggability_tier or "").upper()
    if tier == "A":
        s += 0.10
    elif tier == "B":
        s += 0.05

    # OpenTargets tractability — any modality tractable is a positive signal
    sm = getattr(dt, "tractability_sm", None)
    ab = getattr(dt, "tractability_ab", None)
    protac = getattr(dt, "tractability_protac", None)
    if sm or ab or protac:
        s += 0.10

    # Convergent-pocket proximal — the key novel signal:
    # if the convergent residue sits near the druggable pocket, modulating that
    # pocket directly perturbs the evolved adaptive position.
    proximal = getattr(dt, "convergent_pocket_proximal", None)
    if proximal:
        s += 0.15

    return round(min(s, 1.0), 4)


def safety_score(sf: Optional[SafetyFlag]) -> float:
    """Score in [0, 1]; high hub_risk, essential, large family reduce score (inverted for composite).

    Extended with DepMap CRISPR essentiality, GTEx expression breadth, and PheWAS disease
    association breadth (step 14):

      DepMap:
        - depmap_score > 0.7: very broadly essential → serious toxicity concern (-0.30)
        - depmap_score > 0.5: broadly essential (-0.15)
        - depmap_score > 0.3: moderate essentiality (-0.05)

      GTEx breadth:
        - tissue_count > 40: ubiquitous → off-target risk (-0.20)
        - tissue_count > 25: moderate breadth (-0.10)
        - tissue_count > 15: slight breadth (-0.05)

      GTEx depth:
        - max_tpm > 1000: critical physiological role in at least one tissue (-0.15)
        - max_tpm > 500: high peak expression (-0.08)

      PheWAS (step 14):
        Genes associated with many human diseases carry liability risk — modifying
        them pharmacologically is more likely to cause unintended phenotypes.
        The score is derived from phewas_hits, a dict {disease: association_score}.
        We penalise based on the count of high-confidence associations (score > 0.3)
        and the maximum association score.
          - ≥ 15 disease associations at score > 0.3 → high pleiotropic risk (-0.20)
          - ≥  5 disease associations at score > 0.3 → moderate risk (-0.10)
          - max association score > 0.7 → severe on-target liability (-0.10)
    """
    if sf is None:
        return 0.5
    s = 1.0
    if sf.hub_risk:
        s -= 0.1
    if sf.is_essential:
        s -= 0.2
    if sf.family_size is not None and sf.family_size > 100:
        s -= 0.2
    depmap = getattr(sf, "depmap_score", None)
    if depmap is not None:
        if depmap > 0.7:
            s -= 0.3
        elif depmap > 0.5:
            s -= 0.15
        elif depmap > 0.3:
            s -= 0.05
    gtex_count = getattr(sf, "gtex_tissue_count", None)
    if gtex_count is not None:
        if gtex_count > 40:
            s -= 0.2
        elif gtex_count > 25:
            s -= 0.1
        elif gtex_count > 15:
            s -= 0.05
    gtex_max = getattr(sf, "gtex_max_tpm", None)
    if gtex_max is not None:
        if gtex_max > 1000:
            s -= 0.15
        elif gtex_max > 500:
            s -= 0.08
    # PheWAS (step 14): pleiotropic disease liability from Open Targets disease associations.
    # phewas_hits is stored as JSON {disease_name: association_score} by annotate_genes_phewas().
    phewas = getattr(sf, "phewas_hits", None)
    if phewas and isinstance(phewas, dict):
        high_conf = sum(1 for v in phewas.values() if isinstance(v, (int, float)) and v > 0.3)
        max_assoc = max((v for v in phewas.values() if isinstance(v, (int, float))), default=0.0)
        if high_conf >= 15:
            s -= 0.20   # Highly pleiotropic — many disease associations
        elif high_conf >= 5:
            s -= 0.10   # Moderate pleiotropic risk
        if max_assoc > 0.7:
            s -= 0.10   # Very strong single disease association — on-target liability
    return round(max(s, 0.0), 4)


def regulatory_score(gene_id: str, session) -> float:
    """Score in [0, 1] from Track B AlphaGenome regulatory divergence.

    Aggregates the max promoter divergence across all resilient species,
    boosted by the number of independent lineages with significant signal.
    Returns 0.0 if no RegulatoryDivergence rows exist for this gene.

    Note: This function is for external callers and tests only. The main
    run_scoring() path uses a pre-fetched bulk map (reg_map) to avoid N+1 queries.
    """
    rows = (
        session.query(RegulatoryDivergence)
        .filter(RegulatoryDivergence.gene_id == gene_id)
        .all()
    )
    if not rows:
        return 0.0
    effects = [r.regulatory_score for r in rows if r.regulatory_score is not None]
    if not effects:
        return 0.0
    # lineage_count is stored on each row (same for all rows of a gene)
    lineage_count = max((r.lineage_count or 0) for r in rows)
    max_effect = max(effects)
    return round(min(max_effect + 0.05 * lineage_count, 1.0), 4)


def _compute_nucleotide_divergence_score(
    conservation_score: float,
    regulatory_divergence_count: int,
    regulatory_convergence_count: int,
    promoter_phylo_score: Optional[float],
) -> float:
    """Core nucleotide divergence score computation, reusable from bulk maps or per-gene queries.

    Combines four signals into a single [0, 1] score:
      1. Promoter divergence  (40%) — inverted conservation; lower conservation in resilient species
         implies regulatory divergence.
      2. Divergence count     (30%) — regulatory mutations shared by ≥3 resilient species and
         absent in human + controls. Capped at 10.
      3. PhyloP signal        (20%) — reward accelerated evolution (negative phyloP = faster than
         neutral rate in the promoter).
      4. Convergence bonus    (10%) — same mutation in ≥2 distinct lineage clusters. Capped at 5.
    """
    promoter_divergence = round(1.0 - min(max(conservation_score, 0.0), 1.0), 4)
    div_component = round(min(regulatory_divergence_count / 10.0, 1.0), 4)
    # conv_bonus maps [0, 5] convergence count → [0, 0.2], then normalised to [0, 1] for the
    # weighted slot. The raw bonus is NOT added a second time outside the weighted sum.
    conv_count_norm = round(min(regulatory_convergence_count / 5.0, 1.0), 4)
    phylo_component = 0.0
    if promoter_phylo_score is not None and promoter_phylo_score < 0:
        phylo_component = round(min(-promoter_phylo_score / 3.0, 1.0), 4)

    score = (
        0.40 * promoter_divergence
        + 0.30 * div_component
        + 0.20 * phylo_component
        + 0.10 * conv_count_norm
    )
    return round(min(score, 1.0), 4)


def nucleotide_divergence_score(gene_id: str, session) -> float:
    """Score [0, 1] from NucleotideScore and PhyloConservationScore tables.

    Returns 0.0 if no nucleotide data is available (gracefully absent before step3c).

    Note: This function is for external callers and tests only. The main
    run_scoring() path uses pre-fetched bulk maps via _compute_nucleotide_divergence_score().
    """
    ns_promoter = (
        session.query(NucleotideScore)
        .filter_by(gene_id=gene_id, region_type="promoter")
        .first()
    )
    if ns_promoter is None:
        return 0.0

    pcs = session.get(PhyloConservationScore, gene_id)
    promoter_phylo_score = pcs.promoter_phylo_score if pcs is not None else None

    return _compute_nucleotide_divergence_score(
        conservation_score=ns_promoter.conservation_score or 0.0,
        regulatory_divergence_count=ns_promoter.regulatory_divergence_count or 0,
        regulatory_convergence_count=ns_promoter.regulatory_convergence_count or 0,
        promoter_phylo_score=promoter_phylo_score,
    )


def composite_score(sub_scores: dict[str, float], weights: dict[str, float]) -> float:
    """Weighted sum of sub-scores normalised to [0, 1].

    Only includes a component in the denominator if:
      - Its weight > 0, AND
      - A score is actually available (> 0) OR it is a core evolutionary signal
        (convergence, selection) that should always anchor the composite.

    This prevents Phase 2 disease/druggability weights from diluting the score
    to near-zero when those databases return no annotations (e.g. small test runs).
    """
    # Core signals always counted even if zero (they reflect real absence of signal),
    # EXCEPT expression — when expression_score=0 it typically means no data rather
    # than a genuinely zero score, so exclude it from the denominator in that case
    # to avoid diluting strong evolutionary signal for genes lacking expression data.
    core_keys = {"convergence", "selection"}
    total_weight = 0.0
    weighted_sum = 0.0
    for k, w in weights.items():
        if w <= 0:
            continue
        score = sub_scores.get(k, 0.0)
        if k in core_keys or score > 0:
            total_weight += w
            weighted_sum += score * w
    if total_weight == 0:
        return 0.0
    return round(weighted_sum / total_weight, 4)


def assign_tier(score: float, thresholds: dict[str, float], human_genetics_score: float = 0.0) -> str:
    """Assign tier based on composite score.

    'Validated' is a special designation for genes that have BOTH strong
    evolutionary signal (Tier1/Tier2) AND human genetic evidence (GWAS,
    OpenTargets genetic association, or a rare protective variant matching
    the animal divergence direction in gnomAD).

    This is the highest confidence category for clinical translation.
    """
    if score >= thresholds.get("tier1", 0.70):
        tier = "Tier1"
    elif score >= thresholds.get("tier2", 0.40):
        tier = "Tier2"
    else:
        return "Tier3"
    # Upgrade to Validated if human genetics evidence is present
    if human_genetics_score >= 0.3:
        return "Validated"
    return tier


def human_genetics_score_from_disease(ann: Optional[DiseaseAnnotation]) -> float:
    """Compute a human genetics confidence score from disease annotation data.

    This is NOT the same as disease_score() — it specifically captures whether
    the gene has *genetic* (not just phenotypic) evidence in humans:
    - GWAS association (gwas_pvalue): genetic variant in humans linked to disease
    - gnomAD pLI (gnomad_pli): high pLI means loss of function is intolerant —
      strong signal the gene is functionally critical in humans
    - OpenTargets genetic_evidence (subset of opentargets_score): causal genetic link
    - Rare protective variant matching animal divergence direction (U6):
      the PCSK9 paradigm — a naturally occurring human variant pointing the same
      biochemical direction as the animal adaptation. Immediate 'Validated' upgrade
      if pvalue < 5e-8.

    Returns 0 if no human genetic support; >0.3 qualifies for 'Validated' tier.
    """
    if ann is None:
        return 0.0
    score = 0.0
    # GWAS: direct human genetic evidence
    if ann.gwas_pvalue is not None and 0 < ann.gwas_pvalue < 5e-8:
        score += min(-math.log10(ann.gwas_pvalue) / 15.0, 0.5)
    elif ann.gwas_pvalue is not None and ann.gwas_pvalue < 1e-5:
        score += 0.1
    # gnomAD pLI: functional constraint in human population
    if ann.gnomad_pli is not None:
        score += min(ann.gnomad_pli, 1.0) * 0.3
    # OpenTargets: genetic association evidence
    if ann.opentargets_score is not None and ann.opentargets_score > 0.3:
        score += min(ann.opentargets_score, 1.0) * 0.3
    # U6: Rare protective variant — PCSK9 paradigm
    # A genome-wide significant protective variant in the same region as the animal
    # divergence is the strongest possible human validation signal.
    prot_pvalue = getattr(ann, "protective_variant_pvalue", None)
    prot_count = getattr(ann, "protective_variant_count", None)
    if prot_pvalue is not None and prot_pvalue < 5e-8 and prot_count and prot_count >= 1:
        score += 0.5  # Immediate Validated qualification
    elif prot_pvalue is not None and prot_pvalue < 1e-5 and prot_count and prot_count >= 1:
        score += 0.2
    elif prot_count and prot_count >= 1:
        score += 0.1  # Rare variant present, no GWAS signal yet — still informative
    return round(min(score, 1.0), 4)


# ---------------------------------------------------------------------------
# Rank-product scoring (Phase 1)
# ---------------------------------------------------------------------------

def _rank_values(values: list[float], ascending: bool = True) -> list[int]:
    """Rank values with ties getting the average rank.

    ascending=True: smallest value gets rank 1 (good for p-values).
    ascending=False: largest value gets rank 1 (good for scores).
    """
    n = len(values)
    indexed = list(enumerate(values))
    indexed.sort(key=lambda x: x[1], reverse=not ascending)

    ranks = [0] * n
    i = 0
    while i < n:
        j = i
        while j < n - 1 and indexed[j + 1][1] == indexed[j][1]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = avg_rank
        i = j + 1
    return ranks


def _compute_rank_product(ranks_per_layer: list[list[int]], n_genes: int) -> list[float]:
    """Compute normalized rank product across k evidence layers.

    RP_i = (product of ranks across layers) / n^k
    """
    k = len(ranks_per_layer)
    rp = []
    for i in range(n_genes):
        product = 1.0
        for layer_ranks in ranks_per_layer:
            product *= layer_ranks[i]
        rp.append(product / (n_genes ** k))
    return rp


def _rank_product_pvalues_analytical(rp_observed: list[float], k: int) -> list[float]:
    """Analytical rank-product p-values using the log-normal approximation.

    Under the null hypothesis, each rank r_i is Uniform(1, n), so r_i/n → Uniform(0,1).
    For k independent uniform layers:
        ln(RP) = Σ ln(r_i/n)  ~  N(-k, k)   (CLT; exact for large n)

    Therefore:
        p_i = Φ((ln(RP_i) + k) / √k)

    This gives properly calibrated, continuous p-values for any n and k without
    the 1/n_permutations floor that makes Tier 2 unreachable with n=12k genes.

    Reference: Breitling et al. (2004) FEBS Letters 573:83-92 (log-normal approx);
               Koziol (2010) for exactness of the approximation.
    """
    sqrt_k = math.sqrt(k)
    pvals = []
    for rp in rp_observed:
        if rp <= 0:
            pvals.append(0.0)
            continue
        z = (math.log(rp) + k) / sqrt_k
        # Normal CDF via math.erf — no scipy dependency
        p = (1.0 + math.erf(z / math.sqrt(2.0))) / 2.0
        pvals.append(max(0.0, min(1.0, p)))
    return pvals


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------


def run_scoring(phase: str = "phase1", trait_id: Optional[str] = None) -> None:
    """Compute and persist CandidateScore for all genes with Layer 1/2 data.

    Phase 1: uses rank-product method across three evidence layers:
      1. Selection (MEME p-value)
      2. Convergence (PhyloP score)
      3. Convergent AAs (max count per gene)
    Tiers assigned by FDR-corrected rank-product p-value.

    Phase 2: uses weighted composite (when disease/druggability data available).
    """
    weights = get_scoring_weights(phase)
    tier_thresholds = get_tier_thresholds()
    tid = trait_id if trait_id is not None else ""

    log.info("Assembling composite scores (phase=%s, trait_id=%r)", phase, tid)

    if phase == "phase1":
        _run_scoring_rank_product(tid)
    else:
        _run_scoring_weighted(weights, tier_thresholds, tid)

    apply_bh_to_evolution_scores()
    log.info("Scoring complete.")


def _run_scoring_rank_product(tid: str) -> None:
    """Phase 1 rank-product scoring across evolution + functional evidence layers.

    Evidence layers (added when data is available):
      1. Selection     — PAML branch-site LRT p-value (always present)
      2. Convergence   — PhyloP conservation score
      3. Convergent AAs — max convergent AA count per gene
      4. Functional    — Open Targets + GTEx + DepMap combined score (step 8)
    """
    with get_session() as session:
        genes = session.query(Gene).all()
        gene_ids = [g.id for g in genes]
        n = len(gene_ids)
        log.info("Rank-product scoring for %d genes...", n)

        if n == 0:
            return

        ev_map = {r.gene_id: r for r in session.query(EvolutionScore).all()}

        conv_aa_rows = (
            session.query(
                Ortholog.gene_id,
                func.max(DivergentMotif.convergent_aa_count),
            )
            .join(DivergentMotif, DivergentMotif.ortholog_id == Ortholog.id)
            .group_by(Ortholog.gene_id)
            .all()
        )
        conv_aa_map: dict[str, int] = {gid: cnt or 0 for gid, cnt in conv_aa_rows}

        # Read functional evidence scores written by step 8 (preserved across the update).
        # expression_score is set by pipeline/layer1_sequence/functional_evidence.py;
        # we read it here so it can feed into the rank product.
        existing_cs = {
            r.gene_id: r
            for r in session.query(CandidateScore).filter_by(trait_id=tid).all()
        }
        func_ev_scores: dict[str, float] = {
            gid: float(cs.expression_score)
            for gid, cs in existing_cs.items()
            if cs.expression_score is not None and cs.expression_score > 0
        }

        selection_pvals = []
        convergence_signals = []   # Fix 3+2: convergence_weight or 1−convergence_pval
        convergent_aa_counts = []
        functional_scores = []

        for gid in gene_ids:
            ev = ev_map.get(gid)
            # Accept p-values from both PAML (paml_branch_site) and HyPhy
            # (meme, MEME_episodic, busted_ph, aBSREL) — all write dnds_pvalue.
            # Genes with no recognised model or NULL p-value get pval=1.0 (no signal).
            _VALID_SELECTION_MODELS = {
                "paml_branch_site", "meme", "MEME_episodic", "busted_ph", "aBSREL",
                # Proxy model: computed from sequence divergence when CDS is unavailable.
                # Included so genes without CDS data are not unfairly penalised with pval=1.0.
                "protein_divergence_proxy",
            }
            if (
                ev is not None
                and ev.selection_model in _VALID_SELECTION_MODELS
                and ev.dnds_pvalue is not None
                and 0.0 < ev.dnds_pvalue <= 1.0
            ):
                pval = ev.dnds_pvalue
            else:
                pval = 1.0
            selection_pvals.append(pval)

            # Fix 2+3: use permutation p-value if available; fall back to raw
            # convergence_weight (Fix 3); only fall back to the old phylop_score
            # if neither new column is populated (legacy runs before migration 0023).
            if ev is not None and ev.convergence_pval is not None:
                # Convert p-value to a "higher = more convergent" signal for ranking
                conv_signal = 1.0 - ev.convergence_pval
            elif ev is not None and ev.convergence_weight is not None:
                conv_signal = ev.convergence_weight
            else:
                # Legacy fallback: phylop_score may hold the old convergence weight
                conv_signal = ev.phylop_score if ev and ev.phylop_score is not None else 0.0
            convergence_signals.append(conv_signal)

            aa_count = conv_aa_map.get(gid, 0)
            convergent_aa_counts.append(float(aa_count))

            functional_scores.append(func_ev_scores.get(gid, 0.0))

        sel_ranks  = _rank_values(selection_pvals,      ascending=True)
        conv_ranks = _rank_values(convergence_signals,   ascending=False)
        aa_ranks   = _rank_values(convergent_aa_counts,  ascending=False)
        func_ranks = _rank_values(functional_scores,     ascending=False)

        available_layers = [sel_ranks]
        layer_names = ["selection"]

        has_conv = any(v > 0 for v in convergence_signals)
        if has_conv:
            available_layers.append(conv_ranks)
            layer_names.append("convergence")

        has_aa = any(v > 0 for v in convergent_aa_counts)
        if has_aa:
            available_layers.append(aa_ranks)
            layer_names.append("convergent_aa")

        has_func = any(v > 0 for v in functional_scores)
        if has_func:
            available_layers.append(func_ranks)
            layer_names.append("functional_evidence")

        log.info("Using %d evidence layers: %s", len(available_layers), ", ".join(layer_names))

        rp_values = _compute_rank_product(available_layers, n)

        log.info(
            "Computing analytical rank-product p-values (log-normal, k=%d layers)...",
            len(available_layers),
        )
        rp_pvals = _rank_product_pvalues_analytical(rp_values, k=len(available_layers))

        log.info("Applying BH FDR correction...")
        rp_qvals = apply_bh_correction(rp_pvals)

        cs_map = existing_cs  # reuse already-loaded map

        # Pre-compute per-gene "active layer" counts for the minimum-evidence guard.
        # A gene must have meaningful signal in at least 2 layers to leave Tier3;
        # this prevents a single strong convergence hit from single-handedly pulling
        # a gene into Tier1/2 when all other evidence is absent.
        def _active_layers(idx: int, gid: str) -> int:
            count = 0
            if selection_pvals[idx] < 0.1:
                count += 1
            if convergence_signals[idx] > 0:
                count += 1
            if convergent_aa_counts[idx] > 0:
                count += 1
            if functional_scores[idx] > 0:
                count += 1
            return count

        tier_counts = {"Tier1": 0, "Tier2": 0, "Tier3": 0}
        for idx, gid in enumerate(gene_ids):
            ev = ev_map.get(gid)
            qval = rp_qvals[idx] if rp_qvals[idx] is not None else 1.0

            if qval < 0.05:
                tier = "Tier1"
            elif qval < 0.20:
                tier = "Tier2"
            else:
                tier = "Tier3"

            # Fix 3 — minimum 2-layer evidence guard:
            # Demote to Tier3 if fewer than 2 independent evidence layers show signal.
            # This prevents genes that score well on convergence alone (or any single
            # layer) from appearing in Tier1/2 without corroborating evidence.
            if tier in ("Tier1", "Tier2") and _active_layers(idx, gid) < 2:
                log.debug(
                    "Demoting gene %s from %s to Tier3: only 1 active evidence layer.",
                    gid[:8], tier,
                )
                tier = "Tier3"

            tier_counts[tier] += 1

            composite = round(1.0 - qval, 4)

            conv_score_val = convergence_score(
                ev.convergence_count if ev else 0,
                phylo_weight=ev.convergence_weight if ev and ev.convergence_weight else None,
            )
            sel_score_val = selection_score(
                ev.dnds_ratio if ev else None,
                ev.dnds_pvalue if ev else None,
                fel_sites=ev.fel_sites if ev else None,
                busted_pvalue=ev.busted_pvalue if ev else None,
                relax_k=ev.relax_k if ev else None,
                relax_pvalue=ev.relax_pvalue if ev else None,
                branches_under_selection=ev.branches_under_selection if ev else None,
                selection_model=ev.selection_model if ev else None,
            )

            cs = cs_map.get(gid)
            if cs is None:
                cs = CandidateScore(gene_id=gid, trait_id=tid)
                session.add(cs)
                cs_map[gid] = cs

            cs.convergence_score = conv_score_val
            cs.selection_score = sel_score_val
            # Preserve expression_score written by step 8; do not zero it out.
            if cs.expression_score is None:
                cs.expression_score = 0.0
            cs.disease_score = 0.0
            cs.druggability_score = 0.0
            cs.safety_score = 0.5
            cs.regulatory_score = 0.0
            cs.composite_score = composite
            cs.tier = tier
            cs.updated_at = datetime.now(timezone.utc)

    for tier_name in ["Tier1", "Tier2", "Tier3"]:
        log.info("  %s: %d genes", tier_name, tier_counts[tier_name])


def _load_phenotype_disease_ids(trait_id: str) -> set[str]:
    """Return the set of EFO/MONDO disease IDs for this phenotype from functional_evidence_config.json.

    Used by disease_score() to decide whether a known drug's indication matches the
    current phenotype, avoiding credit for off-phenotype pharmaceutical attention.
    Returns an empty set when the config is missing or the trait is not listed.
    """
    import json as _json
    from pathlib import Path as _Path
    cfg_path = _Path(__file__).parent.parent / "config" / "functional_evidence_config.json"
    try:
        with open(cfg_path) as f:
            cfg = _json.load(f)
        ot_cfg = cfg.get(trait_id, {}).get("opentargets", {})
        return set(ot_cfg.get("disease_ids", []))
    except Exception:
        return set()


def _run_scoring_weighted(weights: dict, tier_thresholds: dict, tid: str) -> None:
    """Phase 2+ weighted composite scoring (original method)."""

    # Load phenotype disease IDs once for phenotype-aware drug phase credit in disease_score().
    phenotype_disease_ids = _load_phenotype_disease_ids(tid) if tid else set()
    if phenotype_disease_ids:
        log.info("disease_score: phenotype-aware drug credit enabled (%d IDs for %r)", len(phenotype_disease_ids), tid)
    else:
        log.info("disease_score: no phenotype disease IDs found for %r — drug credit is phenotype-neutral", tid)

    with get_session() as session:
        genes = session.query(Gene).all()
        log.info("Weighted scoring for %d genes...", len(genes))

        ev_map  = {r.gene_id: r for r in session.query(EvolutionScore).all()}
        ann_map = {r.gene_id: r for r in session.query(DiseaseAnnotation).all()}
        dt_map  = {r.gene_id: r for r in session.query(DrugTarget).all()}
        sf_map  = {r.gene_id: r for r in session.query(SafetyFlag).all()}
        # Structural score pre-computed by Step 9b (compute_gene_structural_score).
        # Read from CandidateScore.structural_score written during step 9b.
        # We use the cs_map built below to avoid a second bulk query here.
        # If the structural weight is zero (phase1) or the column is NULL (step 9b
        # not yet run), the score defaults to 0.0 and is excluded from the weighted
        # composite denominator by composite_score()'s missing-data guard.
        struct_map: dict[str, float] = {}  # populated from cs_map after it is built
        expr_map = {
            r.gene_id: (r.expression_score or 0.0)
            for r in session.query(CandidateScore).filter_by(trait_id="").all()
        }
        reg_rows_all = session.query(RegulatoryDivergence).all()
        reg_map: dict[str, list] = {}
        for r in reg_rows_all:
            reg_map.setdefault(r.gene_id, []).append(r)
        nucl_map = {
            r.gene_id: r
            for r in session.query(NucleotideScore).filter_by(region_type="promoter").all()
        }
        phylo_map = {r.gene_id: r for r in session.query(PhyloConservationScore).all()}
        cs_map = {
            r.gene_id: r
            for r in session.query(CandidateScore).filter_by(trait_id=tid).all()
        }
        # Populate struct_map from the structural_score already stored by step 9b.
        # Defaults to 0.0 for genes where step 9b has not run or found no motifs.
        struct_map = {
            gid: (cs.structural_score or 0.0)
            for gid, cs in cs_map.items()
            if cs.structural_score is not None
        }

        for gene in genes:
            ev  = ev_map.get(gene.id)
            ann = ann_map.get(gene.id)
            dt  = dt_map.get(gene.id)
            sf  = sf_map.get(gene.id)

            # Use the same convergence signal hierarchy as phase 1 scoring:
            # permutation p-value → convergence_weight → phylop_score (legacy fallback)
            if ev and ev.convergence_pval is not None:
                _phylo_w = 1.0 - ev.convergence_pval
            elif ev and ev.convergence_weight is not None:
                _phylo_w = ev.convergence_weight
            else:
                _phylo_w = ev.phylop_score if ev and ev.phylop_score else None
            conv_score_val = convergence_score(
                ev.convergence_count if ev else 0,
                phylo_weight=_phylo_w,
            )
            # Apply control-divergence penalty stored by step9b (compute_control_divergence_fractions).
            # We read the fraction from the *existing* CandidateScore row (pre-loaded in cs_map)
            # before overwriting it below. This ensures run_scoring() is idempotent: if the
            # fraction was stored by step9b, the composite_score correctly reflects it every time.
            # Without this, apply_control_divergence_penalty() modifies convergence_score but
            # the change is overwritten the next time _run_scoring_weighted runs.
            existing_cs = cs_map.get(gene.id)
            if existing_cs is not None and existing_cs.control_divergence_fraction is not None:
                frac = existing_cs.control_divergence_fraction
                if frac > 0.5:
                    penalty = 1.0 - 0.6 * frac
                    conv_score_val = round(max(0.0, conv_score_val * penalty), 4)
                    log.debug(
                        "Control divergence penalty (frac=%.2f → ×%.2f) applied to %s: conv=%.4f",
                        frac, penalty, gene.id, conv_score_val,
                    )
            sel_score_val = selection_score(
                ev.dnds_ratio if ev else None,
                ev.dnds_pvalue if ev else None,
                fel_sites=ev.fel_sites if ev else None,
                busted_pvalue=ev.busted_pvalue if ev else None,
                relax_k=ev.relax_k if ev else None,
                relax_pvalue=ev.relax_pvalue if ev else None,
                branches_under_selection=ev.branches_under_selection if ev else None,
                selection_model=ev.selection_model if ev else None,
            )
            expr_score_val = expr_map.get(gene.id, 0.0)
            dis_score = disease_score(ann, phenotype_disease_ids=phenotype_disease_ids or None) if weights.get("disease", 0) > 0 else 0.0
            drug_score = druggability_score(dt) if weights.get("druggability", 0) > 0 else 0.0
            # Structural score from Step 9b: pre-computed in annotate_convergent_structural_context
            # and stored in CandidateScore.structural_score. Defaults to 0.0 if step 9b has not
            # run (excluded from composite_score denominator when weight > 0 but score == 0.0 by
            # composite_score's missing-data guard for non-core signals).
            struct_score = struct_map.get(gene.id, 0.0) if weights.get("structural", 0) > 0 else 0.0
            # Always compute safety_score even when safety weight is 0.
            # The score is used as a hard multiplicative floor (safety_floor in
            # scoring_weights.json) that zeros out unsafe genes regardless of
            # evolutionary signal.  Never fall back to 0.5 — that would bypass
            # the floor for every gene.
            safe_score = safety_score(sf)
            if weights.get("regulatory", 0) > 0:
                reg_rows_gene = reg_map.get(gene.id, [])
                if reg_rows_gene:
                    effects = [r.regulatory_score for r in reg_rows_gene if r.regulatory_score is not None]
                    lineage_count = max((r.lineage_count or 0) for r in reg_rows_gene)
                    alpha_score = round(min(max(effects) + 0.05 * lineage_count, 1.0), 4) if effects else 0.0
                else:
                    alpha_score = 0.0
                ns = nucl_map.get(gene.id)
                if ns is not None:
                    pcs = phylo_map.get(gene.id)
                    nucl_score = _compute_nucleotide_divergence_score(
                        conservation_score=ns.conservation_score or 0.0,
                        regulatory_divergence_count=ns.regulatory_divergence_count or 0,
                        regulatory_convergence_count=ns.regulatory_convergence_count or 0,
                        promoter_phylo_score=pcs.promoter_phylo_score if pcs is not None else None,
                    )
                else:
                    nucl_score = 0.0
                reg_score = round(0.6 * alpha_score + 0.4 * nucl_score, 4)
            else:
                reg_score = 0.0

            sub_scores = {
                "convergence": conv_score_val,
                "selection": sel_score_val,
                "expression": expr_score_val,
                "disease": dis_score,
                "druggability": drug_score,
                "safety": safe_score,
                "regulatory": reg_score,
                "structural": struct_score,
            }

            comp = composite_score(sub_scores, weights)

            # Safety is a hard floor, not an additive contribution.
            # If safety_score falls below the configured threshold the gene is
            # zeroed out regardless of how strong the evolutionary signal is.
            # This prevents a gene with a fatal safety liability (e.g. cardiac
            # ion channel, ubiquitously essential) from reaching Tier1 on the
            # back of convergence alone.
            safety_floor = weights.get("safety_floor", 0.0)
            if safety_floor > 0 and safe_score < safety_floor:
                comp = 0.0

            hg_score = human_genetics_score_from_disease(ann)
            tier = assign_tier(comp, tier_thresholds, human_genetics_score=hg_score)

            cs = cs_map.get(gene.id)
            if cs is None:
                cs = CandidateScore(gene_id=gene.id, trait_id=tid)
                session.add(cs)
                cs_map[gene.id] = cs

            cs.convergence_score = conv_score_val
            cs.selection_score = sel_score_val
            cs.expression_score = expr_score_val
            cs.disease_score = dis_score
            cs.druggability_score = drug_score
            cs.safety_score = safe_score
            cs.regulatory_score = reg_score
            # Preserve existing structural_score from step 9b — do NOT overwrite it
            # here.  _run_scoring_weighted merely READS the structural_score to fold
            # it into the composite.  The authoritative value is computed by
            # annotate_convergent_structural_context() and must not be reset to 0.
            if weights.get("structural", 0) > 0:
                cs.structural_score = struct_score  # reflects what actually went into composite
            cs.composite_score = comp
            cs.tier = tier
            cs.updated_at = datetime.now(timezone.utc)

    with get_session() as session:
        for tier_name in ["Tier1", "Tier2", "Tier3"]:
            count = session.query(CandidateScore).filter_by(tier=tier_name, trait_id=tid).count()
            log.info("  %s: %d genes", tier_name, count)


def get_top_candidates(n: int = 20, tier: Optional[str] = None, trait_id: Optional[str] = None) -> list[dict]:
    """Return top N candidates sorted by composite_score descending."""
    tid = trait_id if trait_id is not None else ""
    with get_session() as session:
        q = session.query(CandidateScore, Gene).join(Gene, CandidateScore.gene_id == Gene.id).filter(CandidateScore.trait_id == tid)
        if tier:
            q = q.filter(CandidateScore.tier == tier)
        q = q.order_by(CandidateScore.composite_score.desc()).limit(n)

        results = []
        for cs, gene in q:
            results.append({
                "gene_id": gene.id,
                "gene_symbol": gene.gene_symbol,
                "human_protein": gene.human_protein,
                "composite_score": cs.composite_score,
                "tier": cs.tier,
                "convergence_score": cs.convergence_score,
                "selection_score": cs.selection_score,
                "expression_score": cs.expression_score,
                "regulatory_score": cs.regulatory_score,
                "updated_at": cs.updated_at.isoformat() if cs.updated_at else None,
            })
        return results
