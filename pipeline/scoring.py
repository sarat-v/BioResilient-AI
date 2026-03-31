"""Step 9 — Composite score assembly.

Phase 1 uses a rank-product method instead of arbitrary weighted composites:

Evidence layers (Phase 1):
  1. Selection:  MEME p-value (lower = stronger positive selection)
  2. Convergence: PhyloP score (higher = more convergent evolution)
  3. Convergent AAs: max convergent_aa_count per gene (higher = more convergent)

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
import random
from collections import defaultdict
from datetime import datetime, timezone
from typing import Optional

from sqlalchemy import func

from db.models import (
    CandidateScore,
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
)
from db.session import get_session
from pipeline.config import get_scoring_weights, get_tier_thresholds, get_thresholds
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
) -> float:
    """Score in [0, 1] combining MEME dN/dS, FEL, BUSTED, and RELAX signals.

    Combined formula (when all tests present):
        selection_score = 0.40 × meme + 0.25 × fel + 0.20 × busted + 0.15 × relax

    RELAX component rewards significant branch-specific acceleration (k > 1):
        relax_score = min(-log10(p) / 10, 1.0) × min((k - 1) / 2.0, 1.0)  if k > 1
    """
    # --- MEME component ---
    if dnds_ratio is None or dnds_pvalue is None:
        meme_score = 0.0
    else:
        dnds_norm = min(dnds_ratio / 5.0, 1.0)
        if dnds_pvalue <= 0:
            pval_weight = 1.0
        elif dnds_pvalue >= 1:
            pval_weight = 0.0
        else:
            pval_weight = min(-math.log10(dnds_pvalue) / 10.0, 1.0)
        meme_score = round(0.5 * dnds_norm + 0.5 * pval_weight, 4)

    # If no supplementary tests available, return MEME score alone
    if fel_sites is None and busted_pvalue is None and relax_k is None:
        return meme_score

    # --- FEL component ---
    fel_score = round(min(fel_sites / 10.0, 1.0), 4) if fel_sites and fel_sites > 0 else 0.0

    # --- BUSTED component ---
    if busted_pvalue is not None and 0 < busted_pvalue <= 1:
        busted_score = round(min(-math.log10(busted_pvalue) / 10.0, 1.0), 4)
    else:
        busted_score = 0.0

    # --- RELAX component: intensification (k>1) with p-value significance ---
    if relax_k is not None and relax_k > 1.0 and relax_pvalue is not None and 0 < relax_pvalue <= 1:
        k_component = min((relax_k - 1.0) / 2.0, 1.0)
        p_component = min(-math.log10(relax_pvalue) / 10.0, 1.0)
        relax_score = round(k_component * p_component, 4)
    else:
        relax_score = 0.0

    return round(0.40 * meme_score + 0.25 * fel_score + 0.20 * busted_score + 0.15 * relax_score, 4)


def expression_score_from_db(gene_id: str, session) -> float:
    """Retrieve expression_score already stored in CandidateScore (from Step 8).

    Note: This function is for external callers and tests only. The main
    run_scoring() path uses a pre-fetched bulk map (expr_map) to avoid N+1 queries.
    """
    cs = session.query(CandidateScore).filter_by(gene_id=gene_id, trait_id="").first()
    if cs is None:
        return 0.0
    return cs.expression_score or 0.0


def disease_score(ann: Optional[DiseaseAnnotation]) -> float:
    """Score in [0, 1] from OpenTargets, GWAS p-value, gnomAD pLI (disease relevance)."""
    if ann is None:
        return 0.0
    s = 0.0
    if ann.opentargets_score is not None:
        s += min(ann.opentargets_score, 1.0) * 0.4
    if ann.gwas_pvalue is not None and ann.gwas_pvalue > 0:
        s += min(-math.log10(ann.gwas_pvalue) / 15.0, 1.0) * 0.3
    if ann.gnomad_pli is not None:
        s += min(ann.gnomad_pli, 1.0) * 0.3
    return round(min(s, 1.0), 4)


def druggability_score(dt: Optional[DrugTarget]) -> float:
    """Score in [0, 1] from pockets, ChEMBL, CanSAR tier, and P2Rank ML prediction."""
    if dt is None:
        return 0.0
    s = 0.0
    if dt.pocket_count is not None and dt.pocket_count > 0:
        s += min(dt.pocket_count / 5.0, 1.0) * 0.25
    if dt.top_pocket_score is not None:
        s += min(dt.top_pocket_score, 1.0) * 0.25
    if dt.chembl_target_id:
        s += 0.2
    if dt.existing_drugs and len(dt.existing_drugs) > 0:
        s += 0.2
    tier = (dt.druggability_tier or "").upper()
    if tier == "A":
        s += 0.2
    elif tier == "B":
        s += 0.1
    # P2Rank ML prediction — supplementary signal (up to 0.1 bonus)
    p2rank = getattr(dt, "p2rank_score", None)
    if p2rank is not None and p2rank > 0:
        s += min(p2rank, 1.0) * 0.1
    return round(min(s, 1.0), 4)


def safety_score(sf: Optional[SafetyFlag]) -> float:
    """Score in [0, 1]; high hub_risk, essential, large family reduce score (inverted for composite).

    Extended with DepMap CRISPR essentiality and GTEx expression breadth:
      - depmap_score < -0.5: broadly essential across cell lines → significant penalty
      - gtex_tissue_count > 30: ubiquitous expression → moderate off-target risk
      - gtex_max_tpm > 5000: extreme expression in any tissue → modest penalty
    """
    if sf is None:
        return 1.0
    s = 1.0
    if sf.hub_risk:
        s -= 0.3
    if sf.is_essential:
        s -= 0.2
    if sf.family_size is not None and sf.family_size > 100:
        s -= 0.2
    # DepMap: broad essentiality is a red flag for drug targeting
    depmap = getattr(sf, "depmap_score", None)
    if depmap is not None:
        if depmap < -0.7:
            s -= 0.3   # Very broadly essential — serious toxicity concern
        elif depmap < -0.5:
            s -= 0.15  # Broadly essential
        elif depmap < -0.3:
            s -= 0.05  # Moderate essentiality — small penalty
    # GTEx: ubiquitous expression → harder to target without side effects
    gtex_count = getattr(sf, "gtex_tissue_count", None)
    if gtex_count is not None:
        if gtex_count > 40:
            s -= 0.2
        elif gtex_count > 25:
            s -= 0.1
        elif gtex_count > 15:
            s -= 0.05
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
    # Core signals always counted even if zero (they reflect real absence of signal)
    core_keys = {"convergence", "selection", "expression"}
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
    if ann.gwas_pvalue is not None and ann.gwas_pvalue < 5e-8:
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


def _rank_product_pvalues(rp_observed: list[float], ranks_per_layer: list[list[int]],
                          n_genes: int, n_permutations: int = 1000) -> list[float]:
    """Estimate rank-product p-values via permutation.

    For each permutation, shuffle the ranks independently in each layer,
    compute the RP, and count how often the permuted RP is <= the observed RP.
    """
    k = len(ranks_per_layer)
    rng = random.Random(42)
    counts = [0] * n_genes

    for _ in range(n_permutations):
        perm_layers = []
        for layer_ranks in ranks_per_layer:
            shuffled = list(layer_ranks)
            rng.shuffle(shuffled)
            perm_layers.append(shuffled)

        for i in range(n_genes):
            perm_product = 1.0
            for layer in perm_layers:
                perm_product *= layer[i]
            perm_rp = perm_product / (n_genes ** k)
            if perm_rp <= rp_observed[i]:
                counts[i] += 1

    return [c / n_permutations for c in counts]


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
    """Phase 1 rank-product scoring across evolution evidence layers."""

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

        selection_pvals = []
        convergence_phylop = []
        convergent_aa_counts = []

        for gid in gene_ids:
            ev = ev_map.get(gid)
            pval = ev.dnds_pvalue if ev and ev.dnds_pvalue is not None else 1.0
            selection_pvals.append(pval)

            phylop = ev.phylop_score if ev and ev.phylop_score is not None else 0.0
            convergence_phylop.append(phylop)

            aa_count = conv_aa_map.get(gid, 0)
            convergent_aa_counts.append(float(aa_count))

        sel_ranks = _rank_values(selection_pvals, ascending=True)
        conv_ranks = _rank_values(convergence_phylop, ascending=False)
        aa_ranks = _rank_values(convergent_aa_counts, ascending=False)

        available_layers = [sel_ranks]
        layer_names = ["selection"]

        has_conv = any(v > 0 for v in convergence_phylop)
        if has_conv:
            available_layers.append(conv_ranks)
            layer_names.append("convergence")

        has_aa = any(v > 0 for v in convergent_aa_counts)
        if has_aa:
            available_layers.append(aa_ranks)
            layer_names.append("convergent_aa")

        log.info("Using %d evidence layers: %s", len(available_layers), ", ".join(layer_names))

        rp_values = _compute_rank_product(available_layers, n)

        log.info("Computing rank-product p-values (1000 permutations)...")
        rp_pvals = _rank_product_pvalues(rp_values, available_layers, n, n_permutations=1000)

        log.info("Applying BH FDR correction...")
        rp_qvals = apply_bh_correction(rp_pvals)

        cs_map = {
            r.gene_id: r
            for r in session.query(CandidateScore).filter_by(trait_id=tid).all()
        }

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
            tier_counts[tier] += 1

            composite = round(1.0 - qval, 4)

            conv_score_val = convergence_score(
                ev.convergence_count if ev else 0,
                phylo_weight=ev.phylop_score if ev and ev.phylop_score else None,
            )
            sel_score_val = selection_score(
                ev.dnds_ratio if ev else None,
                ev.dnds_pvalue if ev else None,
                fel_sites=ev.fel_sites if ev else None,
                busted_pvalue=ev.busted_pvalue if ev else None,
                relax_k=ev.relax_k if ev else None,
                relax_pvalue=ev.relax_pvalue if ev else None,
            )

            cs = cs_map.get(gid)
            if cs is None:
                cs = CandidateScore(gene_id=gid, trait_id=tid)
                session.add(cs)
                cs_map[gid] = cs

            cs.convergence_score = conv_score_val
            cs.selection_score = sel_score_val
            cs.expression_score = 0.0
            cs.disease_score = 0.0
            cs.druggability_score = 0.0
            cs.safety_score = 1.0
            cs.regulatory_score = 0.0
            cs.composite_score = composite
            cs.tier = tier
            cs.updated_at = datetime.now(timezone.utc)

    for tier_name in ["Tier1", "Tier2", "Tier3"]:
        log.info("  %s: %d genes", tier_name, tier_counts[tier_name])


def _run_scoring_weighted(weights: dict, tier_thresholds: dict, tid: str) -> None:
    """Phase 2+ weighted composite scoring (original method)."""

    with get_session() as session:
        genes = session.query(Gene).all()
        log.info("Weighted scoring for %d genes...", len(genes))

        ev_map  = {r.gene_id: r for r in session.query(EvolutionScore).all()}
        ann_map = {r.gene_id: r for r in session.query(DiseaseAnnotation).all()}
        dt_map  = {r.gene_id: r for r in session.query(DrugTarget).all()}
        sf_map  = {r.gene_id: r for r in session.query(SafetyFlag).all()}
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

        for gene in genes:
            ev  = ev_map.get(gene.id)
            ann = ann_map.get(gene.id)
            dt  = dt_map.get(gene.id)
            sf  = sf_map.get(gene.id)

            conv_score_val = convergence_score(
                ev.convergence_count if ev else 0,
                phylo_weight=ev.phylop_score if ev and ev.phylop_score else None,
            )
            sel_score_val = selection_score(
                ev.dnds_ratio if ev else None,
                ev.dnds_pvalue if ev else None,
                fel_sites=ev.fel_sites if ev else None,
                busted_pvalue=ev.busted_pvalue if ev else None,
                relax_k=ev.relax_k if ev else None,
                relax_pvalue=ev.relax_pvalue if ev else None,
            )
            expr_score_val = expr_map.get(gene.id, 0.0)
            dis_score = disease_score(ann) if weights.get("disease", 0) > 0 else 0.0
            drug_score = druggability_score(dt) if weights.get("druggability", 0) > 0 else 0.0
            safe_score = safety_score(sf) if weights.get("safety", 0) > 0 else 1.0
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
            }

            comp = composite_score(sub_scores, weights)
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
