"""Step 9 — Composite score assembly.

Reads Layer 1 + Layer 2 data from the database, applies scoring weights,
normalises each sub-score to [0, 1], and writes CandidateScore rows.

Phase 1 weights (from config/scoring_weights.json, phase1 key):
  convergence:  0.40
  selection:    0.35
  expression:   0.25
  disease:      0.00  (zero-weighted until Phase 2)
  druggability: 0.00  (zero-weighted until Phase 2)
  safety:       0.00  (zero-weighted until Phase 2)

Tier assignment:
  Tier 1: composite_score ≥ 0.70
  Tier 2: composite_score ≥ 0.40
  Tier 3: composite_score < 0.40
"""

import logging
import math
from datetime import datetime
from typing import Optional

from sqlalchemy import func

from db.models import (
    CandidateScore,
    DiseaseAnnotation,
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
    """Retrieve expression_score already stored in CandidateScore (from Step 8)."""
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


def nucleotide_divergence_score(gene_id: str, session) -> float:
    """Score [0, 1] from NucleotideScore and PhyloConservationScore tables.

    Combines four signals:
      1. Promoter conservation (inverted — lower conservation in resilient species
         implies regulatory divergence, which is the signal of interest):
           promoter_divergence = 1.0 - conservation_score  (from NucleotideScore)
      2. Regulatory divergence count (mutations shared by ≥3 resilient species,
         absent in human + controls), capped at 10 → [0, 1]:
           div_component = min(divergence_count / 10, 1.0)
      3. Regulatory convergence count (same mutation in ≥2 distinct lineage clusters):
           conv_bonus = min(convergence_count / 5, 0.2)   (up to 0.2 bonus)
      4. PhyloP/PhastCons scores (accelerated evolution = negative phyloP score
         in the promoter, which we reward):
           phylo_component = 0 if promoter_phylo_score ≥ 0 else min(-phylo_score / 3, 1.0)

    Returns 0.0 if no nucleotide data is available (gracefully absent before step3c).
    """
    ns_promoter = (
        session.query(NucleotideScore)
        .filter_by(gene_id=gene_id, region_type="promoter")
        .first()
    )
    if ns_promoter is None:
        return 0.0

    # 1. Promoter divergence (inverted conservation)
    cons = ns_promoter.conservation_score or 0.0
    promoter_divergence = round(1.0 - min(cons, 1.0), 4)

    # 2. Regulatory divergence count
    div_count = ns_promoter.regulatory_divergence_count or 0
    div_component = round(min(div_count / 10.0, 1.0), 4)

    # 3. Regulatory convergence count (bonus capped at 0.2)
    conv_count = ns_promoter.regulatory_convergence_count or 0
    conv_bonus = round(min(conv_count / 5.0, 0.2), 4)

    # 4. phyloP promoter score: reward accelerated evolution (negative score = faster than neutral)
    phylo_component = 0.0
    pcs = session.get(PhyloConservationScore, gene_id)
    if pcs is not None and pcs.promoter_phylo_score is not None:
        pps = pcs.promoter_phylo_score
        if pps < 0:
            phylo_component = round(min(-pps / 3.0, 1.0), 4)

    # Weighted combination: divergence 40%, div_count 30%, phylo 20%, conv_bonus 10%
    score = (
        0.40 * promoter_divergence
        + 0.30 * div_component
        + 0.20 * phylo_component
        + 0.10 * (conv_bonus / 0.2)  # normalise conv_bonus to [0, 1] for weight
    )
    return round(min(score + conv_bonus, 1.0), 4)


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
# Pipeline entry point
# ---------------------------------------------------------------------------


def run_scoring(phase: str = "phase1", trait_id: Optional[str] = None) -> None:
    """Compute and persist CandidateScore for all genes with Layer 1/2 data.

    Args:
        phase: "phase1" or "phase2" — determines which weights to apply.
        trait_id: Optional trait identifier (e.g. "cancer_resistance"). Use "" for default/legacy.
    """
    weights = get_scoring_weights(phase)
    tier_thresholds = get_tier_thresholds()
    tid = trait_id if trait_id is not None else ""

    log.info("Assembling composite scores (phase=%s, trait_id=%r, weights=%s)", phase, tid, weights)

    with get_session() as session:
        genes = session.query(Gene).all()
        log.info("Scoring %d genes...", len(genes))

        # Pre-fetch all related tables in bulk to avoid N+1 queries
        ev_map  = {r.gene_id: r for r in session.query(EvolutionScore).all()}
        ann_map = {r.gene_id: r for r in session.query(DiseaseAnnotation).all()}
        dt_map  = {r.gene_id: r for r in session.query(DrugTarget).all()}
        sf_map  = {r.gene_id: r for r in session.query(SafetyFlag).all()}
        # Expression scores stored in CandidateScore trait_id=""
        expr_map = {
            r.gene_id: (r.expression_score or 0.0)
            for r in session.query(CandidateScore).filter_by(trait_id="").all()
        }
        # Regulatory divergence: max promoter score per gene
        reg_rows_all = session.query(RegulatoryDivergence).all()
        reg_map: dict[str, list] = {}
        for r in reg_rows_all:
            reg_map.setdefault(r.gene_id, []).append(r)
        # Nucleotide scores (promoter)
        nucl_map = {
            r.gene_id: r
            for r in session.query(NucleotideScore).filter_by(region_type="promoter").all()
        }
        # PhyloConservationScore
        phylo_map = {r.gene_id: r for r in session.query(PhyloConservationScore).all()}
        # Pre-fetch existing CandidateScore rows for this trait to avoid per-gene SELECT
        cs_map = {
            r.gene_id: r
            for r in session.query(CandidateScore).filter_by(trait_id=tid).all()
        }

        for gene in genes:
            ev  = ev_map.get(gene.id)
            ann = ann_map.get(gene.id)
            dt  = dt_map.get(gene.id)
            sf  = sf_map.get(gene.id)

            conv_score = convergence_score(
                ev.convergence_count if ev else 0,
                phylo_weight=ev.phylop_score if ev and ev.phylop_score else None,
            )
            sel_score = selection_score(
                ev.dnds_ratio if ev else None,
                ev.dnds_pvalue if ev else None,
                fel_sites=ev.fel_sites if ev else None,
                busted_pvalue=ev.busted_pvalue if ev else None,
                relax_k=ev.relax_k if ev else None,
                relax_pvalue=ev.relax_pvalue if ev else None,
            )
            expr_score = expr_map.get(gene.id, 0.0)
            dis_score = disease_score(ann) if weights.get("disease", 0) > 0 else 0.0
            drug_score = druggability_score(dt) if weights.get("druggability", 0) > 0 else 0.0
            safe_score = safety_score(sf) if weights.get("safety", 0) > 0 else 1.0
            if weights.get("regulatory", 0) > 0:
                # Regulatory score from pre-fetched maps
                reg_rows_gene = reg_map.get(gene.id, [])
                if reg_rows_gene:
                    effects = [r.regulatory_score for r in reg_rows_gene if r.regulatory_score is not None]
                    lineage_count = max((r.lineage_count or 0) for r in reg_rows_gene)
                    alpha_score = round(min(max(effects) + 0.05 * lineage_count, 1.0), 4) if effects else 0.0
                else:
                    alpha_score = 0.0
                # Nucleotide divergence score from pre-fetched maps
                ns = nucl_map.get(gene.id)
                if ns is not None:
                    cons = ns.conservation_score or 0.0
                    promoter_divergence = round(1.0 - min(cons, 1.0), 4)
                    div_component = round(min((ns.regulatory_divergence_count or 0) / 10.0, 1.0), 4)
                    conv_bonus = round(min((ns.regulatory_convergence_count or 0) / 5.0, 0.2), 4)
                    pcs = phylo_map.get(gene.id)
                    phylo_component = 0.0
                    if pcs is not None and pcs.promoter_phylo_score is not None:
                        pps = pcs.promoter_phylo_score
                        if pps < 0:
                            phylo_component = round(min(-pps / 3.0, 1.0), 4)
                    nucl_score = round(min(
                        0.4 * promoter_divergence + 0.3 * div_component
                        + 0.2 * phylo_component + conv_bonus, 1.0
                    ), 4)
                else:
                    nucl_score = 0.0
                reg_score = round(0.6 * alpha_score + 0.4 * nucl_score, 4)
            else:
                reg_score = 0.0

            sub_scores = {
                "convergence": conv_score,
                "selection": sel_score,
                "expression": expr_score,
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

            cs.convergence_score = conv_score
            cs.selection_score = sel_score
            cs.expression_score = expr_score
            cs.disease_score = dis_score
            cs.druggability_score = drug_score
            cs.safety_score = safe_score
            cs.regulatory_score = reg_score
            cs.composite_score = comp
            cs.tier = tier
            cs.updated_at = datetime.utcnow()

    # Summary
    with get_session() as session:
        for tier_name in ["Tier1", "Tier2", "Tier3"]:
            count = session.query(CandidateScore).filter_by(tier=tier_name, trait_id=tid).count()
            log.info("  %s: %d genes", tier_name, count)

    # Apply global BH FDR correction across all EvolutionScore p-values
    from pipeline.stats import apply_bh_to_evolution_scores
    apply_bh_to_evolution_scores()

    log.info("Scoring complete.")


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
