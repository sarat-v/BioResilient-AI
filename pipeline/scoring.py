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
    Ortholog,
    RegulatoryDivergence,
    SafetyFlag,
)
from db.session import get_session
from pipeline.config import get_scoring_weights, get_tier_thresholds, get_thresholds

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sub-score computation
# ---------------------------------------------------------------------------


def convergence_score(convergence_count: int, max_lineages: int = 6) -> float:
    """Score in [0, 1] based on number of independent lineages with the motif.

    Sigmoid-shaped: 0 lineages → 0.0, 3 lineages → ~0.73, 6 lineages → ~0.98.
    """
    if convergence_count is None or convergence_count <= 0:
        return 0.0
    # Logistic function centred at 3 lineages
    return round(1.0 / (1.0 + math.exp(-1.5 * (convergence_count - 2))), 4)


def selection_score(dnds_ratio: Optional[float], dnds_pvalue: Optional[float]) -> float:
    """Score in [0, 1] combining dN/dS ratio and statistical significance.

    High score = strong positive selection (dN/dS >> 1) with low p-value.
    Genes under strong purifying selection (dN/dS << 1) score near 0.
    """
    if dnds_ratio is None or dnds_pvalue is None:
        return 0.0

    # Normalise dN/dS: cap at 5, score = min(dnds/5, 1.0)
    dnds_norm = min(dnds_ratio / 5.0, 1.0)

    # p-value weight: -log10(p) normalised to [0, 1], capped at -log10(1e-10)
    if dnds_pvalue <= 0:
        pval_weight = 1.0
    elif dnds_pvalue >= 1:
        pval_weight = 0.0
    else:
        pval_weight = min(-math.log10(dnds_pvalue) / 10.0, 1.0)

    return round(0.5 * dnds_norm + 0.5 * pval_weight, 4)


def expression_score_from_db(gene_id: str, session) -> float:
    """Retrieve expression_score already stored in CandidateScore (from Step 8)."""
    cs = session.get(CandidateScore, gene_id)
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
    """Score in [0, 1] from pockets, ChEMBL, CanSAR tier."""
    if dt is None:
        return 0.0
    s = 0.0
    if dt.pocket_count is not None and dt.pocket_count > 0:
        s += min(dt.pocket_count / 5.0, 1.0) * 0.3
    if dt.top_pocket_score is not None:
        s += min(dt.top_pocket_score, 1.0) * 0.3
    if dt.chembl_target_id:
        s += 0.2
    if dt.existing_drugs and len(dt.existing_drugs) > 0:
        s += 0.2
    tier = (dt.druggability_tier or "").upper()
    if tier == "A":
        s += 0.2
    elif tier == "B":
        s += 0.1
    return round(min(s, 1.0), 4)


def safety_score(sf: Optional[SafetyFlag]) -> float:
    """Score in [0, 1]; high hub_risk, essential, large family reduce score (inverted for composite)."""
    if sf is None:
        return 1.0
    s = 1.0
    if sf.hub_risk:
        s -= 0.3
    if sf.is_essential:
        s -= 0.2
    if sf.family_size is not None and sf.family_size > 100:
        s -= 0.2
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


def assign_tier(score: float, thresholds: dict[str, float]) -> str:
    if score >= thresholds.get("tier1", 0.70):
        return "Tier1"
    if score >= thresholds.get("tier2", 0.40):
        return "Tier2"
    return "Tier3"


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------


def run_scoring(phase: str = "phase1") -> None:
    """Compute and persist CandidateScore for all genes with Layer 1/2 data.

    Args:
        phase: "phase1" or "phase2" — determines which weights to apply.
    """
    weights = get_scoring_weights(phase)
    tier_thresholds = get_tier_thresholds()

    log.info("Assembling composite scores (phase=%s, weights=%s)", phase, weights)

    with get_session() as session:
        genes = session.query(Gene).all()
        log.info("Scoring %d genes...", len(genes))

        for gene in genes:
            ev = session.get(EvolutionScore, gene.id)
            ann = session.get(DiseaseAnnotation, gene.id)
            dt = session.get(DrugTarget, gene.id)
            sf = session.get(SafetyFlag, gene.id)

            conv_score = convergence_score(ev.convergence_count if ev else 0)
            sel_score = selection_score(
                ev.dnds_ratio if ev else None,
                ev.dnds_pvalue if ev else None,
            )
            expr_score = expression_score_from_db(gene.id, session)
            dis_score = disease_score(ann) if weights.get("disease", 0) > 0 else 0.0
            drug_score = druggability_score(dt) if weights.get("druggability", 0) > 0 else 0.0
            safe_score = safety_score(sf) if weights.get("safety", 0) > 0 else 1.0
            reg_score = regulatory_score(gene.id, session) if weights.get("regulatory", 0) > 0 else 0.0

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
            tier = assign_tier(comp, tier_thresholds)

            cs = session.get(CandidateScore, gene.id)
            if cs is None:
                cs = CandidateScore(gene_id=gene.id)
                session.add(cs)

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
            count = session.query(CandidateScore).filter_by(tier=tier_name).count()
            log.info("  %s: %d genes", tier_name, count)

    log.info("Scoring complete.")


def get_top_candidates(n: int = 20, tier: Optional[str] = None) -> list[dict]:
    """Return top N candidates sorted by composite_score descending."""
    with get_session() as session:
        q = session.query(CandidateScore, Gene).join(Gene, CandidateScore.gene_id == Gene.id)
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
