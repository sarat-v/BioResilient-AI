"""Step 17 — Preclinical readiness scoring.

Synthesises data from all previous pipeline steps into a single readiness
assessment per Tier1/Tier2 gene with four components:

  synthesis_score      — Can we synthesise / deliver a therapeutic?
  deliverability_score — How tractable is the target?
  model_score          — Can we validate in a preclinical model?
  assay_score          — Do validated assays already exist?

Readiness tiers: Ready (≥0.7) / Needs Work (0.4–0.7) / Exploratory (<0.4)
"""

import logging
from typing import Optional

from db.models import (
    CandidateScore, DiseaseAnnotation, DrugTarget, GeneTherapyScore,
    PreclinicalReadiness, SafetyFlag,
)
from db.session import get_session

log = logging.getLogger(__name__)

_READY_THRESHOLD    = 0.70
_NEEDS_WORK_THRESHOLD = 0.40


def _synthesis_score(
    therapy: Optional[GeneTherapyScore],
) -> float:
    """0-1 score based on gene therapy deliverability."""
    if therapy is None:
        return 0.3  # unknown
    score = 0.0
    if therapy.aav_compatible:
        score += 0.5
    if therapy.gene_size_bp and therapy.gene_size_bp < 4700:
        score += 0.3
    if therapy.crispr_sites and therapy.crispr_sites > 3:
        score += 0.2
    return min(round(score, 3), 1.0)


def _deliverability_score(
    drug_target: Optional[DrugTarget],
) -> float:
    """0-1 score based on druggability tier and tractability."""
    if drug_target is None:
        return 0.1
    tier_scores = {"A": 1.0, "B": 0.7, "C": 0.4, "undruggable": 0.05}
    base = tier_scores.get(drug_target.druggability_tier or "C", 0.4)

    bonus = 0.0
    if drug_target.tractability_sm:
        bonus += 0.1
    if drug_target.tractability_ab:
        bonus += 0.1
    if drug_target.tractability_protac:
        bonus += 0.05

    return min(round(base + bonus, 3), 1.0)


def _model_score(
    safety: Optional[SafetyFlag],
    disease: Optional[DiseaseAnnotation],
) -> float:
    """0-1 score for preclinical model availability."""
    score = 0.2  # baseline
    if disease and disease.mouse_ko_phenotype:
        score += 0.4
    if safety and safety.depmap_score is not None:
        score += 0.1   # at least has cell line data
        # Non-essential (< 0.5 prob) = safe to KO = better model
        if safety.depmap_score < 0.5:
            score += 0.3
    return min(round(score, 3), 1.0)


def _assay_score(
    drug_target: Optional[DrugTarget],
    disease: Optional[DiseaseAnnotation],
) -> float:
    """0-1 score based on existing assay chemistry."""
    score = 0.0
    existing_drugs = (drug_target.existing_drugs or []) if drug_target else []
    if existing_drugs:
        score += min(len(existing_drugs) * 0.15, 0.6)
    if disease and disease.known_drug_phase:
        score += min(disease.known_drug_phase * 0.1, 0.4)
    return min(round(score, 3), 1.0)


def _readiness_tier(score: float) -> str:
    if score >= _READY_THRESHOLD:
        return "Ready"
    if score >= _NEEDS_WORK_THRESHOLD:
        return "Needs Work"
    return "Exploratory"


def score_preclinical_readiness(gene_ids: list[str]) -> int:
    """Compute and store PreclinicalReadiness rows for given genes."""
    weights = {
        "synthesis":     0.25,
        "deliverability": 0.35,
        "model":         0.25,
        "assay":         0.15,
    }

    with get_session() as session:
        therapy_map = {r.gene_id: r for r in
                       session.query(GeneTherapyScore).filter(
                           GeneTherapyScore.gene_id.in_(gene_ids)).all()}
        drug_map = {r.gene_id: r for r in
                    session.query(DrugTarget).filter(
                        DrugTarget.gene_id.in_(gene_ids)).all()}
        safety_map = {r.gene_id: r for r in
                      session.query(SafetyFlag).filter(
                          SafetyFlag.gene_id.in_(gene_ids)).all()}
        disease_map = {r.gene_id: r for r in
                       session.query(DiseaseAnnotation).filter(
                           DiseaseAnnotation.gene_id.in_(gene_ids)).all()}

    updated = 0
    for gid in gene_ids:
        s_syn  = _synthesis_score(therapy_map.get(gid))
        s_del  = _deliverability_score(drug_map.get(gid))
        s_mod  = _model_score(safety_map.get(gid), disease_map.get(gid))
        s_ass  = _assay_score(drug_map.get(gid), disease_map.get(gid))

        overall = round(
            s_syn  * weights["synthesis"]
            + s_del * weights["deliverability"]
            + s_mod * weights["model"]
            + s_ass * weights["assay"],
            4,
        )
        tier = _readiness_tier(overall)

        notes_parts = []
        if s_syn < 0.3:
            notes_parts.append("gene delivery needs verification")
        if s_del < 0.3:
            notes_parts.append("target tractability low")
        if s_mod < 0.3:
            notes_parts.append("no preclinical model data")
        if s_ass < 0.3:
            notes_parts.append("no validated assay chemistry")

        log.info("Preclinical %s: syn=%.2f del=%.2f mod=%.2f ass=%.2f → %.2f (%s)",
                 gid, s_syn, s_del, s_mod, s_ass, overall, tier)

        with get_session() as session:
            row = session.get(PreclinicalReadiness, gid)
            if row is None:
                row = PreclinicalReadiness(gene_id=gid)
                session.add(row)
            row.synthesis_score      = s_syn
            row.deliverability_score = s_del
            row.model_score          = s_mod
            row.assay_score          = s_ass
            row.overall_readiness_score = overall
            row.readiness_tier       = tier
            row.notes = "; ".join(notes_parts) if notes_parts else None
            session.commit()

        updated += 1

    log.info("Preclinical readiness: %d genes scored.", updated)
    return updated


def run_preclinical_pipeline(gene_ids: list[str]) -> int:
    """Entry point called from orchestrator step17."""
    return score_preclinical_readiness(gene_ids)
