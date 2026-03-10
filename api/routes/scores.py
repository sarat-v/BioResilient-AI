"""GET /scores/{gene_id} — full score breakdown for a gene."""

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel
from typing import Optional

from db.models import CandidateScore, EvolutionScore, ExpressionResult, Gene
from db.session import get_session

router = APIRouter()


class ExpressionEvidenceOut(BaseModel):
    geo_accession: str
    comparison: str | None
    log2fc: float | None
    padj: float | None
    n_samples: int | None


class ScoreBreakdown(BaseModel):
    gene_id: str
    gene_symbol: str
    tier: str
    composite_score: float
    sub_scores: dict[str, float]
    evolution: dict
    expression_evidence: list[ExpressionEvidenceOut] = []


@router.get("/{gene_id}", response_model=ScoreBreakdown)
def get_scores(gene_id: str, trait_id: Optional[str] = Query(None, description="Trait preset id. Omit for default.")):
    """Full score breakdown for a gene, including sub-scores and evolution data."""
    tid = trait_id if trait_id is not None else ""
    with get_session() as session:
        gene = session.get(Gene, gene_id)
        if gene is None:
            raise HTTPException(status_code=404, detail=f"Gene '{gene_id}' not found.")

        cs = session.query(CandidateScore).filter_by(gene_id=gene_id, trait_id=tid).first()
        if cs is None:
            raise HTTPException(status_code=404, detail=f"No scores computed for gene '{gene_id}'.")

        ev = session.get(EvolutionScore, gene_id)
        evolution_data = {}
        if ev:
            evolution_data = {
                "dnds_ratio": ev.dnds_ratio,
                "dnds_pvalue": ev.dnds_pvalue,
                "selection_model": ev.selection_model,
                "branches_under_selection": ev.branches_under_selection or [],
                "convergence_count": ev.convergence_count,
                "phylop_score": ev.phylop_score,
            }

        expr_evidence = []
        for er in session.query(ExpressionResult).filter_by(gene_id=gene_id).limit(50).all():
            expr_evidence.append(ExpressionEvidenceOut(
                geo_accession=er.geo_accession,
                comparison=er.comparison,
                log2fc=er.log2fc,
                padj=er.padj,
                n_samples=er.n_samples,
            ))

        return ScoreBreakdown(
            gene_id=gene.id,
            gene_symbol=gene.gene_symbol,
            tier=cs.tier or "Tier3",
            composite_score=cs.composite_score or 0.0,
            sub_scores={
                "convergence": cs.convergence_score or 0.0,
                "selection": cs.selection_score or 0.0,
                "expression": cs.expression_score or 0.0,
                "disease": cs.disease_score or 0.0,
                "druggability": cs.druggability_score or 0.0,
                "safety": cs.safety_score or 0.0,
                "regulatory": cs.regulatory_score or 0.0,
            },
            evolution=evolution_data,
            expression_evidence=expr_evidence,
        )
