"""GET /scores/{gene_id} — full score breakdown for a gene."""

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel
from typing import Optional

from db.models import (
    CandidateScore,
    DiseaseAnnotation,
    EvolutionScore,
    ExpressionResult,
    Gene,
    NucleotideScore,
    PhyloConservationScore,
)
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
    # U4: FEL + BUSTED supplementary selection
    fel_sites: Optional[int] = None
    busted_pvalue: Optional[float] = None
    # U6: Rare protective variants
    protective_variant_count: Optional[int] = None
    best_protective_trait: Optional[str] = None
    protective_variant_pvalue: Optional[float] = None
    # Item 10: Literature validation
    lit_score: Optional[float] = None
    lit_pmid_count: Optional[int] = None
    # DNA conservation (steps 3c / 3d)
    nucleotide_cds_conservation:      Optional[float] = None
    promoter_conservation:            Optional[float] = None
    downstream_conservation:          Optional[float] = None
    regulatory_divergence_count:      Optional[int]   = None
    regulatory_convergence_count:     Optional[int]   = None
    cds_phylo_score:                  Optional[float] = None
    promoter_phylo_score:             Optional[float] = None
    downstream_phylo_score:           Optional[float] = None


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
                "fel_sites": getattr(ev, "fel_sites", None),
                "busted_pvalue": getattr(ev, "busted_pvalue", None),
            }

        da = session.get(DiseaseAnnotation, gene_id)

        # DNA conservation data (steps 3c / 3d) — may be None if steps haven't run
        ns_cds        = session.query(NucleotideScore).filter_by(gene_id=gene_id, region_type="cds").first()
        ns_promoter   = session.query(NucleotideScore).filter_by(gene_id=gene_id, region_type="promoter").first()
        ns_downstream = session.query(NucleotideScore).filter_by(gene_id=gene_id, region_type="downstream").first()
        pcs           = session.get(PhyloConservationScore, gene_id)

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
            fel_sites=getattr(ev, "fel_sites", None) if ev else None,
            busted_pvalue=getattr(ev, "busted_pvalue", None) if ev else None,
            protective_variant_count=getattr(da, "protective_variant_count", None) if da else None,
            best_protective_trait=getattr(da, "best_protective_trait", None) if da else None,
            protective_variant_pvalue=getattr(da, "protective_variant_pvalue", None) if da else None,
            lit_score=getattr(da, "lit_score", None) if da else None,
            lit_pmid_count=getattr(da, "lit_pmid_count", None) if da else None,
            # DNA conservation
            nucleotide_cds_conservation=ns_cds.conservation_score if ns_cds else None,
            promoter_conservation=ns_promoter.conservation_score if ns_promoter else None,
            downstream_conservation=ns_downstream.conservation_score if ns_downstream else None,
            regulatory_divergence_count=ns_promoter.regulatory_divergence_count if ns_promoter else None,
            regulatory_convergence_count=ns_promoter.regulatory_convergence_count if ns_promoter else None,
            cds_phylo_score=pcs.cds_phylo_score if pcs else None,
            promoter_phylo_score=pcs.promoter_phylo_score if pcs else None,
            downstream_phylo_score=pcs.downstream_phylo_score if pcs else None,
        )
