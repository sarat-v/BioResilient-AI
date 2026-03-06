"""GET /candidates and GET /candidates/{gene_id}."""

from typing import Optional

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel

from db.models import CandidateScore, DivergentMotif, EvolutionScore, Gene, Ortholog
from db.session import get_session

router = APIRouter()


# ---------------------------------------------------------------------------
# Response models
# ---------------------------------------------------------------------------


class MotifOut(BaseModel):
    id: str
    start_pos: int
    end_pos: int
    animal_seq: str
    human_seq: str
    divergence_score: Optional[float]
    esm_distance: Optional[float]


class OrthologOut(BaseModel):
    id: str
    species_id: str
    protein_id: Optional[str]
    sequence_identity_pct: Optional[float]
    orthofinder_og: Optional[str]
    motifs: list[MotifOut] = []


class EvolutionOut(BaseModel):
    dnds_ratio: Optional[float]
    dnds_pvalue: Optional[float]
    selection_model: Optional[str]
    branches_under_selection: Optional[list[str]]
    convergence_count: Optional[int]
    phylop_score: Optional[float]


class CandidateListItem(BaseModel):
    gene_id: str
    gene_symbol: str
    human_protein: Optional[str]
    composite_score: float
    tier: str
    convergence_score: float
    selection_score: float
    expression_score: float
    updated_at: Optional[str]


class CandidateDetail(BaseModel):
    gene_id: str
    gene_symbol: str
    human_protein: Optional[str]
    composite_score: float
    tier: str
    convergence_score: float
    selection_score: float
    expression_score: float
    evolution: Optional[EvolutionOut]
    orthologs: list[OrthologOut] = []
    updated_at: Optional[str]


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------


@router.get("", response_model=list[CandidateListItem])
def list_candidates(
    tier: Optional[str] = Query(None, description="Filter by tier: Tier1, Tier2, Tier3"),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    min_score: Optional[float] = Query(None, description="Minimum composite score"),
):
    """Return ranked candidates sorted by composite score descending.

    Filters: `tier` (Tier1/Tier2/Tier3), `min_score`, pagination with `limit`/`offset`.
    """
    with get_session() as session:
        q = (
            session.query(CandidateScore, Gene)
            .join(Gene, CandidateScore.gene_id == Gene.id)
        )
        if tier:
            q = q.filter(CandidateScore.tier == tier)
        if min_score is not None:
            q = q.filter(CandidateScore.composite_score >= min_score)

        q = q.order_by(CandidateScore.composite_score.desc())
        q = q.offset(offset).limit(limit)

        results = []
        for cs, gene in q:
            results.append(CandidateListItem(
                gene_id=gene.id,
                gene_symbol=gene.gene_symbol,
                human_protein=gene.human_protein,
                composite_score=cs.composite_score or 0.0,
                tier=cs.tier or "Tier3",
                convergence_score=cs.convergence_score or 0.0,
                selection_score=cs.selection_score or 0.0,
                expression_score=cs.expression_score or 0.0,
                updated_at=cs.updated_at.isoformat() if cs.updated_at else None,
            ))
        return results


@router.get("/{gene_id}", response_model=CandidateDetail)
def get_candidate(gene_id: str):
    """Full candidate detail including motifs and orthologs."""
    with get_session() as session:
        gene = session.get(Gene, gene_id)
        if gene is None:
            raise HTTPException(status_code=404, detail=f"Gene '{gene_id}' not found.")

        cs = session.get(CandidateScore, gene_id)
        if cs is None:
            raise HTTPException(status_code=404, detail=f"No score found for gene '{gene_id}'.")

        ev = session.get(EvolutionScore, gene_id)
        evolution_out = None
        if ev:
            evolution_out = EvolutionOut(
                dnds_ratio=ev.dnds_ratio,
                dnds_pvalue=ev.dnds_pvalue,
                selection_model=ev.selection_model,
                branches_under_selection=ev.branches_under_selection,
                convergence_count=ev.convergence_count,
                phylop_score=ev.phylop_score,
            )

        orthologs_out = []
        for orth in session.query(Ortholog).filter_by(gene_id=gene_id).all():
            motifs_out = [
                MotifOut(
                    id=m.id,
                    start_pos=m.start_pos,
                    end_pos=m.end_pos,
                    animal_seq=m.animal_seq,
                    human_seq=m.human_seq,
                    divergence_score=m.divergence_score,
                    esm_distance=m.esm_distance,
                )
                for m in orth.motifs
            ]
            orthologs_out.append(OrthologOut(
                id=orth.id,
                species_id=orth.species_id,
                protein_id=orth.protein_id,
                sequence_identity_pct=orth.sequence_identity_pct,
                orthofinder_og=orth.orthofinder_og,
                motifs=motifs_out,
            ))

        return CandidateDetail(
            gene_id=gene.id,
            gene_symbol=gene.gene_symbol,
            human_protein=gene.human_protein,
            composite_score=cs.composite_score or 0.0,
            tier=cs.tier or "Tier3",
            convergence_score=cs.convergence_score or 0.0,
            selection_score=cs.selection_score or 0.0,
            expression_score=cs.expression_score or 0.0,
            evolution=evolution_out,
            orthologs=orthologs_out,
            updated_at=cs.updated_at.isoformat() if cs.updated_at else None,
        )
