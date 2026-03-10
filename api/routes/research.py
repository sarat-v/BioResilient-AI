"""Research assistant: gene search and narrative generation."""

import json
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel
from sqlalchemy import or_

from db.models import CandidateScore, Gene
from db.session import get_session
from pipeline.research_assistant.narrative import generate_narrative

_CONFIG_DIR = Path(__file__).resolve().parents[2] / "config"
_TRAIT_PRESETS_PATH = _CONFIG_DIR / "trait_presets.json"

router = APIRouter()


class NarrativeRequest(BaseModel):
    gene_id: str
    force: bool = False


class NarrativeResponse(BaseModel):
    gene_id: str
    narrative: str
    cached: bool


class SearchHit(BaseModel):
    gene_id: str
    gene_symbol: str
    human_gene_id: Optional[str]
    human_protein: Optional[str]


class PathwayConvergenceOut(BaseModel):
    pathway_id: str
    pathway_name: Optional[str]
    gene_count: Optional[int]
    candidate_count: Optional[int]
    log_pvalue: Optional[float]
    evolutionary_weight: Optional[float]
    pathway_score: Optional[float]
    gene_symbols: Optional[list[str]]


@router.get("/search", response_model=list[SearchHit])
def search_genes(q: str = Query(..., min_length=1, description="Search by symbol, NCBI Gene ID, or UniProt accession")):
    """Full-text search across gene symbol, human_gene_id, and human_protein."""
    term = (q or "").strip().lower()
    if not term:
        return []
    with get_session() as session:
        genes = (
            session.query(Gene)
            .filter(
                Gene.gene_symbol.ilike(f"%{term}%")
                | Gene.human_gene_id.ilike(f"%{term}%")
                | (Gene.human_protein.isnot(None) & Gene.human_protein.ilike(f"%{term}%"))
            )
            .limit(50)
            .all()
        )
        return [
            SearchHit(
                gene_id=g.id,
                gene_symbol=g.gene_symbol,
                human_gene_id=g.human_gene_id,
                human_protein=g.human_protein,
            )
            for g in genes
        ]


@router.get("/traits", response_model=list)
def get_traits():
    """Return trait presets from config/trait_presets.json for the preset selector."""
    if not _TRAIT_PRESETS_PATH.exists():
        return []
    with open(_TRAIT_PRESETS_PATH) as f:
        return json.load(f)


@router.get("/pathways", response_model=dict)
def get_pathways():
    """Return distinct GO terms and pathway IDs from Tier1/Tier2 genes for filter dropdown."""
    with get_session() as session:
        rows = (
            session.query(Gene)
            .join(CandidateScore, (CandidateScore.gene_id == Gene.id) & (CandidateScore.trait_id == ""))
            .filter(CandidateScore.tier.in_(["Tier1", "Tier2"]))
            .all()
        )
        go_set = set()
        pathway_set = set()
        for g in rows:
            if g.go_terms:
                go_set.update(g.go_terms)
            if g.pathway_ids:
                pathway_set.update(g.pathway_ids)
        return {"go_terms": sorted(go_set)[:200], "pathway_ids": sorted(pathway_set)[:200]}


@router.post("/narrative", response_model=NarrativeResponse)
def post_narrative(body: NarrativeRequest):
    """Generate or return cached plain-English research summary for a gene."""
    try:
        narrative, cached = generate_narrative(body.gene_id, force=body.force)
        return NarrativeResponse(gene_id=body.gene_id, narrative=narrative, cached=cached)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except RuntimeError as e:
        raise HTTPException(status_code=502, detail=str(e))


@router.get("/pathway-convergence", response_model=list[PathwayConvergenceOut])
def get_pathway_convergence(
    top: int = Query(50, ge=1, le=500, description="Maximum number of pathways to return"),
    min_candidates: int = Query(2, ge=1, description="Minimum candidate genes per pathway"),
):
    """Return pathway-level convergence enrichment scores.

    Shows which biological pathways are over-represented in convergent
    candidate genes, helping identify the biological *systems* rather than
    just individual genes that are under selection.
    """
    from db.models import PathwayConvergence

    with get_session() as session:
        rows = (
            session.query(PathwayConvergence)
            .filter(PathwayConvergence.candidate_count >= min_candidates)
            .order_by(PathwayConvergence.pathway_score.desc())
            .limit(top)
            .all()
        )
        return [
            PathwayConvergenceOut(
                pathway_id=r.pathway_id,
                pathway_name=r.pathway_name,
                gene_count=r.gene_count,
                candidate_count=r.candidate_count,
                log_pvalue=r.log_pvalue,
                evolutionary_weight=r.evolutionary_weight,
                pathway_score=r.pathway_score,
                gene_symbols=r.gene_symbols,
            )
            for r in rows
        ]
