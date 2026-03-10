"""GET /candidates and GET /candidates/{gene_id}."""

from typing import Optional

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import Response
from pydantic import BaseModel
from sqlalchemy import or_

from db.models import (
    CandidateScore,
    DiseaseAnnotation,
    DivergentMotif,
    DrugTarget,
    EvolutionScore,
    Gene,
    Ortholog,
)
from db.session import get_session

import csv
import io

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


class DiseaseOut(BaseModel):
    disease_id: Optional[str]
    disease_name: Optional[str]
    opentargets_score: Optional[float]
    gwas_pvalue: Optional[float]
    gnomad_pli: Optional[float]
    mouse_ko_phenotype: Optional[str]


class DrugTargetOut(BaseModel):
    pocket_count: Optional[int]
    top_pocket_score: Optional[float]
    chembl_target_id: Optional[str]
    existing_drugs: Optional[list[str]]
    cansar_score: Optional[float]
    druggability_tier: Optional[str]


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
    disease: Optional[DiseaseOut] = None
    drug_target: Optional[DrugTargetOut] = None
    narrative: Optional[str] = None
    updated_at: Optional[str]


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------


@router.get("", response_model=list[CandidateListItem])
def list_candidates(
    tier: Optional[str] = Query(None, description="Filter by tier: Tier1, Tier2, Tier3"),
    species_id: Optional[str] = Query(None, description="Filter by species (genes with ortholog in this species)"),
    trait_id: Optional[str] = Query(None, description="Trait preset id (e.g. cancer_resistance). Omit for default."),
    pathway: Optional[str] = Query(None, description="Filter by GO term or pathway ID (gene must have this in go_terms or pathway_ids)"),
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    min_score: Optional[float] = Query(None, description="Minimum composite score"),
):
    """Return ranked candidates sorted by composite score descending.

    Filters: `tier` (Tier1/Tier2/Tier3), `species_id`, `trait_id`, `min_score`, pagination with `limit`/`offset`.
    """
    tid = trait_id if trait_id is not None else ""
    with get_session() as session:
        q = (
            session.query(CandidateScore, Gene)
            .join(Gene, CandidateScore.gene_id == Gene.id)
            .filter(CandidateScore.trait_id == tid)
        )
        if species_id:
            q = q.filter(
                Gene.id.in_(session.query(Ortholog.gene_id).filter(Ortholog.species_id == species_id).distinct())
            )
        if pathway:
            q = q.filter(
                or_(
                    Gene.go_terms.contains([pathway]),
                    Gene.pathway_ids.contains([pathway]),
                )
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


@router.get("/export")
def export_candidates(
    trait_id: Optional[str] = Query(None),
    tier: Optional[str] = Query(None),
    limit: int = Query(5000, ge=1, le=50000),
):
    """Export candidates as CSV with scores and Phase 2 data."""
    tid = trait_id if trait_id is not None else ""
    with get_session() as session:
        q = (
            session.query(CandidateScore, Gene)
            .join(Gene, CandidateScore.gene_id == Gene.id)
            .filter(CandidateScore.trait_id == tid)
        )
        if tier:
            q = q.filter(CandidateScore.tier == tier)
        q = q.order_by(CandidateScore.composite_score.desc()).limit(limit)
        rows = q.all()

        da_map = {}
        dt_map = {}
        for gene_id in [g.id for _, g in rows]:
            da = session.get(DiseaseAnnotation, gene_id)
            if da:
                da_map[gene_id] = da
            dt = session.get(DrugTarget, gene_id)
            if dt:
                dt_map[gene_id] = dt

    out = io.StringIO()
    w = csv.writer(out)
    w.writerow([
        "gene_id", "gene_symbol", "human_protein", "tier", "composite_score",
        "convergence_score", "selection_score", "expression_score", "disease_score", "druggability_score", "safety_score", "regulatory_score",
        "disease_name", "opentargets_score", "pocket_count", "top_pocket_score", "druggability_tier", "existing_drugs",
    ])
    for cs, gene in rows:
        da = da_map.get(gene.id)
        dt = dt_map.get(gene.id)
        w.writerow([
            gene.id,
            gene.gene_symbol,
            gene.human_protein or "",
            cs.tier or "",
            f"{cs.composite_score or 0:.4f}",
            f"{cs.convergence_score or 0:.4f}",
            f"{cs.selection_score or 0:.4f}",
            f"{cs.expression_score or 0:.4f}",
            f"{cs.disease_score or 0:.4f}",
            f"{cs.druggability_score or 0:.4f}",
            f"{cs.safety_score or 0:.4f}",
            f"{cs.regulatory_score or 0:.4f}",
            da.disease_name if da else "",
            f"{da.opentargets_score:.4f}" if da and da.opentargets_score is not None else "",
            dt.pocket_count if dt else "",
            f"{dt.top_pocket_score:.4f}" if dt and dt.top_pocket_score is not None else "",
            dt.druggability_tier if dt else "",
            ";".join(dt.existing_drugs or []) if dt else "",
        ])
    return Response(content=out.getvalue(), media_type="text/csv", headers={"Content-Disposition": "attachment; filename=candidates.csv"})


@router.get("/{gene_id}", response_model=CandidateDetail)
def get_candidate(gene_id: str, trait_id: Optional[str] = Query(None, description="Trait preset id. Omit for default.")):
    """Full candidate detail including motifs and orthologs."""
    tid = trait_id if trait_id is not None else ""
    with get_session() as session:
        gene = session.get(Gene, gene_id)
        if gene is None:
            raise HTTPException(status_code=404, detail=f"Gene '{gene_id}' not found.")

        cs = session.query(CandidateScore).filter_by(gene_id=gene_id, trait_id=tid).first()
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

        da = session.get(DiseaseAnnotation, gene_id)
        disease_out = None
        if da:
            disease_out = DiseaseOut(
                disease_id=da.disease_id,
                disease_name=da.disease_name,
                opentargets_score=da.opentargets_score,
                gwas_pvalue=da.gwas_pvalue,
                gnomad_pli=da.gnomad_pli,
                mouse_ko_phenotype=da.mouse_ko_phenotype,
            )

        dt = session.get(DrugTarget, gene_id)
        drug_target_out = None
        if dt:
            drug_target_out = DrugTargetOut(
                pocket_count=dt.pocket_count,
                top_pocket_score=dt.top_pocket_score,
                chembl_target_id=dt.chembl_target_id,
                existing_drugs=dt.existing_drugs or [],
                cansar_score=dt.cansar_score,
                druggability_tier=dt.druggability_tier,
            )

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
            disease=disease_out,
            drug_target=drug_target_out,
            narrative=getattr(gene, "narrative", None),
            updated_at=cs.updated_at.isoformat() if cs.updated_at else None,
        )
