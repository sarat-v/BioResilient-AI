"""GET /species — species registry."""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Optional

from db.models import Species
from db.session import get_session

router = APIRouter()


class SpeciesOut(BaseModel):
    id: str
    taxid: int
    scientific_name: str
    phenotypes: list[str]
    genome_assembly: Optional[str]
    lineage_group: Optional[str]
    proteome_path: Optional[str]


@router.get("", response_model=list[SpeciesOut])
def list_species():
    """Return all registered species."""
    with get_session() as session:
        species_list = session.query(Species).order_by(Species.id).all()
        return [
            SpeciesOut(
                id=s.id,
                taxid=s.taxid,
                scientific_name=s.scientific_name,
                phenotypes=s.phenotypes or [],
                genome_assembly=s.genome_assembly,
                lineage_group=s.lineage_group,
                proteome_path=s.proteome_path,
            )
            for s in species_list
        ]


@router.get("/{species_id}", response_model=SpeciesOut)
def get_species(species_id: str):
    """Return one species by ID."""
    with get_session() as session:
        species = session.get(Species, species_id)
        if species is None:
            raise HTTPException(status_code=404, detail=f"Species '{species_id}' not found.")
        return SpeciesOut(
            id=species.id,
            taxid=species.taxid,
            scientific_name=species.scientific_name,
            phenotypes=species.phenotypes or [],
            genome_assembly=species.genome_assembly,
            lineage_group=species.lineage_group,
            proteome_path=species.proteome_path,
        )
