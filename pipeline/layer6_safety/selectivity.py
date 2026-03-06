"""Selectivity — protein family size and essentiality (gnomAD pLI) for safety."""

import logging
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene, SafetyFlag
from db.session import get_session

log = logging.getLogger(__name__)

UNIPROT_FAMILY = "https://rest.uniprot.org/uniprotkb/search"
ESSENTIAL_PLI_THRESHOLD = 0.9


def fetch_family_size(uniprot_id: str) -> Optional[int]:
    """Fetch number of proteins in same family (similarity) from UniProt. Returns count or None."""
    try:
        r = requests.get(
            UNIPROT_FAMILY,
            params={
                "query": f"(family:{uniprot_id}) OR (reviewed:true AND identity:0.5)",
                "size": 500,
                "format": "json",
                "fields": "accession",
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        results = data.get("results", [])
        return len(results) if results else 0
    except Exception as exc:
        log.debug("UniProt family %s: %s", uniprot_id, exc)
        return None


def annotate_genes_selectivity(gene_ids: list[str]) -> int:
    """Populate SafetyFlag.family_size and is_essential (from DiseaseAnnotation.gnomad_pli)."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene or not gene.human_protein:
                continue
            family_size = fetch_family_size(gene.human_protein.strip())
            ann = session.get(DiseaseAnnotation, gid)
            pli = ann.gnomad_pli if ann else None
            is_essential = pli is not None and pli >= ESSENTIAL_PLI_THRESHOLD

            sf = session.get(SafetyFlag, gid)
            if sf is None:
                sf = SafetyFlag(gene_id=gid)
                session.add(sf)
            if family_size is not None:
                sf.family_size = family_size
            sf.is_essential = is_essential
            updated += 1
    log.info("Selectivity: updated %d genes.", updated)
    return updated
