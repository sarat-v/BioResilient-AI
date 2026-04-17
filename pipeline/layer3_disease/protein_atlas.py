"""Human Protein Atlas API — tissue expression (TPM) per gene."""

import logging
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene
from db.session import get_session

log = logging.getLogger(__name__)

HPA_SEARCH = "https://www.proteinatlas.org/search_download.php"

_SPECIES_SUFFIXES = ("_HUMAN", "_MOUSE", "_RAT", "_BOVIN", "_DANRE", "_YEAST", "_CAEEL", "_DROME")


def _hgnc_symbol(gene_symbol: str) -> str:
    """Strip UniProt species suffix to get an HGNC-style gene symbol (AKT1_HUMAN → AKT1)."""
    upper = gene_symbol.upper()
    for suffix in _SPECIES_SUFFIXES:
        if upper.endswith(suffix):
            return gene_symbol[: -len(suffix)]
    return gene_symbol


def fetch_tissue_expression(gene_symbol: str) -> Optional[dict]:
    """Fetch tissue TPM expression map from Human Protein Atlas. Returns {tissue: tpm_value}."""
    try:
        r = requests.get(
            HPA_SEARCH,
            params={
                "search": gene_symbol,
                "format": "json",
                "columns": "g, Tissue",
            },
            timeout=15,
            headers={"Accept": "application/json"},
        )
        if r.status_code != 200:
            return None
        data = r.json()
        # Response format: list of dicts with Gene, Tissue, TPM etc.
        if not isinstance(data, list) or not data:
            return None
        tissue_tpm = {}
        for row in data:
            if row.get("Gene") and row.get("Gene").upper() == gene_symbol.upper():
                tissue = row.get("Tissue", "unknown")
                tpm = row.get("TPM", row.get("Level", 0))
                try:
                    tissue_tpm[tissue] = float(tpm)
                except (TypeError, ValueError):
                    tissue_tpm[tissue] = 0.0
        return tissue_tpm if tissue_tpm else None
    except Exception as exc:
        log.debug("HPA for %s: %s", gene_symbol, exc)
        return None


def annotate_genes_protein_atlas(gene_ids: list[str]) -> int:
    """Populate DiseaseAnnotation.tissue_expression for the given genes."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            expr = fetch_tissue_expression(_hgnc_symbol(gene.gene_symbol))
            if expr is None:
                continue
            ann = session.get(DiseaseAnnotation, gid)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gid)
                session.add(ann)
            ann.tissue_expression = expr
            updated += 1
    log.info("Human Protein Atlas: updated %d genes.", updated)
    return updated
