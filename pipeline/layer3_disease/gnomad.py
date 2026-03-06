"""gnomAD GraphQL API — pLI (loss-of-function) constraint score per gene."""

import logging
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene
from db.session import get_session

log = logging.getLogger(__name__)

GNOMAD_GRAPHQL = "https://gnomad.broadinstitute.org/api"


def _symbol_to_ensembl_id(gene_symbol: str) -> Optional[str]:
    """Resolve gene symbol to Ensembl gene ID via gnomAD search."""
    try:
        r = requests.post(
            GNOMAD_GRAPHQL,
            json={
                "query": """
                query GeneSearch($query: String!) {
                  gene_search(query: $query) {
                    hits {
                      gene_id
                      symbol
                    }
                  }
                }
                """,
                "variables": {"query": gene_symbol},
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        hits = data.get("data", {}).get("gene_search", {}).get("hits", [])
        for h in hits:
            if h.get("symbol", "").upper() == gene_symbol.upper():
                return h.get("gene_id")
        if hits:
            return hits[0].get("gene_id")
    except Exception as exc:
        log.debug("gnomAD search for %s: %s", gene_symbol, exc)
    return None


def fetch_gnomad_pli(ensembl_gene_id: str) -> Optional[float]:
    """Fetch pLI (probability loss-of-function intolerant) for a gene from gnomAD."""
    try:
        r = requests.post(
            GNOMAD_GRAPHQL,
            json={
                "query": """
                query GeneConstraint($geneId: String!) {
                  gene(gene_id: $geneId) {
                    pli
                  }
                }
                """,
                "variables": {"geneId": ensembl_gene_id},
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        gene = data.get("data", {}).get("gene")
        if gene is None:
            return None
        pli = gene.get("pli")
        return round(float(pli), 4) if pli is not None else None
    except Exception as exc:
        log.debug("gnomAD pli for %s: %s", ensembl_gene_id, exc)
        return None


def annotate_genes_gnomad(gene_ids: list[str]) -> int:
    """Populate DiseaseAnnotation.gnomad_pli for the given genes."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            ensembl_id = _symbol_to_ensembl_id(gene.gene_symbol)
            if not ensembl_id:
                continue
            pli = fetch_gnomad_pli(ensembl_id)
            if pli is None:
                continue
            ann = session.get(DiseaseAnnotation, gid)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gid)
                session.add(ann)
            ann.gnomad_pli = pli
            updated += 1
    log.info("gnomAD: updated %d genes.", updated)
    return updated
