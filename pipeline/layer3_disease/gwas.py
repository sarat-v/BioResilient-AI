"""GWAS Catalog REST API — strongest association p-value per gene."""

import logging
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene
from db.session import get_session

log = logging.getLogger(__name__)

GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"


def fetch_gwas_pvalue(gene_symbol: str) -> Optional[float]:
    """Fetch the strongest (minimum) p-value from GWAS Catalog for the given gene symbol."""
    try:
        r = requests.get(
            f"{GWAS_API}/genes",
            params={"geneName": gene_symbol, "size": 100},
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        embeds = data.get("_embedded", {}).get("genes", [])
        if not embeds:
            return None
        gene_uri = embeds[0].get("_links", {}).get("self", {}).get("href")
        if not gene_uri:
            return None
        r2 = requests.get(f"https://www.ebi.ac.uk{gene_uri}/associations", params={"size": 200}, timeout=15)
        if r2.status_code != 200:
            return None
        assoc = r2.json()
        assoc_list = assoc.get("_embedded", {}).get("associations", [])
        pvalues = []
        for a in assoc_list:
            pval = a.get("pvalue")
            if pval is not None:
                try:
                    pvalues.append(float(pval))
                except (TypeError, ValueError):
                    pass
        return min(pvalues) if pvalues else None
    except Exception as exc:
        log.debug("GWAS for %s: %s", gene_symbol, exc)
        return None


def annotate_genes_gwas(gene_ids: list[str]) -> int:
    """Populate DiseaseAnnotation.gwas_pvalue for the given genes."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            pval = fetch_gwas_pvalue(gene.gene_symbol)
            if pval is None:
                continue
            ann = session.get(DiseaseAnnotation, gid)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gid)
                session.add(ann)
            ann.gwas_pvalue = pval
            updated += 1
    log.info("GWAS Catalog: updated %d genes.", updated)
    return updated
