"""GWAS Catalog REST API — strongest association p-value per gene."""

import logging
import time
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene
from db.session import get_session

log = logging.getLogger(__name__)

GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"

# EBI asks for "reasonable use". 0.25 s between calls = ~4 req/s per IP.
# At 60 genes × 2 calls each that's 120 requests in ~30 s — well within limits.
_RATE_SLEEP = 0.25


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
    """Populate DiseaseAnnotation.gwas_pvalue for the given genes.

    Idempotent: skips genes where gwas_pvalue is already populated.
    Bulk write: all DB writes in one session after all fetches complete.
    """
    with get_session() as session:
        gene_map = {g.id: g for g in session.query(Gene).filter(Gene.id.in_(gene_ids)).all()}
        already_done = {
            r.gene_id
            for r in session.query(DiseaseAnnotation.gene_id)
            .filter(
                DiseaseAnnotation.gene_id.in_(gene_ids),
                DiseaseAnnotation.gwas_pvalue.isnot(None),
            )
            .all()
        }

    todo = [gid for gid in gene_ids if gid not in already_done and gid in gene_map]
    if not todo:
        log.info("GWAS Catalog: all %d genes already annotated, skipping.", len(already_done))
        return 0

    # Fetch all p-values outside DB session (rate-limited to avoid EBI 429)
    fetch_results: dict[str, float] = {}
    for i, gid in enumerate(todo):
        pval = fetch_gwas_pvalue(gene_map[gid].gene_symbol)
        if pval is not None:
            fetch_results[gid] = pval
        if i < len(todo) - 1:
            time.sleep(_RATE_SLEEP)

    if not fetch_results:
        return 0

    # Bulk write in one session
    updated = 0
    with get_session() as session:
        ann_map = {
            r.gene_id: r
            for r in session.query(DiseaseAnnotation)
            .filter(DiseaseAnnotation.gene_id.in_(list(fetch_results)))
            .all()
        }
        for gid, pval in fetch_results.items():
            ann = ann_map.get(gid)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gid)
                session.add(ann)
            ann.gwas_pvalue = pval
            updated += 1

    log.info("GWAS Catalog: updated %d genes.", updated)
    return updated
