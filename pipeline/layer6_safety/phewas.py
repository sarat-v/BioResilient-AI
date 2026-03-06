"""PheWAS / GWAS Catalog phenome-wide associations — safety signal."""

import logging
from typing import Any, Optional

import requests

from db.models import Gene, SafetyFlag
from db.session import get_session

log = logging.getLogger(__name__)

GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"


def fetch_phewas_hits(gene_symbol: str) -> Optional[dict[str, Any]]:
    """Fetch phenotype–pvalue map from GWAS Catalog for the gene (phenome-wide style). Returns {trait: pvalue}."""
    try:
        r = requests.get(
            f"{GWAS_API}/genes",
            params={"geneName": gene_symbol, "size": 1},
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        embeds = data.get("_embedded", {}).get("genes", [])
        if not embeds:
            return None
        self_link = embeds[0].get("_links", {}).get("self", {}).get("href", "")
        if not self_link.startswith("/"):
            return None
        r2 = requests.get(
            f"https://www.ebi.ac.uk{self_link}/associations",
            params={"size": 500},
            timeout=15,
        )
        if r2.status_code != 200:
            return None
        assoc = r2.json()
        assoc_list = assoc.get("_embedded", {}).get("associations", [])
        hits = {}
        for a in assoc_list:
            efo = a.get("efoTraits") or []
            pval = a.get("pvalue")
            for trait in efo:
                name = trait.get("trait", "") if isinstance(trait, dict) else str(trait)
                if name and pval is not None:
                    try:
                        p = float(pval)
                        if name not in hits or p < hits[name]:
                            hits[name] = p
                    except (TypeError, ValueError):
                        pass
        return hits if hits else None
    except Exception as exc:
        log.debug("PheWAS %s: %s", gene_symbol, exc)
        return None


def annotate_genes_phewas(gene_ids: list[str]) -> int:
    """Populate SafetyFlag.phewas_hits for the given genes."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            hits = fetch_phewas_hits(gene.gene_symbol)
            if hits is None:
                continue
            sf = session.get(SafetyFlag, gid)
            if sf is None:
                sf = SafetyFlag(gene_id=gid)
                session.add(sf)
            sf.phewas_hits = hits
            updated += 1
    log.info("PheWAS: updated %d genes.", updated)
    return updated
