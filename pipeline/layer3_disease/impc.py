"""IMPC (International Mouse Phenotyping Consortium) — mouse KO phenotype per gene."""

import logging
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene
from db.session import get_session

log = logging.getLogger(__name__)

IMPC_SOLR = "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select"


def fetch_impc_phenotype(gene_symbol: str) -> Optional[str]:
    """Fetch top mouse KO phenotype string from IMPC for the given gene symbol."""
    try:
        r = requests.get(
            IMPC_SOLR,
            params={
                "q": f"marker_symbol:{gene_symbol}",
                "wt": "json",
                "rows": 1,
                "fl": "top_level_mp_term_name,mp_term_name",
                "sort": "phenotype_count desc",
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        response = data.get("response", {})
        docs = response.get("docs", [])
        if not docs:
            return None
        d = docs[0]
        top = d.get("top_level_mp_term_name")
        mp = d.get("mp_term_name")
        if isinstance(top, list):
            top = top[0] if top else ""
        if isinstance(mp, list):
            mp = mp[0] if mp else ""
        return f"{top}; {mp}".strip("; ") if top or mp else None
    except Exception as exc:
        log.debug("IMPC for %s: %s", gene_symbol, exc)
        return None


def annotate_genes_impc(gene_ids: list[str]) -> int:
    """Populate DiseaseAnnotation.mouse_ko_phenotype for the given genes."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            pheno = fetch_impc_phenotype(gene.gene_symbol)
            if pheno is None:
                continue
            ann = session.get(DiseaseAnnotation, gid)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gid)
                session.add(ann)
            ann.mouse_ko_phenotype = pheno[:2000] if len(pheno) > 2000 else pheno
            updated += 1
    log.info("IMPC: updated %d genes.", updated)
    return updated
