"""AAV compatibility — gene size and packaging limit, tissue tropism."""

import logging
import time
from typing import Optional

import requests

from db.models import Gene, GeneTherapyScore
from db.session import get_session
from pipeline.config import get_ncbi_api_key, get_ncbi_email

log = logging.getLogger(__name__)

AAV_PACKAGING_LIMIT_BP = 4700  # ~4.7 kb payload limit for AAV
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

# Serotype → primary tissue (for recommendation)
AAV_TROPISM = {
    "AAV9": ["CNS", "heart", "liver", "skeletal muscle"],
    "AAV8": ["liver", "pancreas"],
    "AAV2": ["CNS", "retina", "muscle"],
    "AAV5": ["lung", "CNS", "retina"],
    "AAV6": ["heart", "lung", "muscle"],
    "AAV1": ["muscle", "CNS"],
}


def _ncbi_gene_size(gene_symbol: str) -> Optional[int]:
    """Get gene length (bp) from NCBI Gene for human."""
    try:
        r = requests.get(
            NCBI_ESEARCH,
            params={
                "db": "gene",
                "term": f"{gene_symbol}[Gene Name] AND 9606[Taxonomy ID]",
                "retmax": 1,
                "retmode": "json",
                "api_key": get_ncbi_api_key(),
                "email": get_ncbi_email(),
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        id_list = data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return None
        time.sleep(0.12)
        r2 = requests.get(
            NCBI_ESUMMARY,
            params={
                "db": "gene",
                "id": id_list[0],
                "retmode": "json",
                "api_key": get_ncbi_api_key(),
                "email": get_ncbi_email(),
            },
            timeout=15,
        )
        if r2.status_code != 200:
            return None
        summ = r2.json()
        result = summ.get("result", {}).get(id_list[0], {})
        genomicinfo = result.get("genomicinfo", [{}])
        if not genomicinfo:
            return None
        start = genomicinfo[0].get("chrstart")
        end = genomicinfo[0].get("chrstop")
        if start is not None and end is not None:
            return int(end) - int(start)
    except Exception as exc:
        log.debug("NCBI gene size %s: %s", gene_symbol, exc)
    return None


def annotate_genes_aav(gene_ids: list[str]) -> int:
    """Populate GeneTherapyScore.gene_size_bp, aav_compatible, tissue_tropism."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            size = _ncbi_gene_size(gene.gene_symbol)
            gs = session.get(GeneTherapyScore, gid)
            if gs is None:
                gs = GeneTherapyScore(gene_id=gid)
                session.add(gs)
            gs.gene_size_bp = size
            gs.aav_compatible = (size is not None and size <= AAV_PACKAGING_LIMIT_BP)
            gs.tissue_tropism = AAV_TROPISM.get("AAV9", ["CNS", "liver"])  # default suggestion
            updated += 1
    log.info("AAV: updated %d genes.", updated)
    return updated
