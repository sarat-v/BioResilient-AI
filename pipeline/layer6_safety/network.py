"""STRING DB — network degree and hub risk."""

import logging
from typing import Optional

import requests

from db.models import Gene, SafetyFlag
from db.session import get_session

log = logging.getLogger(__name__)

STRING_API = "https://string-db.org/api"
SPECIES_HUMAN = "9606"
HUB_DEGREE_THRESHOLD = 50


def _hgnc_symbol(raw: str) -> str:
    """Convert OMA entry name (AKT1_HUMAN) to bare HGNC symbol (AKT1)."""
    if raw and "_" in raw:
        return raw.split("_")[0]
    return raw or ""


def fetch_string_degree(gene_symbol: str, species: str = SPECIES_HUMAN) -> Optional[int]:
    """Fetch interaction count (degree) for a gene from STRING. Returns None on failure."""
    try:
        r = requests.get(
            f"{STRING_API}/json/network",
            params={
                "identifiers": gene_symbol,
                "species": species,
                "required_score": 700,  # high-confidence interactions only (docs specify >700)
            },
            timeout=15,
        )
        if r.status_code != 200:
            log.warning("STRING %s: HTTP %d", gene_symbol, r.status_code)
            return None
        data = r.json()
        if isinstance(data, list):
            return len(data)
        if isinstance(data, dict):
            return data.get("number_of_edges", data.get("count", 0))
        return 0
    except Exception as exc:
        log.warning("STRING %s: %s", gene_symbol, exc)
        return None


def annotate_genes_network(gene_ids: list[str]) -> int:
    """Populate SafetyFlag.network_degree and hub_risk."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            symbol = _hgnc_symbol(gene.gene_symbol)
            log.info("STRING querying: %s (raw: %s)", symbol, gene.gene_symbol)
            degree = fetch_string_degree(symbol)
            sf = session.get(SafetyFlag, gid)
            if sf is None:
                sf = SafetyFlag(gene_id=gid)
                session.add(sf)
            sf.network_degree = degree if degree is not None else 0
            sf.hub_risk = (degree is not None and degree > HUB_DEGREE_THRESHOLD)
            log.info("STRING %s: degree=%s  hub_risk=%s", symbol, degree, sf.hub_risk)
            updated += 1
        session.commit()
    log.info("STRING network: updated %d genes.", updated)
    return updated
