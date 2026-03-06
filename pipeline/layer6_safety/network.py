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


def fetch_string_degree(gene_symbol: str, species: str = SPECIES_HUMAN) -> Optional[int]:
    """Fetch interaction count (degree) for a gene from STRING. Returns None on failure."""
    try:
        r = requests.get(
            f"{STRING_API}/json/network",
            params={
                "identifiers": gene_symbol,
                "species": species,
                "required_score": 400,
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        if isinstance(data, list):
            return len(data)
        if isinstance(data, dict):
            return data.get("number_of_edges", data.get("count", 0))
        return 0
    except Exception as exc:
        log.debug("STRING %s: %s", gene_symbol, exc)
        return None


def annotate_genes_network(gene_ids: list[str]) -> int:
    """Populate SafetyFlag.network_degree and hub_risk."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            degree = fetch_string_degree(gene.gene_symbol)
            if degree is None:
                continue
            sf = session.get(SafetyFlag, gid)
            if sf is None:
                sf = SafetyFlag(gene_id=gid)
                session.add(sf)
            sf.network_degree = degree
            sf.hub_risk = degree > HUB_DEGREE_THRESHOLD
            updated += 1
    log.info("STRING network: updated %d genes.", updated)
    return updated
