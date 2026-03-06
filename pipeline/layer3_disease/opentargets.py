"""OpenTargets Genetics API — disease associations and scores for genes."""

import logging
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene
from db.session import get_session

log = logging.getLogger(__name__)

OPENTARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

QUERY_ASSOCIATIONS = """
query AssociationScore($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    associatedDiseases {
      count
      rows {
        disease {
          id
          name
        }
        score
      }
    }
  }
}
"""


def _symbol_to_ensembl(session, gene_id: str) -> Optional[str]:
    """Resolve gene to Ensembl ID. OpenTargets uses Ensembl; we may have symbol or UniProt."""
    gene = session.get(Gene, gene_id)
    if not gene:
        return None
    # OpenTargets search by symbol - we need Ensembl ID. Use API to search.
    try:
        r = requests.post(
            OPENTARGETS_GRAPHQL,
            json={
                "query": """
                query SearchTarget($queryString: String!) {
                  search(queryString: $queryString, entityNames: ["target"]) {
                    hits {
                      id
                      entity
                      name
                    }
                  }
                }
                """,
                "variables": {"queryString": gene.gene_symbol},
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        hits = data.get("data", {}).get("search", {}).get("hits", [])
        for h in hits:
            if h.get("entity") == "target" and h.get("id", "").startswith("ENSG"):
                return h["id"]
    except Exception as exc:
        log.debug("OpenTargets search for %s: %s", gene.gene_symbol, exc)
    return None


def fetch_opentargets_score(gene_id: str, gene_symbol: str, ensembl_id: Optional[str] = None) -> Optional[float]:
    """Fetch OpenTargets association score for a gene. Returns max disease association score or None."""
    if not ensembl_id:
        with get_session() as session:
            ensembl_id = _symbol_to_ensembl(session, gene_id)
    if not ensembl_id:
        return None
    try:
        r = requests.post(
            OPENTARGETS_GRAPHQL,
            json={
                "query": QUERY_ASSOCIATIONS,
                "variables": {"ensemblId": ensembl_id},
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        target = data.get("data", {}).get("target")
        if not target:
            return None
        rows = target.get("associatedDiseases", {}).get("rows", [])
        if not rows:
            return None
        scores = [r.get("score") for r in rows if r.get("score") is not None]
        return round(max(scores), 4) if scores else None
    except Exception as exc:
        log.debug("OpenTargets for %s: %s", gene_symbol, exc)
        return None


def annotate_genes_opentargets(gene_ids: list[str]) -> int:
    """Populate DiseaseAnnotation.opentargets_score for the given genes."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            score = fetch_opentargets_score(gid, gene.gene_symbol)
            if score is None:
                continue
            ann = session.get(DiseaseAnnotation, gid)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gid)
                session.add(ann)
            ann.opentargets_score = score
            updated += 1
    log.info("OpenTargets: updated %d genes.", updated)
    return updated
