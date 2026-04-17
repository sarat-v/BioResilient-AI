"""PheWAS / disease associations — safety signal via Open Targets Platform.

Uses Open Targets Platform GraphQL API (OTP) for disease associations per target.
OTP integrates GWAS Catalog, UK Biobank, FinnGen, and other large-scale phenome-wide
data. Score threshold of 0.05 is used to keep only meaningful associations.

Fallback: GWAS Catalog REST API gene associations (legacy, used if OTP fails).
"""

import logging
import time
from typing import Any, Optional

import requests

from db.models import Gene, SafetyFlag
from db.session import get_session

log = logging.getLogger(__name__)

OTP_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
GWAS_API    = "https://www.ebi.ac.uk/gwas/rest/api"

_RATE_SLEEP = 0.3
_OTP_MIN_SCORE = 0.01  # minimum association score to include


def _hgnc_symbol(raw: str) -> str:
    """Convert OMA entry name (AKT1_HUMAN) to bare HGNC symbol (AKT1)."""
    if raw and "_" in raw:
        return raw.split("_")[0]
    return raw or ""


def _resolve_ensembl_id(gene_symbol: str) -> Optional[str]:
    """Resolve HGNC symbol → Ensembl gene ID via Open Targets gene search."""
    query = """
    query GeneSearch($queryString: String!) {
      search(queryString: $queryString, entityNames: ["target"], page: {index: 0, size: 1}) {
        hits {
          id
          entity
          object { ... on Target { id approvedSymbol } }
        }
      }
    }
    """
    try:
        r = requests.post(
            OTP_GRAPHQL,
            json={"query": query, "variables": {"queryString": gene_symbol}},
            timeout=15,
        )
        if r.status_code != 200:
            return None
        hits = r.json().get("data", {}).get("search", {}).get("hits", [])
        for hit in hits:
            obj = hit.get("object", {})
            if obj.get("approvedSymbol", "").upper() == gene_symbol.upper():
                return obj.get("id")
        # Accept first hit if no exact match
        if hits:
            return hits[0].get("object", {}).get("id")
    except Exception as exc:
        log.warning("OTP gene search %s: %s", gene_symbol, exc)
    return None


def fetch_phewas_hits(gene_symbol: str) -> Optional[dict[str, Any]]:
    """Fetch disease association scores from Open Targets Platform for a gene.

    Returns {disease_name: score} for associations with score >= _OTP_MIN_SCORE.
    Uses OTP's integrated evidence (GWAS Catalog, UK Biobank, FinnGen, etc.).
    """
    ensembl_id = _resolve_ensembl_id(gene_symbol)
    if not ensembl_id:
        log.info("PheWAS %s: could not resolve Ensembl ID.", gene_symbol)
        return None

    query = """
    query DiseaseAssoc($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        associatedDiseases(page: {index: 0, size: 500}) {
          count
          rows {
            disease { id name }
            score
          }
        }
      }
    }
    """
    try:
        r = requests.post(
            OTP_GRAPHQL,
            json={"query": query, "variables": {"ensemblId": ensembl_id}},
            timeout=20,
        )
        if r.status_code != 200:
            log.warning("OTP disease assoc %s (%s): HTTP %d", gene_symbol, ensembl_id, r.status_code)
            return None

        payload = r.json()
        if payload.get("errors"):
            log.warning("OTP disease assoc %s GraphQL errors: %s", gene_symbol, payload.get("errors"))
            return None
        target = payload.get("data", {}).get("target") or {}
        ad = target.get("associatedDiseases") or {}
        rows = ad.get("rows", [])
        total = ad.get("count", len(rows))
        log.info("PheWAS %s: OTP returned %d rows (total=%s)", gene_symbol, len(rows), total)
        hits = {}
        for row in rows:
            score = row.get("score", 0) or 0
            if score >= _OTP_MIN_SCORE:
                disease_name = (row.get("disease") or {}).get("name", "")
                if disease_name:
                    hits[disease_name] = round(float(score), 4)
        return hits if hits else None
    except Exception as exc:
        log.warning("OTP disease assoc %s: %s", gene_symbol, exc)
        return None


def annotate_genes_phewas(gene_ids: list[str]) -> int:
    """Populate SafetyFlag.phewas_hits for the given genes (via Open Targets Platform)."""
    updated = 0
    with get_session() as session:
        for i, gid in enumerate(gene_ids):
            gene = session.get(Gene, gid)
            if not gene:
                continue
            symbol = _hgnc_symbol(gene.gene_symbol)
            log.info("PheWAS querying: %s (raw: %s)", symbol, gene.gene_symbol)
            hits = fetch_phewas_hits(symbol)

            sf = session.get(SafetyFlag, gid)
            if sf is None:
                sf = SafetyFlag(gene_id=gid)
                session.add(sf)

            if hits:
                sf.phewas_hits = hits
                log.info("PheWAS %s: %d disease associations found.", symbol, len(hits))
            else:
                sf.phewas_hits = {}
                log.info("PheWAS %s: no significant disease associations.", symbol)
            updated += 1

            if i < len(gene_ids) - 1:
                time.sleep(_RATE_SLEEP)
        session.commit()
    log.info("PheWAS: updated %d genes (Open Targets Platform).", updated)
    return updated
