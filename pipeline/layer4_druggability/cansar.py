"""CanSAR API — druggability score and tier (requires registration)."""

import logging
import os
from typing import Optional

import requests

from db.models import DrugTarget, Gene
from db.session import get_session

log = logging.getLogger(__name__)

CANSAR_API = "https://cansar.ai/api"


def _get_cansar_key() -> str:
    return os.environ.get("CANSAR_API_KEY", "")


def fetch_cansar_score(uniprot_id: str) -> tuple[Optional[float], Optional[str]]:
    """Fetch CanSAR druggability score and tier for a UniProt ID. Returns (score, tier)."""
    key = _get_cansar_key()
    if not key:
        log.debug("CanSAR: no CANSAR_API_KEY set.")
        return None, None
    try:
        r = requests.get(
            f"{CANSAR_API}/target",
            params={"uniprot_id": uniprot_id},
            headers={"Authorization": f"Bearer {key}"},
            timeout=15,
        )
        if r.status_code != 200:
            return None, None
        data = r.json()
        score = data.get("druggability_score")
        tier = data.get("druggability_tier") or data.get("tier")
        if score is not None:
            return round(float(score), 4), (str(tier) if tier else None)
        return None, None
    except Exception as exc:
        log.debug("CanSAR %s: %s", uniprot_id, exc)
        return None, None


def annotate_genes_cansar(gene_ids: list[str]) -> int:
    """Populate DrugTarget.cansar_score and druggability_tier for the given genes."""
    if not _get_cansar_key():
        return 0
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene or not gene.human_protein:
                continue
            score, tier = fetch_cansar_score(gene.human_protein.strip())
            if score is None:
                continue
            dt = session.get(DrugTarget, gid)
            if dt is None:
                dt = DrugTarget(gene_id=gid)
                session.add(dt)
            dt.cansar_score = score
            if tier:
                dt.druggability_tier = tier[:32]
            updated += 1
    log.info("CanSAR: updated %d genes.", updated)
    return updated
