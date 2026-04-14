"""Open Targets tractability — replaces CanSAR (no public API).

Open Targets Platform is free, requires no registration, and provides richer
small-molecule tractability data than CanSAR.  Results are stored in the same
DB columns (cansar_score, druggability_tier) so no schema change is needed.

Tier mapping:
  A  — Clinical Precedence (approved drug exists for this target)
  B  — Discovery Precedence (active drug discovery programme)
  C  — Predicted Tractable (ML prediction, high or medium confidence)
  D  — low-confidence predicted or biotherapeutic only
  None / undruggable — no tractability evidence

Score mapping:  A→1.0  B→0.75  C→0.5  D→0.25  None→0.0
"""

import json
import logging
import time
from pathlib import Path
from typing import Optional

import requests

from db.models import DrugTarget, Gene
from db.session import get_session
from pipeline.config import get_local_storage_root

log = logging.getLogger(__name__)

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

# ---------------------------------------------------------------------------
# Shared symbol → Ensembl ID disk cache
# ---------------------------------------------------------------------------
# opentargets.py (step 11) and cansar.py (step 12) both resolve gene symbols to
# Ensembl IDs via the Open Targets GraphQL search endpoint.  This cache avoids
# duplicate HTTP calls: step 11 populates it, step 12 reads it.

_ENSEMBL_CACHE: dict[str, Optional[str]] = {}
_ENSEMBL_CACHE_DIRTY = False


def _ensembl_cache_path() -> Path:
    return Path(get_local_storage_root()) / "ensembl_id_cache.json"


def _load_ensembl_cache() -> None:
    p = _ensembl_cache_path()
    if p.exists():
        try:
            _ENSEMBL_CACHE.update(json.loads(p.read_text()))
        except Exception:
            pass


def _save_ensembl_cache() -> None:
    if not _ENSEMBL_CACHE_DIRTY:
        return
    p = _ensembl_cache_path()
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(_ENSEMBL_CACHE, indent=2))
    except Exception:
        pass

_TIER_SCORE = {"A": 1.0, "B": 0.75, "C": 0.5, "D": 0.25}

# Small-molecule tractability labels from Open Targets v4, in priority order.
# Modality="SM" labels and their tier mapping:
_SM_TIERS: list[tuple[str, str]] = [
    ("Approved Drug", "A"),           # approved drug targets this protein
    ("Advanced Clinical", "A"),       # phase 3 / NDA / BLA compound
    ("Phase 1 Clinical", "B"),        # phase 1 clinical compound
    ("Structure with Ligand", "B"),   # co-crystal structure with drug-like ligand
    ("High-Quality Ligand", "C"),     # high-quality hit compound
    ("High-Quality Pocket", "C"),     # well-defined druggable pocket
    ("Druggable Family", "C"),        # member of a known druggable family
    ("Med-Quality Pocket", "D"),      # medium-confidence pocket
    ("Literature", "D"),              # literature-based evidence only
]


def _hgnc_symbol(gene_symbol: str) -> str:
    """Strip UniProt species suffix (e.g. AKT1_HUMAN → AKT1)."""
    for suffix in ("_HUMAN", "_MOUSE", "_RAT", "_DOG", "_CAT"):
        if gene_symbol.upper().endswith(suffix):
            return gene_symbol[: -len(suffix)]
    return gene_symbol


def _resolve_ensembl_id(hgnc_symbol: str) -> Optional[str]:
    """Search Open Targets for a gene symbol and return its Ensembl ID.

    Results are cached to disk so step 11 (opentargets.py) and step 12 (cansar.py)
    share the same symbol→Ensembl mapping without duplicate HTTP requests.
    """
    if not _ENSEMBL_CACHE:
        _load_ensembl_cache()

    key = hgnc_symbol.upper()
    if key in _ENSEMBL_CACHE:
        return _ENSEMBL_CACHE[key]

    query = """
    query SearchTarget($term: String!) {
      search(queryString: $term, entityNames: ["target"], page: {index: 0, size: 1}) {
        hits { id entity name }
      }
    }"""
    result: Optional[str] = None
    try:
        r = requests.post(
            OT_GRAPHQL,
            json={"query": query, "variables": {"term": hgnc_symbol}},
            timeout=15,
        )
        if r.status_code == 200:
            hits = r.json().get("data", {}).get("search", {}).get("hits", [])
            for hit in hits:
                if hit.get("entity") == "target" and hit.get("name", "").upper() == key:
                    result = hit["id"]
                    break
            if result is None:
                for hit in hits:
                    if hit.get("entity") == "target":
                        result = hit["id"]
                        break
    except Exception as exc:
        log.debug("OT symbol lookup %s: %s", hgnc_symbol, exc)

    global _ENSEMBL_CACHE_DIRTY
    _ENSEMBL_CACHE[key] = result
    _ENSEMBL_CACHE_DIRTY = True
    _save_ensembl_cache()
    return result


def fetch_opentargets_tractability(ensembl_id: str) -> tuple[Optional[float], Optional[str]]:
    """Fetch small-molecule tractability from Open Targets. Returns (score, tier)."""
    query = """
    query Tractability($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        tractability {
          label
          modality
          value
        }
      }
    }"""
    try:
        r = requests.post(
            OT_GRAPHQL,
            json={"query": query, "variables": {"ensemblId": ensembl_id}},
            timeout=15,
        )
        if r.status_code != 200:
            return None, None
        target_data = r.json().get("data", {}).get("target")
        if not target_data:
            return None, None
        tractability = target_data.get("tractability") or []
        # Find best small-molecule (SM) tractability tier
        sm_labels: set[str] = {
            t["label"]
            for t in tractability
            if t.get("modality") == "SM" and t.get("value") is True
        }
        for label, tier in _SM_TIERS:
            if label in sm_labels:
                return _TIER_SCORE[tier], tier
        # Check if any AB (antibody) tier exists as fallback
        ab_labels: set[str] = {
            t["label"]
            for t in tractability
            if t.get("modality") == "AB" and t.get("value") is True
        }
        if "Clinical Precedence" in ab_labels:
            return _TIER_SCORE["B"], "B"
        return None, None
    except Exception as exc:
        log.debug("OT tractability %s: %s", ensembl_id, exc)
    return None, None


def annotate_genes_cansar(gene_ids: list[str]) -> int:
    """Populate DrugTarget.cansar_score and druggability_tier using Open Targets.

    Kept as 'annotate_genes_cansar' so orchestrator.py needs no change.
    """
    updated = 0
    with get_session() as session:
        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()

    for gene in genes:
        symbol = _hgnc_symbol(gene.gene_symbol or "")
        if not symbol:
            continue

        ensembl_id = _resolve_ensembl_id(symbol)
        if not ensembl_id:
            log.debug("OT: no Ensembl ID for %s", symbol)
            continue

        time.sleep(0.1)  # gentle rate limit (~10 req/s)

        score, tier = fetch_opentargets_tractability(ensembl_id)
        if score is None:
            continue

        with get_session() as session:
            dt = session.get(DrugTarget, gene.id)
            if dt is None:
                dt = DrugTarget(gene_id=gene.id)
                session.add(dt)
            dt.cansar_score = score
            if tier:
                dt.druggability_tier = tier[:32]
            session.commit()
        updated += 1
        log.debug("OT tractability %s (%s): tier=%s score=%s", symbol, ensembl_id, tier, score)

    log.info("Open Targets tractability: updated %d genes.", updated)
    return updated
