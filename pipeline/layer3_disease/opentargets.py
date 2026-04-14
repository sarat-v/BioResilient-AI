"""OpenTargets Platform GraphQL API — disease associations, tractability, known drugs, safety.

Queries three things per gene:
  1. associatedDiseases — max association score (existing behaviour)
  2. tractability       — SM / AB / PROTAC modality flags
  3. knownDrugs         — best-phase approved or clinical-stage drug
  4. safetyLiabilities  — known adverse event labels

Tractability and known-drug results are stored on DrugTarget (tractability_sm,
tractability_ab, tractability_protac) and DiseaseAnnotation (known_drug_name,
known_drug_phase, ot_safety_liability) respectively.
"""

import logging
import time
from typing import Optional

import requests

# OpenTargets GraphQL: no official rate limit, but 0.1 s gap avoids CloudFlare
# throttling when many genes are looked up from the same AWS IP.
_RATE_SLEEP = 0.1

from db.models import DiseaseAnnotation, DrugTarget, Gene
from db.session import get_session

log = logging.getLogger(__name__)

OPENTARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

_SPECIES_SUFFIXES = ("_HUMAN", "_MOUSE", "_RAT", "_BOVIN", "_DANRE", "_YEAST", "_CAEEL", "_DROME")


def _hgnc_symbol(gene_symbol: str) -> str:
    """Strip UniProt species suffix to get an HGNC-style gene symbol.

    OMA/UniProt entry names like AKT1_HUMAN → AKT1.
    Symbols without a known suffix are returned unchanged.
    """
    upper = gene_symbol.upper()
    for suffix in _SPECIES_SUFFIXES:
        if upper.endswith(suffix):
            return gene_symbol[: -len(suffix)]
    return gene_symbol

# Combined query: associations + tractability + known drugs + safety in one round-trip
QUERY_FULL = """
query TargetAnnotation($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    associatedDiseases(page: {index: 0, size: 50}, orderByScore: "score") {
      rows {
        disease { id name }
        score
      }
    }
    tractability {
      label
      modality
      value
    }
    drugAndClinicalCandidates {
      rows {
        drug {
          name
          maximumClinicalStage
        }
        maxClinicalStage
      }
    }
    safetyLiabilities {
      event
      effects {
        direction
        dosing
      }
    }
  }
}
"""


def _symbol_to_ensembl_by_symbol(gene_symbol: str) -> Optional[str]:
    """Resolve HGNC gene symbol to Ensembl ID via OpenTargets search (no DB session needed)."""
    try:
        r = requests.post(
            OPENTARGETS_GRAPHQL,
            json={
                "query": """
                query SearchTarget($queryString: String!) {
                  search(queryString: $queryString, entityNames: ["target"]) {
                    hits { id entity name }
                  }
                }
                """,
                "variables": {"queryString": gene_symbol},
            },
            timeout=15,
        )
        if r.status_code != 200:
            log.debug("OpenTargets search HTTP %s for %s", r.status_code, gene_symbol)
            return None
        data = r.json()
        hits = (data.get("data") or {}).get("search", {}).get("hits", [])
        for h in hits:
            if h.get("entity") == "target" and h.get("id", "").startswith("ENSG"):
                return h["id"]
    except Exception as exc:
        log.debug("OpenTargets search for %s: %s", gene_symbol, exc)
    return None


def _symbol_to_ensembl(session, gene_id: str) -> Optional[str]:
    """Resolve gene (by DB id) to Ensembl ID — kept for backward compatibility."""
    gene = session.get(Gene, gene_id)
    if not gene:
        return None
    return _symbol_to_ensembl_by_symbol(_hgnc_symbol(gene.gene_symbol))


def _parse_tractability(tractability_list: list) -> dict[str, bool]:
    """Parse OT tractability list into {sm, ab, protac} bool flags.

    OTP v4 uses short modality codes: 'SM', 'AB', 'PR', 'OC'.
    Older OTP versions used full strings like 'Small molecule', 'PROTAC'.
    We handle both forms. 'value' is True when the target has evidence for
    that modality label (any single True entry sets the flag).
    """
    flags = {"sm": False, "ab": False, "protac": False}
    for item in tractability_list or []:
        modality = (item.get("modality") or "").lower()
        value = item.get("value", False)
        if not value:
            continue
        if modality in ("sm", "small_molecule") or "small" in modality:
            flags["sm"] = True
        elif modality in ("ab", "antibody"):
            flags["ab"] = True
        elif modality in ("pr", "protac") or "degrader" in modality or "protac" in modality:
            flags["protac"] = True
    return flags


_STAGE_TO_INT = {
    "APPROVED": 4, "PHASE_4": 4,
    "PHASE_3B": 3, "PHASE_3": 3,
    "PHASE_2B": 2, "PHASE_2": 2,
    "PHASE_1B": 1, "PHASE_1": 1,
}


def _stage_int(stage_str) -> int:
    """Convert OT v4 maximumClinicalStage string (e.g. 'PHASE_3') to int 1-4."""
    if stage_str is None:
        return 0
    if isinstance(stage_str, (int, float)):
        return int(stage_str)
    return _STAGE_TO_INT.get(str(stage_str).upper(), 0)


def _best_known_drug(known_drugs_rows: list) -> tuple[Optional[str], Optional[int], list[str]]:
    """Return (drug_name, max_phase, indication_ids) for the highest-phase known drug.

    indication_ids is a list of EFO/MONDO disease IDs the best-phase drug is approved/
    studied for. Used by disease_score() to determine whether the drug is relevant to
    the current phenotype, avoiding data-volume bias from off-phenotype drugs.

    Handles both legacy int phase and OT v4 string stage (e.g. 'PHASE_3', 'APPROVED').
    """
    best_name: Optional[str] = None
    best_phase: Optional[int] = None
    best_indication_ids: list[str] = []
    for row in known_drugs_rows or []:
        drug = row.get("drug") or {}
        # OT v4: maxClinicalStage on the row, maximumClinicalStage on the drug
        row_phase = _stage_int(row.get("maxClinicalStage") or row.get("phase"))
        drug_phase = _stage_int(drug.get("maximumClinicalStage") or drug.get("maximumClinicalTrialPhase"))
        effective_phase = max(row_phase, drug_phase)
        if best_phase is None or effective_phase > best_phase:
            best_phase = effective_phase if effective_phase > 0 else None
            best_name = drug.get("name")
            # Collect disease IDs from the indication list on this drug row
            indications = row.get("indications") or row.get("disease") or {}
            if isinstance(indications, dict):
                # Single disease object (OT v4 row-level)
                did = indications.get("id")
                best_indication_ids = [did] if did else []
            elif isinstance(indications, list):
                best_indication_ids = [
                    ind.get("disease", {}).get("id") or ind.get("id")
                    for ind in indications
                    if (ind.get("disease", {}).get("id") or ind.get("id"))
                ]
    return best_name, best_phase, best_indication_ids


def _safety_liability_summary(safety_list: list) -> Optional[str]:
    """Collect distinct safety event labels into a comma-separated string."""
    events = sorted({(item.get("event") or "").strip() for item in safety_list or [] if item.get("event")})
    return ", ".join(events) if events else None


def fetch_opentargets_full(gene_id: str, gene_symbol: str, ensembl_id: Optional[str] = None) -> Optional[dict]:
    """Fetch full OT annotation for a gene. Returns dict with all parsed fields, or None on failure."""
    if not ensembl_id:
        with get_session() as session:
            ensembl_id = _symbol_to_ensembl(session, gene_id)
    if not ensembl_id:
        return None
    try:
        r = requests.post(
            OPENTARGETS_GRAPHQL,
            json={"query": QUERY_FULL, "variables": {"ensemblId": ensembl_id}},
            timeout=20,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        target = (data.get("data") or {}).get("target")
        if not target:
            return None

        # Association score + top-scoring disease name/ID
        rows = (target.get("associatedDiseases") or {}).get("rows") or []
        top_disease_id: Optional[str] = None
        top_disease_name: Optional[str] = None
        top_score: Optional[float] = None
        scores = []
        for row in rows:
            s = row.get("score")
            if s is None:
                continue
            scores.append(s)
            if top_score is None or s > top_score:
                top_score = s
                disease = row.get("disease") or {}
                top_disease_id = disease.get("id")
                top_disease_name = disease.get("name")
        assoc_score = round(max(scores), 4) if scores else None

        # Tractability
        tractability = _parse_tractability(target.get("tractability") or [])

        # Known drugs — OT v4 renamed knownDrugs → drugAndClinicalCandidates
        drug_rows = (target.get("drugAndClinicalCandidates") or {}).get("rows") or []
        known_drug_name, known_drug_phase, known_drug_indication_ids = _best_known_drug(drug_rows)

        # Safety liabilities
        safety_text = _safety_liability_summary(target.get("safetyLiabilities") or [])

        return {
            "ensembl_id": ensembl_id,
            "assoc_score": assoc_score,
            "disease_id": top_disease_id,
            "disease_name": top_disease_name,
            "tractability_sm": tractability["sm"],
            "tractability_ab": tractability["ab"],
            "tractability_protac": tractability["protac"],
            "known_drug_name": known_drug_name,
            "known_drug_phase": known_drug_phase,
            "known_drug_indication_ids": known_drug_indication_ids or [],
            "ot_safety_liability": safety_text,
        }
    except Exception as exc:
        log.debug("OpenTargets full query for %s: %s", gene_symbol, exc)
        return None


# Kept for backward compatibility — delegates to the full query
def fetch_opentargets_score(gene_id: str, gene_symbol: str, ensembl_id: Optional[str] = None) -> Optional[float]:
    result = fetch_opentargets_full(gene_id, gene_symbol, ensembl_id)
    return result["assoc_score"] if result else None


def annotate_genes_opentargets(gene_ids: list[str]) -> int:
    """Populate DiseaseAnnotation and DrugTarget with extended OT data for the given genes.

    Idempotent: skips genes already annotated (opentargets_score is not null).
    Batch write: all DB writes happen in a single session at the end (not N sessions).
    """
    # Load genes and existing annotation state in one query
    with get_session() as session:
        gene_map: dict[str, Gene] = {
            g.id: g for g in session.query(Gene).filter(Gene.id.in_(gene_ids)).all()
        }
        # Skip genes where BOTH opentargets_score AND disease_name are already populated.
        # If disease_name is NULL but score exists, we still re-fetch to backfill the name.
        already_done = {
            r.gene_id
            for r in session.query(DiseaseAnnotation.gene_id)
            .filter(
                DiseaseAnnotation.gene_id.in_(gene_ids),
                DiseaseAnnotation.opentargets_score.isnot(None),
                DiseaseAnnotation.disease_name.isnot(None),
            )
            .all()
        }

    todo = [gid for gid in gene_ids if gid not in already_done and gid in gene_map]
    if not todo:
        log.info("OpenTargets: all %d genes already annotated, skipping.", len(already_done))
        return 0

    log.info("OpenTargets: fetching %d genes (%d already done)...", len(todo), len(already_done))

    # Resolve Ensembl IDs outside any DB session to avoid long-lived transactions.
    # Small sleep avoids CloudFlare/CDN rate-limiting from AWS IPs.
    ensembl_map: dict[str, str] = {}
    for i, gid in enumerate(todo):
        gene = gene_map[gid]
        eid = _symbol_to_ensembl_by_symbol(_hgnc_symbol(gene.gene_symbol))
        if eid:
            ensembl_map[gid] = eid
        if i < len(todo) - 1:
            time.sleep(_RATE_SLEEP)

    log.info("OpenTargets: resolved Ensembl IDs for %d / %d genes.", len(ensembl_map), len(todo))

    # Fetch full annotation — one call per gene with small sleep
    results: dict[str, dict] = {}
    resolved = list(ensembl_map.items())
    for i, (gid, eid) in enumerate(resolved):
        gene = gene_map[gid]
        r = fetch_opentargets_full(gid, gene.gene_symbol, ensembl_id=eid)
        if r is not None:
            results[gid] = r
        if i < len(resolved) - 1:
            time.sleep(_RATE_SLEEP)

    if not results:
        log.info("OpenTargets: no results returned for %d genes.", len(todo))
        return 0

    # Bulk write — single session for all genes
    updated = 0
    with get_session() as session:
        ann_map = {
            r.gene_id: r
            for r in session.query(DiseaseAnnotation)
            .filter(DiseaseAnnotation.gene_id.in_(list(results)))
            .all()
        }
        dt_map = {
            r.gene_id: r
            for r in session.query(DrugTarget)
            .filter(DrugTarget.gene_id.in_(list(results)))
            .all()
        }

        for gid, result in results.items():
            gene = gene_map[gid]

            ann = ann_map.get(gid)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gid)
                session.add(ann)
            if result["assoc_score"] is not None:
                ann.opentargets_score = result["assoc_score"]
            if result["disease_id"] is not None:
                ann.disease_id = result["disease_id"]
            if result["disease_name"] is not None:
                ann.disease_name = result["disease_name"]
            if result["known_drug_name"] is not None:
                ann.known_drug_name = result["known_drug_name"]
            if result["known_drug_phase"] is not None:
                ann.known_drug_phase = result["known_drug_phase"]
            if result.get("known_drug_indication_ids"):
                ann.known_drug_indication_ids = result["known_drug_indication_ids"]
            if result["ot_safety_liability"] is not None:
                ann.ot_safety_liability = result["ot_safety_liability"]

            dt = dt_map.get(gid)
            if dt is None:
                dt = DrugTarget(gene_id=gid)
                session.add(dt)
            dt.tractability_sm = result["tractability_sm"]
            dt.tractability_ab = result["tractability_ab"]
            dt.tractability_protac = result["tractability_protac"]

            updated += 1
            log.debug(
                "OT %s: score=%.3f sm=%s ab=%s protac=%s drug=%s(phase%s) safety=%s",
                gene.gene_symbol,
                result["assoc_score"] or 0,
                result["tractability_sm"],
                result["tractability_ab"],
                result["tractability_protac"],
                result["known_drug_name"],
                result["known_drug_phase"],
                result["ot_safety_liability"],
            )

    log.info("OpenTargets: updated %d genes (associations + tractability + drugs + safety).", updated)
    return updated
