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
from typing import Optional

import requests

from db.models import DiseaseAnnotation, DrugTarget, Gene
from db.session import get_session

log = logging.getLogger(__name__)

OPENTARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

# Combined query: associations + tractability + known drugs + safety in one round-trip
QUERY_FULL = """
query TargetAnnotation($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    associatedDiseases(page: {index: 0, size: 20}) {
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
    knownDrugs(size: 5) {
      count
      rows {
        drug {
          name
          maximumClinicalTrialPhase
        }
        disease { name }
        phase
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


def _symbol_to_ensembl(session, gene_id: str) -> Optional[str]:
    """Resolve gene to Ensembl ID via OpenTargets search."""
    gene = session.get(Gene, gene_id)
    if not gene:
        return None
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


def _parse_tractability(tractability_list: list) -> dict[str, bool]:
    """Parse OT tractability list into {sm, ab, protac} bool flags.

    OT modality labels: 'Small molecule', 'Antibody', 'PROTAC'.
    'value' is True if the target has evidence for that modality.
    """
    flags = {"sm": False, "ab": False, "protac": False}
    for item in tractability_list or []:
        modality = (item.get("modality") or "").lower()
        value = item.get("value", False)
        if not value:
            continue
        if "small" in modality or modality == "sm":
            flags["sm"] = True
        elif "antibody" in modality or modality in ("ab", "antibody"):
            flags["ab"] = True
        elif "protac" in modality or "degrader" in modality:
            flags["protac"] = True
    return flags


def _best_known_drug(known_drugs_rows: list) -> tuple[Optional[str], Optional[int]]:
    """Return (drug_name, max_phase) for the highest-phase known drug."""
    best_name: Optional[str] = None
    best_phase: Optional[int] = None
    for row in known_drugs_rows or []:
        drug = row.get("drug") or {}
        phase = row.get("phase")
        drug_max = drug.get("maximumClinicalTrialPhase")
        # Use the higher of row phase and drug-level max phase
        effective_phase = max(
            (int(phase) if phase is not None else 0),
            (int(drug_max) if drug_max is not None else 0),
        )
        if best_phase is None or effective_phase > best_phase:
            best_phase = effective_phase
            best_name = drug.get("name")
    return best_name, best_phase


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

        # Association score
        rows = (target.get("associatedDiseases") or {}).get("rows") or []
        scores = [row.get("score") for row in rows if row.get("score") is not None]
        assoc_score = round(max(scores), 4) if scores else None

        # Tractability
        tractability = _parse_tractability(target.get("tractability") or [])

        # Known drugs
        drug_rows = (target.get("knownDrugs") or {}).get("rows") or []
        known_drug_name, known_drug_phase = _best_known_drug(drug_rows)

        # Safety liabilities
        safety_text = _safety_liability_summary(target.get("safetyLiabilities") or [])

        return {
            "ensembl_id": ensembl_id,
            "assoc_score": assoc_score,
            "tractability_sm": tractability["sm"],
            "tractability_ab": tractability["ab"],
            "tractability_protac": tractability["protac"],
            "known_drug_name": known_drug_name,
            "known_drug_phase": known_drug_phase,
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
    """Populate DiseaseAnnotation and DrugTarget with extended OT data for the given genes."""
    updated = 0
    with get_session() as session:
        # Pre-resolve Ensembl IDs in one pass to avoid repeated symbol lookups
        ensembl_map: dict[str, str] = {}
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            eid = _symbol_to_ensembl(session, gid)
            if eid:
                ensembl_map[gid] = eid

    for gid in gene_ids:
        with get_session() as session:
            gene = session.get(Gene, gid)
            if not gene:
                continue

        result = fetch_opentargets_full(gid, gene.gene_symbol, ensembl_id=ensembl_map.get(gid))
        if result is None:
            continue

        with get_session() as session:
            # DiseaseAnnotation — association score + drug + safety
            ann = session.get(DiseaseAnnotation, gid)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gid)
                session.add(ann)
            if result["assoc_score"] is not None:
                ann.opentargets_score = result["assoc_score"]
            if result["known_drug_name"] is not None:
                ann.known_drug_name = result["known_drug_name"]
            if result["known_drug_phase"] is not None:
                ann.known_drug_phase = result["known_drug_phase"]
            if result["ot_safety_liability"] is not None:
                ann.ot_safety_liability = result["ot_safety_liability"]

            # DrugTarget — tractability flags
            dt = session.get(DrugTarget, gid)
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
