"""Step 16 — Translational status lookup.

For each Tier1/Tier2 gene, queries:
  1. ClinicalTrials.gov API v2 — find trials targeting this gene by name
  2. ChEMBL drug_indication — get maximum approved clinical phase for known drugs

Results update DiseaseAnnotation.known_drug_phase and best_protective_trait.
"""

import logging
import time
from typing import Optional

import requests

from db.models import DiseaseAnnotation
from db.session import get_session

log = logging.getLogger(__name__)

CLINICALTRIALS_API = "https://clinicaltrials.gov/api/v2/studies"
CHEMBL_DRUG_INDICATION_API = "https://www.ebi.ac.uk/chembl/api/data/drug_indication.json"
REQUEST_TIMEOUT = 15

_SPECIES_SUFFIXES = ("_HUMAN", "_MOUSE", "_RAT", "_DOG", "_CAT", "_WHALE", "_BAT", "_SHARK")


def _hgnc(symbol: str) -> str:
    for suf in _SPECIES_SUFFIXES:
        if symbol.upper().endswith(suf):
            return symbol[: -len(suf)]
    return symbol


def _fetch_clinicaltrials(gene_symbol: str) -> tuple[int, Optional[str]]:
    """Return (trial_count, top_condition_name) for trials mentioning this gene."""
    try:
        r = requests.get(
            CLINICALTRIALS_API,
            params={
                "query.term": gene_symbol,
                "pageSize": 10,
                "format": "json",
            },
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            return 0, None
        data = r.json()
        total = int(data.get("totalCount", 0))
        studies = data.get("studies", [])
        # Extract first condition name from most relevant study
        condition = None
        for study in studies:
            protocol = study.get("protocolSection", {})
            conditions = (protocol.get("conditionsModule") or {}).get("conditions", [])
            if conditions:
                condition = conditions[0]
                break
        return total, condition
    except Exception as exc:
        log.debug("ClinicalTrials fetch for %s: %s", gene_symbol, exc)
        return 0, None


def _fetch_chembl_max_phase(gene_symbol: str) -> int:
    """Return maximum clinical phase from ChEMBL drug_indication for this target."""
    try:
        r = requests.get(
            CHEMBL_DRUG_INDICATION_API,
            params={
                "target_chembl_id__gene_symbol": gene_symbol,
                "limit": 20,
                "format": "json",
            },
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            return 0
        results = r.json().get("drug_indications", [])
        phases = [
            int(d.get("max_phase_for_ind") or 0)
            for d in results
            if d.get("max_phase_for_ind") is not None
        ]
        return max(phases) if phases else 0
    except Exception as exc:
        log.debug("ChEMBL phase fetch for %s: %s", gene_symbol, exc)
        return 0


def annotate_translational_status(gene_ids: list[str]) -> int:
    """Populate DiseaseAnnotation with ClinicalTrials + ChEMBL translational data."""
    from db.models import Gene

    with get_session() as session:
        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()
        gene_map = {g.id: g for g in genes}

    updated = 0
    for gene_id in gene_ids:
        gene = gene_map.get(gene_id)
        if not gene or not gene.gene_symbol:
            continue

        symbol = _hgnc(gene.gene_symbol)

        trial_count, condition = _fetch_clinicaltrials(symbol)
        time.sleep(0.15)

        max_phase = _fetch_chembl_max_phase(symbol)
        time.sleep(0.15)

        log.info(
            "Translational %s: trials=%d max_phase=%d condition=%s",
            symbol, trial_count, max_phase, condition or "—",
        )

        with get_session() as session:
            ann = session.get(DiseaseAnnotation, gene_id)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gene_id)
                session.add(ann)

            if max_phase > (ann.known_drug_phase or 0):
                ann.known_drug_phase = max_phase

            if condition and not ann.best_protective_trait:
                ann.best_protective_trait = condition

            session.commit()

        updated += 1

    log.info("Translational status: %d genes annotated.", updated)
    return updated


def run_translational_pipeline(gene_ids: list[str]) -> int:
    """Entry point called from orchestrator step16."""
    return annotate_translational_status(gene_ids)
