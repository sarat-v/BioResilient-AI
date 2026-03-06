"""ChEMBL REST API — target and approved drugs by UniProt."""

import logging
from typing import Optional

import requests

from db.models import DrugTarget, Gene
from db.session import get_session

log = logging.getLogger(__name__)

CHEMBL_TARGET_BY_UNIPROT = "https://www.ebi.ac.uk/chembl/api/data/target/search"
CHEMBL_DRUGS_FOR_TARGET = "https://www.ebi.ac.uk/chembl/api/data/drug_indication.json"


def fetch_chembl_target_and_drugs(uniprot_id: str) -> tuple[Optional[str], list[str]]:
    """Return (chembl_target_id, list of approved drug names) for a UniProt ID."""
    try:
        r = requests.get(
            CHEMBL_TARGET_BY_UNIPROT,
            params={"q": uniprot_id, "format": "json"},
            timeout=15,
        )
        if r.status_code != 200:
            return None, []
        data = r.json()
        targets = data.get("targets", [])
        if not targets:
            return None, []
        target = targets[0]
        chembl_id = target.get("target_chembl_id")
        if not chembl_id:
            return None, []

        # Fetch molecules (approved drugs) for this target
        r2 = requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
            params={"target_chembl_id": chembl_id, "max_approved_phase": 4},
            timeout=15,
        )
        drugs = []
        if r2.status_code == 200:
            mol_data = r2.json()
            for m in mol_data.get("molecules", [])[:50]:
                if m.get("first_approval"):
                    drugs.append(m.get("pref_name", m.get("molecule_chembl_id", "")))
        return chembl_id, list(dict.fromkeys(drugs))[:20]
    except Exception as exc:
        log.debug("ChEMBL %s: %s", uniprot_id, exc)
        return None, []


def annotate_genes_chembl(gene_ids: list[str]) -> int:
    """Populate DrugTarget.chembl_target_id and existing_drugs for the given genes."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene or not gene.human_protein:
                continue
            chembl_id, drugs = fetch_chembl_target_and_drugs(gene.human_protein.strip())
            if chembl_id is None:
                continue
            dt = session.get(DrugTarget, gid)
            if dt is None:
                dt = DrugTarget(gene_id=gid)
                session.add(dt)
            dt.chembl_target_id = chembl_id
            dt.existing_drugs = drugs if drugs else None
            updated += 1
    log.info("ChEMBL: updated %d genes.", updated)
    return updated
