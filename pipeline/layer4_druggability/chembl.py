"""ChEMBL REST API — target and approved drugs by UniProt."""

import logging
from typing import Optional

import requests

from db.models import DrugTarget, Gene
from db.session import get_session

log = logging.getLogger(__name__)

CHEMBL_TARGET_BY_UNIPROT = "https://www.ebi.ac.uk/chembl/api/data/target/search"
CHEMBL_DRUGS_FOR_TARGET = "https://www.ebi.ac.uk/chembl/api/data/drug_indication.json"


def _resolve_uniprot_accession(human_protein: str) -> Optional[str]:
    """Resolve OMA-style human_protein field to a UniProt accession.

    Accepts 'human|AKT1_HUMAN' or 'AKT1_HUMAN'; returns 'P31749'.
    ChEMBL and CanSAR require a proper UniProt accession, not an entry name.
    """
    entry_name = human_protein.split("|", 1)[-1] if "|" in human_protein else human_protein
    try:
        r = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={"query": f"id:{entry_name}", "fields": "accession", "format": "json"},
            timeout=15,
        )
        if r.status_code == 200:
            results = r.json().get("results", [])
            if results:
                return results[0].get("primaryAccession")
    except Exception as exc:
        log.debug("UniProt accession lookup %s: %s", entry_name, exc)
    return None


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

        # Use mechanism endpoint: returns drugs and their action types for a target.
        # This is the correct endpoint for "what drugs target this protein?"
        r2 = requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/mechanism.json",
            params={"target_chembl_id": chembl_id, "limit": 100},
            timeout=20,
        )
        drugs = []
        if r2.status_code == 200:
            mech_data = r2.json()
            mol_ids = list({
                m.get("molecule_chembl_id")
                for m in mech_data.get("mechanisms", [])
                if m.get("molecule_chembl_id")
            })
            # Fetch pref_names and phases for each molecule in batch
            if mol_ids:
                r3 = requests.get(
                    "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
                    params={"molecule_chembl_id__in": ",".join(mol_ids[:50]),
                            "max_phase__gte": 3, "limit": 50},
                    timeout=20,
                )
                if r3.status_code == 200:
                    for mol in r3.json().get("molecules", []):
                        name = mol.get("pref_name") or mol.get("molecule_chembl_id", "")
                        if name:
                            drugs.append(name)
        return chembl_id, list(dict.fromkeys(drugs))[:20]
    except Exception as exc:
        log.debug("ChEMBL %s: %s", uniprot_id, exc)
        return None, []


def annotate_genes_chembl(gene_ids: list[str]) -> int:
    """Populate DrugTarget.chembl_target_id and existing_drugs for the given genes.

    Resolves OMA-style human_protein field (human|AKT1_HUMAN) to UniProt
    accession (P31749) before querying ChEMBL.
    """
    updated = 0
    with get_session() as session:
        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()

    for gene in genes:
        if not gene.human_protein:
            continue
        accession = _resolve_uniprot_accession(gene.human_protein.strip())
        if not accession:
            log.debug("ChEMBL: no UniProt accession for %s", gene.human_protein)
            continue
        chembl_id, drugs = fetch_chembl_target_and_drugs(accession)
        if chembl_id is None:
            continue
        with get_session() as session:
            dt = session.get(DrugTarget, gene.id)
            if dt is None:
                dt = DrugTarget(gene_id=gene.id)
                session.add(dt)
            dt.chembl_target_id = chembl_id
            dt.existing_drugs = drugs if drugs else None
            session.commit()
        updated += 1
    log.info("ChEMBL: updated %d genes.", updated)
    return updated
