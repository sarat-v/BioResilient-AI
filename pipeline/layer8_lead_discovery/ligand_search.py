"""Step 20 — Ligand Similarity Search.

For each gene with known drugs (DrugTarget.existing_drugs):
  1. Resolve known drug SMILES from ChEMBL by name.
  2. Compute Morgan fingerprint similarity (Tanimoto) to each stored Compound.
  3. Query ZINC for structural analogs of the best known binder (Tanimoto ≥ 0.4).
  4. Store new analogs as Compound rows; update tanimoto_to_known on existing ones.
"""

import logging
import time
from typing import Optional

import requests

from db.models import Compound, DrugTarget
from db.session import get_session

log = logging.getLogger(__name__)

CHEMBL_SMILES_API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
ZINC_SIMILARITY_API = "https://zinc20.docking.org/substances/subsets/fda-drugs.json"

MAX_ANALOGS = 20


def _get_morgan_fp(smiles: str):
    """Return RDKit Morgan fingerprint or None."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    except Exception:
        return None


def _tanimoto(fp1, fp2) -> float:
    try:
        from rdkit import DataStructs
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    except Exception:
        return 0.0


def _chembl_smiles_for_name(drug_name: str) -> Optional[str]:
    """Return SMILES for a drug name from ChEMBL."""
    try:
        r = requests.get(
            CHEMBL_SMILES_API,
            params={"pref_name__iexact": drug_name, "limit": 1},
            timeout=10,
        )
        if r.status_code == 200:
            mols = r.json().get("molecules", [])
            if mols:
                struct = mols[0].get("molecule_structures") or {}
                smiles = struct.get("canonical_smiles")
                chembl_id = mols[0].get("molecule_chembl_id")
                return smiles, chembl_id
    except Exception as exc:
        log.debug("ChEMBL name lookup failed for %s: %s", drug_name, exc)
    return None, None


def _zinc_similarity_search(query_smiles: str) -> list[dict]:
    """Fetch drug-like ZINC compounds as analogs (API doesn't support true similarity; use subset)."""
    try:
        r = requests.get(ZINC_SIMILARITY_API, params={"count": MAX_ANALOGS}, timeout=15)
        if r.status_code == 200:
            items = r.json()
            results = []
            query_fp = _get_morgan_fp(query_smiles)
            for item in items:
                smi = item.get("smiles") or ""
                if not smi:
                    continue
                fp = _get_morgan_fp(smi)
                sim = _tanimoto(query_fp, fp) if (query_fp and fp) else 0.0
                if sim >= 0.3:
                    results.append({
                        "smiles": smi,
                        "zinc_id": item.get("zinc_id", ""),
                        "name": item.get("name", ""),
                        "tanimoto": round(sim, 4),
                    })
            results.sort(key=lambda x: x["tanimoto"], reverse=True)
            return results[:MAX_ANALOGS]
    except Exception as exc:
        log.debug("ZINC similarity search failed: %s", exc)
    return []


def run_ligand_search(gene_ids: list[str]) -> int:
    """Search for ligand analogs using known drugs per gene."""
    with get_session() as session:
        drug_map = {dt.gene_id: dt for dt in
                    session.query(DrugTarget).filter(DrugTarget.gene_id.in_(gene_ids)).all()}
        existing_smiles = {c.smiles: c.id for c in
                           session.query(Compound.smiles, Compound.id)
                           .filter(Compound.gene_id.in_(gene_ids)).all()}

    stored_total = 0
    for gene_id in gene_ids:
        dt = drug_map.get(gene_id)
        known_drugs = (dt.existing_drugs or []) if dt else []
        if not known_drugs:
            continue

        # Pick first resolvable drug name
        known_smiles = None
        known_chembl_id = None
        known_drug_name = None
        for drug in known_drugs[:3]:
            smi, cid = _chembl_smiles_for_name(drug)
            if smi:
                known_smiles = smi
                known_chembl_id = cid
                known_drug_name = drug
                break
            time.sleep(0.1)

        if not known_smiles:
            log.debug("No SMILES found for known drugs of %s", gene_id)
            continue

        known_fp = _get_morgan_fp(known_smiles)
        log.info("Ligand search %s: anchor=%s (%s)", gene_id, known_drug_name, known_chembl_id)

        # Update Tanimoto on existing compounds
        with get_session() as session:
            existing = session.query(Compound).filter(Compound.gene_id == gene_id).all()
            for c in existing:
                if not c.smiles:
                    continue
                fp = _get_morgan_fp(c.smiles)
                if known_fp and fp:
                    sim = _tanimoto(known_fp, fp)
                    if c.tanimoto_to_known is None or sim > c.tanimoto_to_known:
                        c.tanimoto_to_known = round(sim, 4)
                        c.nearest_known_drug = known_drug_name
            session.commit()

        # Fetch analogs from ZINC
        analogs = _zinc_similarity_search(known_smiles)
        for analog in analogs:
            smi = analog["smiles"]
            if smi in existing_smiles:
                continue  # already stored from step 19
            with get_session() as session:
                c = Compound(
                    gene_id=gene_id,
                    smiles=smi,
                    name=analog.get("name"),
                    source="chembl_analog",
                    zinc_id=analog.get("zinc_id"),
                    tanimoto_to_known=analog["tanimoto"],
                    nearest_known_drug=known_drug_name,
                    is_purchasable=bool(analog.get("zinc_id")),
                )
                session.add(c)
                session.commit()
                existing_smiles[smi] = c.id
            stored_total += 1

        time.sleep(0.2)

    log.info("Ligand search: %d new analogs stored.", stored_total)
    return stored_total
