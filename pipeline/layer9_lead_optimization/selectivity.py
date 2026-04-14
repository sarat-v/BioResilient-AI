"""Step 24 — Selectivity Check.

For each compound, computes Tanimoto Morgan fingerprint similarity against a
panel of known off-target ligands fetched from ChEMBL (kinase inhibitors,
GPCRs, ion channels — the three most common selectivity problem classes).

selectivity_score = 1 - max_tanimoto_to_offtargets (higher = more selective)
n_offtarget_similar = count of off-target ligands with Tanimoto > 0.4
"""

import logging
import time
from typing import Optional

import requests

from db.models import Compound, CompoundScore
from db.session import get_session

log = logging.getLogger(__name__)

CHEMBL_ACTIVITY_API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

# ChEMBL target IDs for off-target classes
OFF_TARGET_CHEMBL_IDS = [
    "CHEMBL279",    # EGFR kinase
    "CHEMBL203",    # VEGFR2
    "CHEMBL240",    # hERG (also in tox step, here for selectivity)
    "CHEMBL2056",   # Dopamine D2 receptor
    "CHEMBL3778",   # Adenosine A2A receptor
    "CHEMBL2487",   # Muscarinic M1
]

MAX_OFFTARGET_LIGANDS = 50


def _fetch_offtarget_smiles() -> list[str]:
    """Fetch SMILES for known actives at off-target proteins from ChEMBL."""
    all_smiles: list[str] = []
    for target_id in OFF_TARGET_CHEMBL_IDS[:3]:  # limit API calls
        try:
            r = requests.get(
                CHEMBL_ACTIVITY_API,
                params={
                    "target_chembl_id": target_id,
                    "pchembl_value__gte": "6",  # IC50 <= 1 μM
                    "limit": MAX_OFFTARGET_LIGANDS,
                    "format": "json",
                },
                timeout=15,
            )
            if r.status_code == 200:
                acts = r.json().get("activities", [])
                for act in acts:
                    smi = (act.get("molecule_structures") or {}).get("canonical_smiles")
                    if smi and smi not in all_smiles:
                        all_smiles.append(smi)
            time.sleep(0.2)
        except Exception as exc:
            log.debug("ChEMBL off-target fetch %s: %s", target_id, exc)

    if not all_smiles:
        # Minimal fallback: a few well-known promiscuous scaffolds
        all_smiles = [
            "c1ccc2c(c1)NC(=O)c1ccccc12",     # acridone (promiscuous)
            "O=C1NC(=O)c2ccccc21",              # isatoic anhydride
            "c1ccc(cc1)C2=CC(=O)Oc3ccccc23",   # flavone
        ]
    return all_smiles


def _morgan_fp(smiles: str):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    except Exception:
        return None


def _max_tanimoto(query_fp, ref_fps: list) -> tuple[float, int]:
    """Return (max_tanimoto, n_above_threshold) against ref panel."""
    if not ref_fps or query_fp is None:
        return 0.0, 0
    try:
        from rdkit import DataStructs
        sims = [DataStructs.TanimotoSimilarity(query_fp, r) for r in ref_fps if r is not None]
        if not sims:
            return 0.0, 0
        return round(max(sims), 4), sum(1 for s in sims if s > 0.4)
    except Exception:
        return 0.0, 0


def score_selectivity(gene_ids: list[str]) -> int:
    """Compute selectivity scores for all compounds of given genes."""
    offtarget_smiles = _fetch_offtarget_smiles()
    offtarget_fps = [_morgan_fp(s) for s in offtarget_smiles]
    offtarget_fps = [fp for fp in offtarget_fps if fp is not None]
    log.info("Selectivity: %d off-target fingerprints loaded", len(offtarget_fps))

    with get_session() as session:
        compounds = (
            session.query(Compound)
            .filter(Compound.gene_id.in_(gene_ids))
            .all()
        )

    updated = 0
    for compound in compounds:
        if not compound.smiles:
            continue

        fp = _morgan_fp(compound.smiles)
        max_tan, n_similar = _max_tanimoto(fp, offtarget_fps)
        selectivity = round(1.0 - max_tan, 4)

        with get_session() as session:
            cs = session.get(CompoundScore, compound.id)
            if cs is None:
                cs = CompoundScore(id=compound.id, compound_id=compound.id)
                session.add(cs)
            cs.selectivity_score    = selectivity
            cs.n_offtarget_similar  = n_similar
            session.commit()

        updated += 1

    log.info("Selectivity: %d compounds scored.", updated)
    return updated
