"""Step 25 — Synthesizability Scoring.

Computes two signals per compound:

  sa_score        — RDKit Synthetic Accessibility (SA) Score [1=trivial, 10=nearly impossible]
                    Ertl & Schuffenhauer 2009. Lower is better.
  sa_normalized   — (10 - sa_score) / 9, so 1.0 = trivial, 0.0 = impossible
  zinc_purchasable — True if exact SMILES found in ZINC (via InChIKey lookup)

After computing these, the final overall_compound_score and compound_tier are
calculated from the full CompoundScore profile.

compound_tier assignment:
  "Lead"     — overall ≥ 0.65
  "Hit"      — overall ≥ 0.40
  "Fragment" — overall ≥ 0.25
  "Reject"   — overall <  0.25
"""

import logging
import time
from typing import Optional

import requests

from db.models import Compound, CompoundScore
from db.session import get_session

log = logging.getLogger(__name__)

ZINC_INCHIKEY_API = "https://zinc20.docking.org/substances/subsets/fda-drugs.json"

# Weight vector for overall_compound_score
_WEIGHTS = {
    "hit_score":         0.30,   # from step 21
    "lipinski":          0.15,   # Lipinski pass
    "selectivity":       0.20,   # 1 - max off-target sim
    "tox_penalty":       0.20,   # 1 - tox21_score
    "sa_normalized":     0.15,   # synthesizability
}

_TIER_LEAD     = 0.65
_TIER_HIT      = 0.40
_TIER_FRAGMENT = 0.25


def _compute_sa_score(smiles: str) -> Optional[float]:
    """Return RDKit SA Score [1-10] or None."""
    try:
        from rdkit import Chem
        from rdkit.Chem.QED import qed as _qed  # ensure rdkit available
        from rdkit.Chem import RDConfig
        import sys
        import os
        # RDKit SA Score is in Contrib
        sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_path not in sys.path:
            sys.path.insert(0, sa_path)
        import sascorer
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return round(sascorer.calculateScore(mol), 3)
    except Exception as exc:
        log.debug("SA Score failed for %r: %s", smiles[:30], exc)
        return None


def _zinc_purchasable(smiles: str) -> bool:
    """Approximate purchasability check via ZINC InChIKey lookup."""
    try:
        from rdkit import Chem
        from rdkit.Chem.inchi import MolToInchiKey
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        inchikey = MolToInchiKey(mol)
        r = requests.get(
            f"https://zinc20.docking.org/substances/{inchikey}.json",
            timeout=8,
        )
        return r.status_code == 200
    except Exception:
        return False


def _overall_score(c: Compound, cs: CompoundScore) -> float:
    """Compute weighted overall compound score from all CompoundScore fields."""
    hit_s  = (c.hit_score or 0.0)
    lip_s  = 1.0 if cs.lipinski_pass else 0.0
    sel_s  = (cs.selectivity_score or 0.5)
    tox_s  = 1.0 - (cs.tox21_score or 0.0)
    sa_s   = (cs.sa_normalized or 0.5)

    score = (
        _WEIGHTS["hit_score"]   * hit_s
        + _WEIGHTS["lipinski"]  * lip_s
        + _WEIGHTS["selectivity"] * sel_s
        + _WEIGHTS["tox_penalty"]  * tox_s
        + _WEIGHTS["sa_normalized"] * sa_s
    )
    return round(score, 4)


def score_synthesizability(gene_ids: list[str]) -> int:
    """Populate sa_score, zinc_purchasable, and final compound_tier."""
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

        sa = _compute_sa_score(compound.smiles)
        sa_norm = round((10.0 - sa) / 9.0, 4) if sa is not None else None

        purchasable = compound.is_purchasable or _zinc_purchasable(compound.smiles)
        time.sleep(0.05)   # gentle rate limit

        with get_session() as session:
            cs = session.get(CompoundScore, compound.id)
            if cs is None:
                cs = CompoundScore(id=compound.id, compound_id=compound.id)
                session.add(cs)
                session.flush()

            cs.sa_score        = sa
            cs.sa_normalized   = sa_norm
            cs.zinc_purchasable = purchasable

            # Reload compound with fresh hit_score from session
            c = session.get(Compound, compound.id)
            overall = _overall_score(c, cs)
            cs.overall_compound_score = overall

            if overall >= _TIER_LEAD:
                tier = "Lead"
            elif overall >= _TIER_HIT:
                tier = "Hit"
            elif overall >= _TIER_FRAGMENT:
                tier = "Fragment"
            else:
                tier = "Reject"
            cs.compound_tier = tier

            session.commit()

        updated += 1

    log.info("Synthesizability: %d compounds scored; final tiers assigned.", updated)
    return updated
