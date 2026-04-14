"""Step 23 — Toxicity Screening.

Computes two toxicity signals per compound:

  tox21_score    — Aggregate structural alert score covering the 12 Tox21 endpoint
                   categories (NR-AR, NR-AhR, SR-ARE, etc.) using curated SMARTS
                   patterns from literature. Each alert contributes 1/12 to score.
  herg_risk      — hERG K+ channel liability estimate from known pharmacophore
                   patterns (basic aromatic amine + positive charge predictor).
  pains_alerts   — Count of Pan Assay Interference Compounds (PAINS) filters.
  structural_alerts — Count of standard medicinal chemistry structural alerts
                      (Murcko / Brenk / SMARTS catalogue).

All logic uses RDKit SMARTS only — no model weights required.
"""

import logging
from typing import Optional

from db.models import Compound, CompoundScore
from db.session import get_session

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Tox21 structural alert SMARTS (one representative pattern per endpoint)
# These cover the most common scaffolds; not exhaustive.
# ---------------------------------------------------------------------------
_TOX21_SMARTS = {
    "NR-AR":     "[#6]-c1ccc(cc1)O",                        # phenol (androgen receptor)
    "NR-AhR":    "c1ccc2ccccc2c1",                          # naphthalene (AhR ligand)
    "NR-ER":     "c1ccc(cc1)O",                              # 4-hydroxyphenyl (ER)
    "NR-ER-LBD": "c1ccc2c(c1)ccc(c2)O",                    # naphthol (ER LBD)
    "NR-PPAR-g": "c1ccc(cc1)OCC(=O)O",                     # phenoxy acetic acid
    "SR-ARE":    "[N;X3]([#6])[#6]=O",                      # N-acylamine (ARE)
    "SR-ATAD5":  "O=C1NC(=O)NC(=O)1",                      # barbituric acid
    "SR-HSE":    "c1ccc(cc1)Cl",                            # chlorobenzene (HSE stress)
    "SR-MMP":    "N(=O)=O",                                 # nitro (MMP)
    "SR-p53":    "[Hg,Cd,As,Pb,Cr,Ni]",                    # heavy metal
    "NR-Ar-LBD": "c1ccncc1",                                # pyridine scaffold
    "NR-AR-LBD": "[#6]-C1=CC(=O)CC[C@@H]1C",               # cyclohexenone
}

# hERG pharmacophore alerts
_HERG_SMARTS = [
    "c1ccccc1CC[NH2+]",                 # phenethylamine
    "[NH+]1CCCCC1",                     # protonated piperidine
    "c1ccccc1C(F)(F)F",                 # trifluoromethylbenzene
    "c1ccc(cc1)C2=NCCN2",              # benzimidazole
]

# PAINS filters (curated subset of Baell & Holloway 2010 PAINS-A)
_PAINS_SMARTS = [
    "[#6]1(=O)[#6]=,:[#6][#6](=O)[#6]=,:[#6]1",  # quinone
    "O=C1c2ccccc2C(=O)c3ccccc13",                  # anthraquinone
    "c1cc([OH])ccc1S(=O)(=O)N",                    # sulfonamide phenol
    "[nH]1cccc1-c1ccccc1",                          # pyrrole biaryl
    "S1C=CC(=N1)c1ccccc1",                          # thiazolidine
]

# Standard medicinal chemistry structural alerts (Brenk 2008 / SMARTS catalogue)
_STRUCT_ALERTS = [
    "[CX4][F,Cl,Br,I]",               # alkyl halide
    "O=CS",                             # thioester
    "N=[N+]=[N-]",                      # azide
    "C=C-C(=O)O",                       # Michael acceptor
    "[c,n]1ccccc1N",                    # aniline
    "O=N-c1ccccc1",                     # nitrosobenzene
    "CC#N",                             # nitrile
    "S(=O)(=O)Cl",                      # sulfonyl chloride
]


def _count_smarts_matches(smiles: str, smarts_list: list[str]) -> int:
    """Return number of distinct SMARTS patterns that match the molecule."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0
        count = 0
        for smarts in smarts_list:
            try:
                pat = Chem.MolFromSmarts(smarts)
                if pat and mol.HasSubstructMatch(pat):
                    count += 1
            except Exception:
                pass
        return count
    except ImportError:
        return 0
    except Exception:
        return 0


def _tox21_score(smiles: str) -> float:
    """Return Tox21 aggregate alert score [0-1]."""
    hits = _count_smarts_matches(smiles, list(_TOX21_SMARTS.values()))
    return round(hits / len(_TOX21_SMARTS), 4)


def _herg_risk(smiles: str) -> float:
    """Return hERG liability estimate [0-1]."""
    hits = _count_smarts_matches(smiles, _HERG_SMARTS)
    return round(min(hits * 0.25, 1.0), 4)


def score_toxicity(gene_ids: list[str]) -> int:
    """Populate CompoundScore tox fields for all compounds of given genes."""
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

        tox21  = _tox21_score(compound.smiles)
        herg   = _herg_risk(compound.smiles)
        pains  = _count_smarts_matches(compound.smiles, _PAINS_SMARTS)
        struct = _count_smarts_matches(compound.smiles, _STRUCT_ALERTS)

        with get_session() as session:
            cs = session.get(CompoundScore, compound.id)
            if cs is None:
                cs = CompoundScore(id=compound.id, compound_id=compound.id)
                session.add(cs)
            cs.tox21_score       = tox21
            cs.herg_risk         = herg
            cs.pains_alerts      = pains
            cs.structural_alerts = struct
            session.commit()

        updated += 1

    log.info("Toxicity: %d compounds scored.", updated)
    return updated
