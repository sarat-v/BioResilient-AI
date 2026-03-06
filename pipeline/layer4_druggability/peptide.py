"""Peptide tractability — Boman index, length, synthesisable flag for DivergentMotif."""

import logging
from typing import Optional

from db.models import DivergentMotif, Ortholog
from db.session import get_session

log = logging.getLogger(__name__)

# Boman index (kcal/mol) for amino acids — used to estimate stability
BOMAN = {
    "L": 4.92, "I": 4.92, "V": 4.04, "M": 2.35, "W": 2.33, "F": 2.98,
    "A": 1.81, "C": 1.28, "G": 0.94, "T": 2.57, "S": 1.60, "Y": 2.00,
    "P": 2.67, "H": 2.33, "E": 1.56, "Q": 1.46, "D": 1.23, "N": 0.84,
    "K": 2.80, "R": 1.76, "X": 1.0,
}

PEPTIDE_LEN_MIN = 10
PEPTIDE_LEN_MAX = 20


def boman_index(seq: str) -> float:
    """Compute Boman index (average residue potential) for a peptide. Higher = more likely soluble/stable."""
    if not seq:
        return 0.0
    total = sum(BOMAN.get(aa.upper(), 1.0) for aa in seq)
    return round(total / len(seq), 4)


def estimate_half_life_min(boman: float, length: int) -> Optional[float]:
    """Rough estimate of half-life (minutes) from Boman and length. Heuristic."""
    if length < 5:
        return None
    # Higher Boman and shorter length → more stable
    base = 30.0
    factor = (boman / 2.0) * (1.0 + (PEPTIDE_LEN_MAX - min(length, PEPTIDE_LEN_MAX)) / 20.0)
    return round(base * factor, 1)


def is_synthesisable(animal_seq: str, human_seq: str, length_ok: bool = True) -> bool:
    """True if motif is in tractable length and not obviously problematic."""
    seq = (animal_seq or human_seq or "").replace("-", "")
    if len(seq) < PEPTIDE_LEN_MIN or len(seq) > PEPTIDE_LEN_MAX:
        return False
    if not length_ok:
        return False
    # Exclude if too many cysteines (disulfide complexity)
    if seq.count("C") + seq.count("c") > 3:
        return False
    return True


def annotate_motifs_peptide(motif_ids: Optional[list[str]] = None) -> int:
    """Populate DivergentMotif.half_life_min and synthesisable from sequence."""
    updated = 0
    with get_session() as session:
        q = session.query(DivergentMotif)
        if motif_ids is not None:
            q = q.filter(DivergentMotif.id.in_(motif_ids))
        for m in q:
            seq = (m.animal_seq or m.human_seq or "").replace("-", "")
            length_ok = PEPTIDE_LEN_MIN <= len(seq) <= PEPTIDE_LEN_MAX
            m.synthesisable = is_synthesisable(m.animal_seq or "", m.human_seq or "", length_ok)
            b = boman_index(seq)
            m.half_life_min = estimate_half_life_min(b, len(seq)) if seq else None
            updated += 1
    log.info("Peptide tractability: updated %d motifs.", updated)
    return updated
