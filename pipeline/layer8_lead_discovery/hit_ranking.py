"""Step 21 — Hit Ranking.

Aggregates docking_score, is_purchasable, and tanimoto_to_known into a
single hit_score per compound, then assigns hit_tier.

Scoring formula:
  docking_component  = normalised docking score within gene's compound set
  novelty_component  = 1 - tanimoto_to_known (closer to known drug = less novel but safer)
  purchasability     = 1.0 if purchasable else 0.4
  hit_score = 0.5 * docking_component + 0.3 * purchasability + 0.2 * novelty_component

Tiers:
  "A" — hit_score ≥ 0.70  (top leads)
  "B" — hit_score ≥ 0.45  (hits for follow-up)
  "C" — hit_score <  0.45  (fragments / low confidence)
"""

import logging

from db.models import Compound
from db.session import get_session

log = logging.getLogger(__name__)

_TIER_A = 0.70
_TIER_B = 0.45


def _normalise(values: list[float], higher_is_better: bool = True) -> list[float]:
    """Min-max normalise a list. Returns [0, 1]."""
    if not values:
        return []
    mn, mx = min(values), max(values)
    if mx == mn:
        return [0.5] * len(values)
    if higher_is_better:
        return [(v - mn) / (mx - mn) for v in values]
    return [(mx - v) / (mx - mn) for v in values]


def rank_hits(gene_ids: list[str]) -> int:
    """Compute hit_score and hit_tier for all compounds of these genes."""
    ranked = 0
    for gene_id in gene_ids:
        with get_session() as session:
            compounds = (
                session.query(Compound)
                .filter(Compound.gene_id == gene_id)
                .all()
            )
            if not compounds:
                continue

            # Docking scores — lower (more negative) is better for Vina; higher for DiffDock
            raw_scores = []
            for c in compounds:
                if c.docking_score is not None:
                    raw_scores.append(c.docking_score)
                else:
                    raw_scores.append(0.0)

            # Determine if higher or lower docking_score is better
            method = next((c.docking_method for c in compounds if c.docking_method), "estimated")
            higher_better = method in ("diffdock", "estimated")
            norm_dock = _normalise(raw_scores, higher_is_better=higher_better)

            for i, c in enumerate(compounds):
                dock_comp = norm_dock[i]
                purch_comp = 1.0 if c.is_purchasable else 0.4
                tan = c.tanimoto_to_known or 0.0
                novelty = 1.0 - tan   # 0 = identical to known drug; 1 = completely novel

                score = round(
                    0.50 * dock_comp
                    + 0.30 * purch_comp
                    + 0.20 * novelty,
                    4,
                )
                tier = "A" if score >= _TIER_A else ("B" if score >= _TIER_B else "C")
                c.hit_score = score
                c.hit_tier  = tier
                ranked += 1

            session.commit()

    log.info("Hit ranking: %d compounds ranked.", ranked)
    return ranked
