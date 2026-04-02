"""Statistical utilities shared across pipeline layers.

Currently provides:
  - apply_bh_correction(): Benjamini-Hochberg FDR correction over any list of p-values
  - apply_bh_to_meme(): convenience wrapper to BH-correct MEME site p-values across all orthogroups
  - apply_bh_to_evolution_scores(): update EvolutionScore.meme_qvalue in the DB after scoring
"""

import logging
import math
from typing import Optional

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Core BH implementation (no scipy dependency)
# ---------------------------------------------------------------------------

def apply_bh_correction(pvalues: list[float]) -> list[float]:
    """Return Benjamini-Hochberg adjusted q-values for a list of raw p-values.

    Handles NaN / None by preserving them in-place.  Returns a list of the
    same length where each value is the BH-adjusted q-value (clamped to [0,1]).

    The BH procedure:
        rank p-values ascending; q_i = p_i * n / rank_i;
        then enforce monotonicity from right to left.
    """
    n = len(pvalues)
    if n == 0:
        return []

    indexed = [(i, p) for i, p in enumerate(pvalues)]
    valid = [(i, p) for i, p in indexed if p is not None and not math.isnan(p)]
    invalid_indices = {i for i, p in indexed if p is None or math.isnan(p)}

    if not valid:
        return [p for _, p in indexed]

    valid_sorted = sorted(valid, key=lambda x: x[1])
    m = len(valid_sorted)

    qvalues_sorted = [None] * m
    for rank_0based, (orig_idx, p) in enumerate(valid_sorted):
        rank = rank_0based + 1
        qvalues_sorted[rank_0based] = p * m / rank

    # Enforce monotonicity (step-down adjustment)
    running_min = qvalues_sorted[-1]
    for k in range(m - 2, -1, -1):
        running_min = min(running_min, qvalues_sorted[k])
        qvalues_sorted[k] = running_min

    # Clamp to [0, 1]
    qvalues_sorted = [min(1.0, max(0.0, q)) for q in qvalues_sorted]

    # Scatter back to original positions
    result: list[Optional[float]] = [None] * n
    for rank_0based, (orig_idx, _) in enumerate(valid_sorted):
        result[orig_idx] = qvalues_sorted[rank_0based]

    return result


# ---------------------------------------------------------------------------
# MEME site-level FDR across orthogroups
# ---------------------------------------------------------------------------

def apply_bh_to_meme_results(meme_results: dict[str, dict]) -> dict[str, dict]:
    """Apply BH correction to MEME site p-values across ALL orthogroups jointly.

    Args:
        meme_results: {og_id: {"sites": [{"pvalue": float, ...}, ...], ...}}

    Returns:
        Same structure with "qvalue" added to each site dict.
    """
    # Flatten all site p-values with (og_id, site_index) keys
    flat: list[tuple[str, int, float]] = []
    for og_id, result in meme_results.items():
        for site_idx, site in enumerate(result.get("sites", [])):
            p = site.get("pvalue")
            if p is not None:
                flat.append((og_id, site_idx, p))

    if not flat:
        return meme_results

    pvalues_only = [p for _, _, p in flat]
    qvalues = apply_bh_correction(pvalues_only)

    for i, (og_id, site_idx, _) in enumerate(flat):
        meme_results[og_id]["sites"][site_idx]["qvalue"] = qvalues[i]

    n_sig = sum(1 for q in qvalues if q is not None and q < 0.05)
    log.info("MEME BH correction: %d / %d sites remain significant at q < 0.05", n_sig, len(flat))
    return meme_results


# ---------------------------------------------------------------------------
# DB-level FDR: update EvolutionScore.meme_qvalue after scoring
# ---------------------------------------------------------------------------

def apply_bh_to_evolution_scores() -> None:
    """Compute BH-adjusted q-values for all EvolutionScore rows and store them.

    Corrects the gene-level MEME p-value (dnds_pvalue) across all genes using
    BH FDR. The corrected q-value is stored in meme_qvalue for downstream
    tier assignment (q < 0.05 for Tier 1).
    """
    from db.models import EvolutionScore
    from db.session import get_session

    with get_session() as session:
        rows = session.query(EvolutionScore).all()
        if not rows:
            log.info("No EvolutionScore rows to BH-correct.")
            return

        pvalues = [r.dnds_pvalue for r in rows]
        qvalues = apply_bh_correction(pvalues)

        updated = 0
        for row, q in zip(rows, qvalues):
            if q is not None:
                row.meme_qvalue = round(q, 6)
                updated += 1
        session.commit()

    log.info("BH-corrected %d EvolutionScore rows (dnds_pvalue → meme_qvalue).", updated)
