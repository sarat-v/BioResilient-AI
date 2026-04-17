"""CRISPR tractability — guide RNA sites, on-target efficiency (Rule Set 1 proxy),
and off-target risk assessment.

When no DNA sequence is available, crispr_sites and offtarget_risk are stored as NULL
(not 0 / 'unknown') so that downstream callers can distinguish "not analysed" from
"analysed and found 0 sites". Storing fabricated zeros would silently bias scoring.

Guide scoring uses a Rule Set 1–inspired proxy (GC content 40–70%, seed region quality,
no poly-T) to estimate on-target efficiency without requiring CRISPick API access.
"""

import logging
import re
from typing import Optional

from db.models import GeneTherapyScore
from db.session import get_session

log = logging.getLogger(__name__)

# NGG PAM pattern for SpCas9
PAM_PATTERN = re.compile(r"(?=([ACGT]{20}GG))", re.I)

# Bases associated with lower efficiency in the seed region (positions 1–4, 3′ end)
_SEED_UNFAVORABLE = {"T"}


def _rule_set1_score(guide: str) -> float:
    """Lightweight Rule Set 1–inspired on-target efficiency score (0-1).

    Based on Doench et al. 2014 heuristics:
      - GC content 40–70%: +0.4
      - No T in seed positions (guide[16:20]): +0.2
      - No poly-T (TTTT) anywhere: +0.2
      - Purine at guide[19] (PAM-proximal): +0.1
      - No G at guide[0] (disfavored): +0.1
    """
    g = guide.upper()
    if len(g) < 20:
        return 0.0

    score = 0.0
    gc = (g.count("G") + g.count("C")) / len(g)
    if 0.40 <= gc <= 0.70:
        score += 0.4

    seed = g[16:20]
    if "T" not in seed:
        score += 0.2

    if "TTTT" not in g:
        score += 0.2

    if g[19] in ("A", "G"):   # purine at PAM-proximal
        score += 0.1

    if g[0] != "G":
        score += 0.1

    return round(min(score, 1.0), 3)


def count_pam_sites(seq: str, max_len: int = 10000) -> int:
    """Count NGG PAM sites in a DNA sequence."""
    if not seq or len(seq) > max_len:
        return 0
    return len(PAM_PATTERN.findall(seq.upper()))


def score_guides(seq: str, max_guides: int = 50) -> tuple[int, float, str]:
    """Score guides in a sequence and return (count, mean_efficiency, risk_level).

    Returns:
        (sites, mean_efficiency, risk_level) where:
          - sites: number of valid NGG PAM guides found
          - mean_efficiency: mean Rule Set 1 score across top guides (0-1)
          - risk_level: 'low' / 'medium' / 'high' based on site count
    """
    seq_upper = seq.upper()
    guides = PAM_PATTERN.findall(seq_upper)
    # Score all guides, keep top max_guides by efficiency
    scored = sorted(
        [(_rule_set1_score(g[:20]), g[:20]) for g in guides if len(g) >= 20],
        reverse=True,
    )[:max_guides]

    sites = len(guides)
    mean_eff = round(sum(s for s, _ in scored) / len(scored), 3) if scored else 0.0

    if sites > 50:
        risk = "high"
    elif sites > 15:
        risk = "medium"
    else:
        risk = "low"

    return sites, mean_eff, risk


def run_crispor_offtarget(seq_region: str) -> tuple[int, str]:
    """Backward-compat wrapper. Returns (sites, risk_level)."""
    sites, _eff, risk = score_guides(seq_region)
    return sites, risk


def annotate_genes_crispr(gene_ids: list[str], sequence_provider=None) -> int:
    """Populate GeneTherapyScore.crispr_sites and offtarget_risk.

    sequence_provider: optional callable gene_id -> DNA sequence string (e.g. from NCBI).

    When sequence_provider is None or returns an empty string, crispr_sites and
    offtarget_risk are stored as NULL — meaning "not analysed", which is scientifically
    honest and distinguishable from a real result of 0 sites found.
    """
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gs = session.get(GeneTherapyScore, gid)
            if gs is None:
                gs = GeneTherapyScore(gene_id=gid)
                session.add(gs)
            seq = (sequence_provider(gid) if sequence_provider else None) or ""
            if seq:
                sites, mean_eff, risk = score_guides(seq)
                gs.crispr_sites = sites
                gs.offtarget_risk = risk
                gs.crispr_efficiency = mean_eff
                log.info("CRISPR %s: sites=%d  eff=%.3f  risk=%s", gid, sites, mean_eff, risk)
            else:
                # No sequence — store NULL rather than fabricating 0 / 'unknown'.
                gs.crispr_sites = None
                gs.offtarget_risk = None
                gs.crispr_efficiency = None
                log.info("CRISPR %s: no sequence provider — storing NULL (not analysed).", gid)
            updated += 1
        session.commit()
    log.info("CRISPR: updated %d genes.", updated)
    return updated
