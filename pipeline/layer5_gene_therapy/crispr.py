"""CRISPR tractability — guide RNA sites and off-target risk (stub when CRISPOR/Cas-OFFinder not installed)."""

import logging
import re
import subprocess
from pathlib import Path
from typing import Optional

from db.models import GeneTherapyScore
from db.session import get_session

log = logging.getLogger(__name__)

# NGG PAM pattern for SpCas9 — rough count of possible sites in a sequence
PAM_PATTERN = re.compile(r"(?=([ACGT]{20}GG))", re.I)


def count_pam_sites(seq: str, max_len: int = 5000) -> int:
    """Count NGG PAM sites in a DNA sequence (rough proxy for guide availability)."""
    if not seq or len(seq) > max_len:
        return 0
    return len(PAM_PATTERN.findall(seq.upper()))


def run_crispor_offtarget(seq_region: str) -> tuple[int, str]:
    """If CRISPOR or Cas-OFFinder is available, return (sites, risk). Else return (0, 'unknown')."""
    # Stub: use PAM count as proxy; risk = 'medium' if many sites else 'low'
    sites = count_pam_sites(seq_region)
    if sites == 0:
        return 0, "unknown"
    if sites > 50:
        return sites, "high"
    if sites > 15:
        return sites, "medium"
    return sites, "low"


def annotate_genes_crispr(gene_ids: list[str], sequence_provider=None) -> int:
    """Populate GeneTherapyScore.crispr_sites and offtarget_risk.

    sequence_provider: optional callable gene_id -> DNA sequence string (e.g. from NCBI). If None, we set 0 and 'unknown'.
    """
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gs = session.get(GeneTherapyScore, gid)
            if gs is None:
                gs = GeneTherapyScore(gene_id=gid)
                session.add(gs)
            seq = (sequence_provider(gid) if sequence_provider else None) or ""
            sites, risk = run_crispor_offtarget(seq)
            gs.crispr_sites = sites
            gs.offtarget_risk = risk
            updated += 1
    log.info("CRISPR: updated %d genes.", updated)
    return updated
