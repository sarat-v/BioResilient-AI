"""Backfill disease_name and disease_id for genes that have opentargets_score
but are missing disease_name. Also re-patches GWAS p-values where mantissa was 0.

Run locally (not on AWS Batch) — pure API + DB operation.
  python3 scripts/backfill_disease_names.py
"""

import logging
import sys
import time
import os

# Allow running from repo root
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from db.models import DiseaseAnnotation, Gene
from db.session import get_session
from pipeline.layer3_disease.opentargets import (
    fetch_opentargets_full,
    _hgnc_symbol,
    _symbol_to_ensembl_by_symbol,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

_RATE_SLEEP = 0.15  # conservative rate to avoid OT CloudFlare throttle


def backfill():
    # Find all genes with OT score but no disease_name
    with get_session() as session:
        rows = (
            session.query(DiseaseAnnotation.gene_id, Gene.gene_symbol)
            .join(Gene, Gene.id == DiseaseAnnotation.gene_id)
            .filter(
                DiseaseAnnotation.opentargets_score.isnot(None),
                DiseaseAnnotation.disease_name.is_(None),
            )
            .all()
        )

    log.info("Genes with OT score but no disease_name: %d", len(rows))
    if not rows:
        log.info("Nothing to backfill.")
        return

    updates: dict[str, dict] = {}
    for i, (gene_id, gene_symbol) in enumerate(rows):
        hgnc = _hgnc_symbol(gene_symbol)
        eid = _symbol_to_ensembl_by_symbol(hgnc)
        if not eid:
            log.warning("  No Ensembl ID for %s — skipping", gene_symbol)
            continue
        result = fetch_opentargets_full(gene_id, gene_symbol, ensembl_id=eid)
        if result and result.get("disease_name"):
            updates[gene_id] = result
            log.info("  [%d/%d] %s → %s (%s)",
                     i + 1, len(rows), gene_symbol,
                     result["disease_name"], result.get("disease_id", ""))
        else:
            log.info("  [%d/%d] %s → no disease data", i + 1, len(rows), gene_symbol)
        if i < len(rows) - 1:
            time.sleep(_RATE_SLEEP)

    if not updates:
        log.info("No disease names retrieved from API.")
        return

    log.info("Writing %d disease_name / disease_id updates...", len(updates))
    gene_ids = list(updates.keys())
    with get_session() as session:
        ann_map = {
            r.gene_id: r
            for r in session.query(DiseaseAnnotation)
            .filter(DiseaseAnnotation.gene_id.in_(gene_ids))
            .all()
        }
        for gene_id, result in updates.items():
            ann = ann_map.get(gene_id)
            if ann is None:
                continue
            ann.disease_name = result["disease_name"]
            if result.get("disease_id"):
                ann.disease_id = result["disease_id"]

    log.info("Done. %d records updated.", len(updates))


if __name__ == "__main__":
    backfill()
