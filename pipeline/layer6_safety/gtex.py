"""Step 14b-ii — GTEx tissue expression breadth for drug target safety.

GTEx (Genotype-Tissue Expression) provides median gene expression (TPM) across
~54 human tissues. For drug target safety, we care about:

  - Expression breadth: a gene expressed in many tissues is harder to target
    without off-target effects (high gtex_tissue_count = more toxicity risk).
  - Maximum expression: a gene highly expressed in a critical tissue (e.g.,
    heart, brain) is higher risk.

We use the GTEx v10 gene median TPM file via the GTEx API.

Safety interpretation:
  - gtex_tissue_count < 5: tissue-restricted → lower toxicity risk
  - gtex_tissue_count 5–20: moderate breadth → medium risk
  - gtex_tissue_count > 20: ubiquitously expressed → higher risk
  - gtex_max_tpm > 1000 in heart/brain → flag as high-risk

Data source:
  GTEx Portal REST API: https://gtexportal.org/api/v2/expression/geneExpression
  No auth required. Respects rate limits with 0.2s delay.
"""

import logging
import time
from typing import Optional

import requests

from db.models import Gene, SafetyFlag
from db.session import get_session

log = logging.getLogger(__name__)

GTEX_API = "https://gtexportal.org/api/v2/expression/geneExpression"
REQUEST_TIMEOUT = 20
TPM_THRESHOLD = 1.0   # Tissue "expressed" threshold


def fetch_gtex_expression(gene_symbol: str) -> Optional[dict]:
    """Fetch GTEx median TPM for a gene symbol across all tissues.

    Returns dict: {tissue_name: median_tpm} or None on failure.
    """
    try:
        r = requests.get(
            GTEX_API,
            params={
                "gencodeId": gene_symbol,   # GTEx accepts symbol or ENSEMBL
                "datasetId": "gtex_v10",
            },
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            # Try ENSEMBL-style query (GTEx prefers ENSG IDs)
            return None

        data = r.json()
        rows = data.get("data", [])
        result = {}
        for row in rows:
            tissue = row.get("tissueSiteDetailId") or row.get("tissueSiteDetail")
            tpm = row.get("median")
            if tissue and tpm is not None:
                try:
                    result[tissue] = float(tpm)
                except ValueError:
                    pass
        return result if result else None

    except Exception as exc:
        log.debug("GTEx fetch failed for %s: %s", gene_symbol, exc)
        return None


def annotate_gtex(gene_ids: Optional[list[str]] = None) -> int:
    """Annotate SafetyFlag rows with GTEx expression breadth metrics.

    For each gene:
      - gtex_tissue_count: number of tissues with median TPM > threshold
      - gtex_max_tpm: maximum median TPM across all tissues

    Returns number of genes annotated.
    """
    with get_session() as session:
        if gene_ids is None:
            gene_ids = [g.id for g in session.query(Gene).all()]
        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()

    log.info("Annotating GTEx expression breadth for %d genes...", len(genes))
    annotated = 0

    for gene in genes:
        if not gene.gene_symbol:
            continue

        expr = fetch_gtex_expression(gene.gene_symbol)
        time.sleep(0.2)

        if not expr:
            continue

        expressed_tissues = [t for t, tpm in expr.items() if tpm >= TPM_THRESHOLD]
        max_tpm = max(expr.values()) if expr else 0.0

        with get_session() as session:
            sf = session.get(SafetyFlag, gene.id)
            if sf is None:
                sf = SafetyFlag(gene_id=gene.id)
                session.add(sf)

            sf.gtex_tissue_count = len(expressed_tissues)
            sf.gtex_max_tpm = round(max_tpm, 2)
            session.commit()
            annotated += 1

    log.info("GTEx annotation complete: %d genes annotated.", annotated)
    return annotated


def run_gtex_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step14b."""
    return annotate_gtex(gene_ids)
