"""Step 14b-ii — GTEx tissue expression breadth for drug target safety.

GTEx (Genotype-Tissue Expression) provides median gene expression (TPM) across
~54 human tissues. For drug target safety, we care about:

  - Expression breadth: a gene expressed in many tissues is harder to target
    without off-target effects (high gtex_tissue_count = more toxicity risk).
  - Maximum expression: a gene highly expressed in a critical tissue (e.g.,
    heart, brain) is higher risk.

We use the GTEx v10 gene median TPM file via the GTEx Portal REST API.

Safety interpretation:
  - gtex_tissue_count < 5: tissue-restricted → lower toxicity risk
  - gtex_tissue_count 5–20: moderate breadth → medium risk
  - gtex_tissue_count > 20: ubiquitously expressed → higher risk
  - gtex_max_tpm > 1000 in heart/brain → flag as high-risk

Data source:
  GTEx Portal REST API v2: https://gtexportal.org/api/v2/expression/geneExpression
  No auth required. Gene symbol query via 'geneId' parameter.
"""

import json
import logging
import statistics
import time
from pathlib import Path
from typing import Optional

import requests

from db.models import Gene, SafetyFlag
from db.session import get_session
from pipeline.config import get_local_storage_root

log = logging.getLogger(__name__)

GTEX_GENE_SEARCH = "https://gtexportal.org/api/v2/reference/gene"
GTEX_EXPRESSION  = "https://gtexportal.org/api/v2/expression/geneExpression"
REQUEST_TIMEOUT  = 20
TPM_THRESHOLD    = 1.0   # Tissue "expressed" threshold

# ---------------------------------------------------------------------------
# Shared GTEx expression disk cache
# ---------------------------------------------------------------------------
# aav.py (step 13) and gtex.py (step 14b) both call fetch_gtex_expression for
# the same gene symbols.  Caching to disk eliminates the duplicate API calls
# (~54 HTTP requests per gene → only done once per gene per run).

_GTEX_EXPR_CACHE: dict[str, Optional[dict]] = {}
_GTEX_EXPR_CACHE_DIRTY = False


def _gtex_cache_path() -> Path:
    return Path(get_local_storage_root()) / "gtex_expression_cache.json"


def _load_gtex_cache() -> None:
    p = _gtex_cache_path()
    if p.exists():
        try:
            _GTEX_EXPR_CACHE.update(json.loads(p.read_text()))
        except Exception:
            pass


def _save_gtex_cache() -> None:
    if not _GTEX_EXPR_CACHE_DIRTY:
        return
    p = _gtex_cache_path()
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(_GTEX_EXPR_CACHE, indent=2))
    except Exception:
        pass


def _hgnc_symbol(raw: str) -> str:
    """Convert OMA entry name (AKT1_HUMAN) to bare HGNC symbol (AKT1)."""
    if raw and "_" in raw:
        return raw.split("_")[0]
    return raw or ""


def _get_gencode_id(gene_symbol: str) -> Optional[str]:
    """Resolve HGNC gene symbol → Gencode ENSG ID via GTEx gene search API.

    GTEx expression endpoint requires a versioned gencode ID (e.g. ENSG00000142208.17).
    We get it from the GTEx /reference/gene endpoint which accepts gene symbols.
    """
    try:
        r = requests.get(
            GTEX_GENE_SEARCH,
            params={"geneId": gene_symbol, "pageSize": 1},
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            log.warning("GTEx gene search %s: HTTP %d", gene_symbol, r.status_code)
            return None
        data = r.json()
        genes = data.get("data", [])
        if not genes:
            log.info("GTEx gene search %s: no result.", gene_symbol)
            return None
        return genes[0].get("gencodeId")
    except Exception as exc:
        log.warning("GTEx gene search %s: %s", gene_symbol, exc)
        return None


def fetch_gtex_expression(gene_symbol: str) -> Optional[dict]:
    """Fetch GTEx median TPM for a gene symbol across all tissues.

    Results are cached to disk so step 13 (aav.py tissue tropism) and
    step 14b (gtex.py safety breadth) share results without duplicate HTTP calls.

    Returns dict: {tissue_name: median_tpm} or None on failure.
    Uses gtex_v8 (gtex_v10 returns empty for most genes in the v2 API).

    Response format: list of rows, each with:
      - 'tissueSiteDetailId': tissue name
      - 'data': list of per-sample TPM floats (we compute median)
    """
    if not _GTEX_EXPR_CACHE:
        _load_gtex_cache()

    key = gene_symbol.upper()
    if key in _GTEX_EXPR_CACHE:
        cached = _GTEX_EXPR_CACHE[key]
        log.debug("GTEx cache hit for %s", gene_symbol)
        return cached  # may be None (previously unsuccessful lookup)

    gencode_id = _get_gencode_id(gene_symbol)
    if not gencode_id:
        global _GTEX_EXPR_CACHE_DIRTY
        _GTEX_EXPR_CACHE[key] = None
        _GTEX_EXPR_CACHE_DIRTY = True
        _save_gtex_cache()
        return None

    result: Optional[dict] = None
    try:
        r = requests.get(
            GTEX_EXPRESSION,
            params={
                "gencodeId": gencode_id,
                "datasetId": "gtex_v8",
            },
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            log.warning("GTEx expression %s (%s): HTTP %d", gene_symbol, gencode_id, r.status_code)
        else:
            data = r.json()
            rows = data.get("data", [])
            expr: dict[str, float] = {}
            for row in rows:
                tissue = row.get("tissueSiteDetailId") or row.get("tissueSiteDetail")
                tpm_values = row.get("data", [])
                # Compute median from per-sample values (GTEx v2 API returns sample-level data)
                if tissue and tpm_values:
                    try:
                        median_tpm = statistics.median(tpm_values)
                        expr[tissue] = round(float(median_tpm), 3)
                    except (TypeError, ValueError, statistics.StatisticsError):
                        pass
                # Fallback: direct 'median' field (some API versions)
                elif tissue:
                    tpm = row.get("median")
                    if tpm is not None:
                        try:
                            expr[tissue] = float(tpm)
                        except ValueError:
                            pass
            result = expr if expr else None

    except Exception as exc:
        log.warning("GTEx expression %s: %s", gene_symbol, exc)

    _GTEX_EXPR_CACHE[key] = result
    _GTEX_EXPR_CACHE_DIRTY = True
    _save_gtex_cache()
    return result


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

        symbol = _hgnc_symbol(gene.gene_symbol)
        log.info("GTEx querying: %s (raw: %s)", symbol, gene.gene_symbol)
        expr = fetch_gtex_expression(symbol)
        time.sleep(0.2)

        if not expr:
            log.info("GTEx %s: no expression data found.", symbol)
            continue

        expressed_tissues = [t for t, tpm in expr.items() if tpm >= TPM_THRESHOLD]
        max_tpm = max(expr.values()) if expr else 0.0
        log.info("GTEx %s: tissues=%d  max_tpm=%.1f", symbol, len(expressed_tissues), max_tpm)

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
