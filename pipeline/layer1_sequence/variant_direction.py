"""Variant direction inference: GoF / LoF / functional_shift / neutral per divergent motif.

Combines three existing signals that are already computed and stored on DivergentMotif:
  - AlphaMissense consequence_score  (0 = benign, 1 = pathogenic in human)
  - ESM-1v LLR score                  (negative = destabilising, positive = neutral/gain)
  - LOEUF from gnomAD                 (fetched once per human gene; low = LoF intolerant)

Classification logic:
  ┌──────────────────────────────────────────────────────────────────────────────────┐
  │ destabilising (ESM-1v LLR < -2)  + high AM score (> 0.6)                       │
  │   + LoF-tolerant LOEUF > 0.5     → loss_of_function (protective LoF)           │
  │ destabilising                    + LoF-intolerant LOEUF ≤ 0.5                  │
  │                                   → functional_shift (hypothesis only — human   │
  │                                     ML tools flag this as pathogenic but the    │
  │                                     variant is fixed in a resilient species;    │
  │                                     may represent adaptive functional change)   │
  │ neutral/stabilising ESM-1v        + high AM  → gain_of_function                 │
  │ all other                          → neutral                                     │
  └──────────────────────────────────────────────────────────────────────────────────┘

IMPORTANT: AlphaMissense and ESM-1v were trained on human clinical variants.
They assume the human reference is optimal and flag deviations as harmful.
In cross-species evolutionary analysis the opposite may be true — a variant fixed
in a long-lived, cancer-resistant species is likely adaptive, not pathogenic.
All direction labels are *hypotheses for experimental validation*, not definitive
mechanism calls. The 'functional_shift' label in particular should be read as
"function likely altered; direction unclear given human-centric model training".

Result stored as DivergentMotif.motif_direction (string enum: 'gain_of_function',
'loss_of_function', 'functional_shift', 'neutral').
"""

import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

import requests

from db.models import DivergentMotif, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_tool_config

log = logging.getLogger(__name__)

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_GENE_QUERY = """
{
  gene(gene_symbol: "%s", reference_genome: GRCh38) {
    gnomad_constraint {
      pLI
      oe_lof_upper
    }
  }
}
"""

# Classification thresholds — named constants so they can be found and changed together.
_AM_PATHOGENIC_THRESHOLD = 0.6   # AlphaMissense: scores above this indicate high pathogenicity in human
_ESM_DESTAB_THRESHOLD = -2.0     # ESM-1v LLR: scores below this indicate protein destabilisation
_LOEUF_INTOLERANT_THRESHOLD = 0.5  # gnomAD LOEUF: scores at or below this indicate LoF intolerance
_GNOMAD_RATE_LIMIT_SLEEP = 0.25  # seconds between gnomAD API requests

# ---------------------------------------------------------------------------
# gnomAD LOEUF fetch (cached in-process for a single run)
# ---------------------------------------------------------------------------

_loeuf_cache: dict[str, Optional[float]] = {}


def _gnomad_workers() -> int:
    """Return number of parallel gnomAD API workers from config (default 20)."""
    return int(get_tool_config().get("gnomad_workers", 20))


def _fetch_loeuf(gene_symbol: str) -> Optional[float]:
    """Fetch LOEUF (LoF o/e upper CI) for a human gene from gnomAD GraphQL API.

    Returns a float in [0, 1+]. Lower = more LoF intolerant.
    Caches results to avoid repeated requests.

    Falls back to None (not pLI) when LOEUF is unavailable because pLI and LOEUF
    are on incompatible scales — converting between them introduces spurious precision.
    """
    if gene_symbol in _loeuf_cache:
        return _loeuf_cache[gene_symbol]

    try:
        r = requests.post(
            GNOMAD_API,
            json={"query": GNOMAD_GENE_QUERY % gene_symbol},
            timeout=15,
        )
        if r.status_code == 200:
            data = r.json()
            gene_data = data.get("data", {}).get("gene", {})
            constraint = gene_data.get("gnomad_constraint") or {}
            loeuf = constraint.get("oe_lof_upper")
            if loeuf is not None:
                _loeuf_cache[gene_symbol] = float(loeuf)
                return _loeuf_cache[gene_symbol]
            # No LOEUF available — do not fall back to pLI since the scales are
            # incompatible. Return None and let the caller treat unknown as neutral.
            log.debug("No LOEUF for %s (pLI available but not used as proxy)", gene_symbol)
        else:
            log.debug("gnomAD API returned %d for gene %s", r.status_code, gene_symbol)
    except Exception as exc:
        log.debug("gnomAD LOEUF fetch error for %s: %s", gene_symbol, exc)

    _loeuf_cache[gene_symbol] = None
    return None


def _fetch_loeuf_bulk(gene_symbols: set[str]) -> dict[str, Optional[float]]:
    """Fetch LOEUF for many gene symbols in parallel using a thread pool.

    Returns {gene_symbol: loeuf_or_None} for every input symbol.
    Already-cached symbols are returned immediately without a network call.
    """
    result: dict[str, Optional[float]] = {}
    to_fetch = [gs for gs in gene_symbols if gs not in _loeuf_cache]
    # Return cached hits immediately
    for gs in gene_symbols:
        if gs in _loeuf_cache:
            result[gs] = _loeuf_cache[gs]

    if not to_fetch:
        return result

    workers = _gnomad_workers()
    log.info("  Fetching gnomAD LOEUF for %d unique gene symbols (%d workers)...",
             len(to_fetch), workers)

    def _worker(gs: str) -> tuple[str, Optional[float]]:
        return gs, _fetch_loeuf(gs)

    done = 0
    total = len(to_fetch)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_worker, gs): gs for gs in to_fetch}
        for future in as_completed(futures):
            gs, loeuf = future.result()
            result[gs] = loeuf
            done += 1
            if done % 500 == 0 or done == total:
                log.info("  gnomAD LOEUF: %d / %d symbols fetched (%.0f%%).",
                         done, total, 100 * done / total)

    return result


# ---------------------------------------------------------------------------
# Direction classifier
# ---------------------------------------------------------------------------

def classify_motif_direction(
    esm1v_score: Optional[float],
    consequence_score: Optional[float],
    loeuf: Optional[float],
) -> str:
    """Return 'gain_of_function', 'loss_of_function', 'functional_shift', or 'neutral'.

    Args:
        esm1v_score: ESM-1v log-likelihood ratio (negative = destabilising)
        consequence_score: AlphaMissense score in [0,1] (high = pathogenic in human)
        loeuf: gnomAD LOEUF value (low = LoF intolerant in human population)

    When loeuf is None (no gnomAD data), constraint is treated as unknown and the
    label defaults to 'neutral' rather than assuming either tolerant or intolerant.

    Note: These scores are derived from human-centric ML models. A variant
    flagged as 'functional_shift' may actually represent successful adaptation
    in a resilient species. Treat all labels as hypotheses, not conclusions.
    """
    am_high = (consequence_score is not None and consequence_score > _AM_PATHOGENIC_THRESHOLD)
    esm_destab = (esm1v_score is not None and esm1v_score < _ESM_DESTAB_THRESHOLD)

    if esm_destab and am_high:
        if loeuf is None:
            # Unknown constraint: cannot confidently assign LoF vs functional_shift
            return "neutral"
        if loeuf > _LOEUF_INTOLERANT_THRESHOLD:
            return "loss_of_function"
        else:
            return "functional_shift"
    elif am_high and not esm_destab:
        return "gain_of_function"
    else:
        return "neutral"


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def annotate_variant_directions(gene_ids: Optional[list[str]] = None) -> int:
    """Classify each DivergentMotif and store motif_direction.

    Args:
        gene_ids: Optional list of gene IDs to limit processing.

    Returns:
        Number of motifs annotated.
    """
    updated = 0
    with get_session() as session:
        q = (
            session.query(DivergentMotif)
            .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
            .join(Gene, Ortholog.gene_id == Gene.id)
        )
        if gene_ids:
            q = q.filter(Gene.id.in_(gene_ids))

        motifs = q.all()
        log.info("Classifying variant direction for %d motifs...", len(motifs))

        # Pre-fetch LOEUF values per gene symbol before the session closes
        gene_symbols_needed: set[str] = set()
        for motif in motifs:
            gs = motif.ortholog.gene.gene_symbol if motif.ortholog and motif.ortholog.gene else None
            if gs:
                gene_symbols_needed.add(gs)

        symbol_map: dict[str, Optional[float]] = _fetch_loeuf_bulk(gene_symbols_needed)

        for motif in motifs:
            gs = (motif.ortholog.gene.gene_symbol
                  if motif.ortholog and motif.ortholog.gene else None)
            loeuf = symbol_map.get(gs) if gs else None
            direction = classify_motif_direction(
                esm1v_score=getattr(motif, "esm1v_score", None),
                consequence_score=getattr(motif, "consequence_score", None),
                loeuf=loeuf,
            )
            motif.motif_direction = direction
            updated += 1

        session.commit()

        # Read distribution from committed values while session is still open
        dist: dict[str, int] = {
            "gain_of_function": 0,
            "loss_of_function": 0,
            "functional_shift": 0,
            "neutral": 0,
        }
        for motif in motifs:
            key = motif.motif_direction or "neutral"
            if key in dist:
                dist[key] += 1
            else:
                dist[key] = 1

    log.info("Variant direction distribution: %s", dist)
    return updated


def run_variant_direction_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    return annotate_variant_directions(gene_ids)
