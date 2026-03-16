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

import psycopg2
import requests

from db.models import DivergentMotif, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_db_url, get_tool_config

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

_CLASSIFY_BATCH_SIZE = 50_000  # rows per DB commit — keeps transaction size manageable


def annotate_variant_directions(gene_ids: Optional[list[str]] = None) -> int:
    """Classify each DivergentMotif and store motif_direction.

    Processes motifs in streaming batches to avoid loading ~1M ORM objects into
    memory at once and to keep individual DB transactions small.

    Args:
        gene_ids: Optional list of gene IDs to limit processing.

    Returns:
        Number of motifs annotated.
    """
    # ------------------------------------------------------------------
    # 1. Determine which gene symbols we need LOEUF for.
    #    Use a lightweight query — just the join keys + gene_symbol.
    # ------------------------------------------------------------------
    with get_session() as session:
        q = (
            session.query(Gene.gene_symbol)
            .join(Ortholog, Ortholog.gene_id == Gene.id)
            .join(DivergentMotif, DivergentMotif.ortholog_id == Ortholog.id)
            .distinct()
        )
        if gene_ids:
            q = q.filter(Gene.id.in_(gene_ids))
        gene_symbols_needed: set[str] = {row[0] for row in q if row[0]}

    symbol_map: dict[str, Optional[float]] = _fetch_loeuf_bulk(gene_symbols_needed)

    # ------------------------------------------------------------------
    # 2. Build a {motif_id → direction} map entirely in Python, loading
    #    only the columns we need (no full ORM hydration of 1M objects).
    # ------------------------------------------------------------------
    log.info("Classifying variant directions (batch size %d)...", _CLASSIFY_BATCH_SIZE)

    conn = psycopg2.connect(get_db_url())
    try:
        with conn.cursor() as cur:
            # Fetch only the columns needed for classification + the gene symbol
            sql = """
                SELECT dm.id,
                       dm.esm1v_score,
                       dm.consequence_score,
                       g.gene_symbol
                FROM   divergent_motif dm
                JOIN   ortholog o  ON o.id = dm.ortholog_id
                JOIN   gene     g  ON g.id = o.gene_id
            """
            params: list = []
            if gene_ids:
                sql += " WHERE g.id = ANY(%s)"
                params.append(gene_ids)

            cur.execute(sql, params or None)
            rows = cur.fetchall()

        log.info("Classifying variant direction for %d motifs...", len(rows))

        dist: dict[str, int] = {
            "gain_of_function": 0,
            "loss_of_function": 0,
            "functional_shift": 0,
            "neutral": 0,
        }

        updates: list[tuple[str, str]] = []  # (direction, motif_id)
        for motif_id, esm1v_score, consequence_score, gene_symbol in rows:
            loeuf = symbol_map.get(gene_symbol) if gene_symbol else None
            direction = classify_motif_direction(
                esm1v_score=esm1v_score,
                consequence_score=consequence_score,
                loeuf=loeuf,
            )
            updates.append((direction, str(motif_id)))
            dist[direction] = dist.get(direction, 0) + 1

        # ------------------------------------------------------------------
        # 3. Write back in batches so no single transaction holds 1M rows.
        # ------------------------------------------------------------------
        total = len(updates)
        committed = 0
        with conn.cursor() as cur:
            for i in range(0, total, _CLASSIFY_BATCH_SIZE):
                batch = updates[i : i + _CLASSIFY_BATCH_SIZE]
                cur.executemany(
                    "UPDATE divergent_motif SET motif_direction = %s WHERE id = %s::uuid",
                    batch,
                )
                conn.commit()
                committed += len(batch)
                log.info("  Variant direction: %d / %d motifs written (%.0f%%).",
                         committed, total, 100 * committed / total)

    finally:
        conn.close()

    log.info("Variant direction distribution: %s", dist)
    return len(updates)


def run_variant_direction_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    return annotate_variant_directions(gene_ids)
