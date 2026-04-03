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

Performance design:
  - LOEUF is fetched once from the gnomAD TSV (downloaded/cached on first run).
  - LOEUF values are bulk-loaded into a temp table via COPY.
  - The full classification is executed as a single SQL UPDATE … CASE statement,
    joined against the LOEUF temp table — no Python row iteration at write time.
  - Optional scatter mode: caller passes gene_ids to process a subset; multiple
    Nextflow tasks run in parallel and each handles a non-overlapping gene batch.
  - SQL COPY + UPDATE on 1.69M rows completes in ~30–60s vs ~3h for executemany.
"""

import csv
import gzip
import io
import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import psycopg2
import requests

from db.models import DivergentMotif, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_db_url, get_local_storage_root, get_tool_config

log = logging.getLogger(__name__)

# gnomAD v4.1 bulk constraint TSV — one download covers all ~20k genes, no API needed.
# Amazon S3 mirror is used on AWS EC2 for zero-cost, high-speed transfer within us-east-1.
_GNOMAD_TSV_URL_S3 = (
    "https://gnomad-public-us-east-1.s3.amazonaws.com"
    "/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
)
_GNOMAD_TSV_URL_GCS = (
    "https://storage.googleapis.com/gcp-public-data--gnomad"
    "/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
)
# GraphQL fallback for genes absent from the TSV (rare, e.g. very new gene symbols)
_GNOMAD_API = "https://gnomad.broadinstitute.org/api"
_GNOMAD_GENE_QUERY = """
{
  gene(gene_symbol: "%s", reference_genome: GRCh38) {
    gnomad_constraint {
      oe_lof_upper
    }
  }
}
"""

# Classification thresholds — named constants so they can be found and changed together.
_AM_PATHOGENIC_THRESHOLD = 0.6     # AlphaMissense: scores above this indicate high pathogenicity in human
_ESM_DESTAB_THRESHOLD = -2.0       # ESM-1v LLR: scores below this indicate protein destabilisation
_LOEUF_INTOLERANT_THRESHOLD = 0.5  # gnomAD LOEUF: scores at or below this indicate LoF intolerance

# ---------------------------------------------------------------------------
# gnomAD constraint TSV — download once, parse into memory, cache to disk
# ---------------------------------------------------------------------------


def _tsv_cache_path() -> Path:
    return Path(get_local_storage_root()) / "gnomad" / "gnomad.v4.1.constraint_metrics.tsv.gz"


def _ensure_gnomad_tsv() -> Path:
    """Download the gnomAD v4.1 constraint TSV if not already cached locally.

    Uses the S3 mirror first (fast, free on EC2), falls back to GCS.
    File is stored compressed to save ~75% disk space.
    """
    dest = _tsv_cache_path()
    if dest.exists():
        log.info("gnomAD constraint TSV already cached at %s", dest)
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".tmp")

    for url in (_GNOMAD_TSV_URL_S3, _GNOMAD_TSV_URL_GCS):
        log.info("Downloading gnomAD constraint TSV from %s ...", url)
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                total_mb = int(r.headers.get("content-length", 0)) / 1_048_576
                downloaded = 0
                with gzip.open(tmp, "wt", encoding="utf-8") as fh:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        fh.write(chunk.decode("utf-8", errors="replace"))
                        downloaded += len(chunk)
                        if total_mb:
                            log.info(
                                "  gnomAD TSV: %.1f / %.1f MB downloaded (%.0f%%).",
                                downloaded / 1_048_576,
                                total_mb,
                                100 * downloaded / (total_mb * 1_048_576),
                            )
            tmp.rename(dest)
            log.info("gnomAD constraint TSV saved to %s", dest)
            return dest
        except Exception as exc:
            log.warning("Failed to download gnomAD TSV from %s: %s", url, exc)
            if tmp.exists():
                tmp.unlink()

    raise RuntimeError(
        "Could not download gnomAD constraint TSV from S3 or GCS. "
        "Check network access or place the file manually at: " + str(dest)
    )


def _load_loeuf_from_tsv() -> dict[str, float]:
    """Parse the gnomAD v4.1 constraint TSV and return {gene_symbol: loeuf}.

    Only canonical transcripts are used (one row per gene).
    Genes without a valid numeric lof.oe_ci.upper are omitted — callers
    treat missing entries as unknown (neutral classification).
    """
    path = _ensure_gnomad_tsv()
    loeuf_map: dict[str, float] = {}

    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # Keep only canonical transcripts to get one row per gene
            if row.get("canonical", "").lower() != "true":
                continue
            gene = row.get("gene", "").strip()
            raw = row.get("lof.oe_ci.upper", "").strip()
            if not gene or not raw or raw in ("NA", ""):
                continue
            try:
                loeuf_map[gene] = float(raw)
            except ValueError:
                continue

    log.info("gnomAD TSV loaded: LOEUF available for %d genes.", len(loeuf_map))
    return loeuf_map


def _fetch_loeuf_api(gene_symbol: str) -> Optional[float]:
    """GraphQL fallback for a single gene not found in the TSV (rare)."""
    try:
        r = requests.post(
            _GNOMAD_API,
            json={"query": _GNOMAD_GENE_QUERY % gene_symbol},
            timeout=15,
        )
        if r.status_code == 200:
            constraint = (
                r.json().get("data", {}).get("gene", {}).get("gnomad_constraint") or {}
            )
            val = constraint.get("oe_lof_upper")
            if val is not None:
                return float(val)
    except Exception as exc:
        log.debug("gnomAD API fallback error for %s: %s", gene_symbol, exc)
    return None


def _resolve_uniprot_to_symbols(accessions: set[str]) -> dict[str, str]:
    """Map UniProt accessions to HGNC gene symbols via the UniProt REST API.

    Gene.gene_symbol is populated with bare UniProt accessions (e.g. 'P04637')
    rather than HGNC symbols (e.g. 'TP53') during Phase 1. This function
    translates them before the gnomAD TSV/API lookup.

    Returns {accession: hgnc_symbol}. Unresolvable accessions are absent.
    """
    _UNI_PAT = re.compile(
        r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$"          # reviewed Swiss-Prot format
        r"|^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$"  # TrEMBL 6-char
        r"|^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z0-9]{2}$",  # 10-char isoform base
        re.IGNORECASE,
    )
    valid = [a for a in accessions if _UNI_PAT.match(a)]
    if not valid:
        log.debug("No UniProt accessions detected among %d gene symbols.", len(accessions))
        return {}

    log.info("Resolving %d UniProt accessions → HGNC symbols via UniProt REST…", len(valid))
    result: dict[str, str] = {}
    _BATCH = 500
    for i in range(0, len(valid), _BATCH):
        batch = valid[i : i + _BATCH]
        try:
            r = requests.get(
                "https://rest.uniprot.org/uniprotkb/accessions",
                params={
                    "accessions": ",".join(batch),
                    "fields": "accession,gene_names",
                    "size": str(_BATCH),
                },
                headers={"Accept": "application/json"},
                timeout=60,
            )
            r.raise_for_status()
            for entry in r.json().get("results", []):
                acc = entry.get("primaryAccession", "")
                genes = entry.get("genes", [])
                if genes:
                    sym = genes[0].get("geneName", {}).get("value", "")
                    if sym:
                        result[acc] = sym
        except Exception as exc:
            log.warning("UniProt symbol resolution batch %d failed: %s", i // _BATCH, exc)

    resolved_pct = 100 * len(result) / len(valid) if valid else 0
    log.info(
        "UniProt → HGNC: %d / %d accessions resolved (%.0f%%).",
        len(result), len(valid), resolved_pct,
    )
    return result


def _fetch_loeuf_bulk(gene_symbols: set[str]) -> dict[str, Optional[float]]:
    """Return {gene_symbol: loeuf_or_None} for every requested symbol.

    Primary strategy: parse the gnomAD v4.1 constraint TSV (one download,
    ~20 MB, covers all ~20k genes, zero per-gene API calls).

    Because Gene.gene_symbol is stored as a UniProt accession during Phase 1
    (e.g. 'P04637' for TP53), we first resolve accessions → HGNC symbols via
    UniProt REST, look those up in the TSV, then fall back to the gnomAD
    GraphQL API for any symbols still missing.

    The returned dict is keyed by the *original* gene_symbol value (accession
    or HGNC symbol) so callers don't need to change.
    """
    # Step 0: resolve UniProt accessions → HGNC symbols
    acc_to_sym: dict[str, str] = _resolve_uniprot_to_symbols(gene_symbols)
    # Invert to find the original key for each resolved HGNC symbol
    sym_to_acc: dict[str, str] = {v: k for k, v in acc_to_sym.items()}

    # Effective lookup symbols: resolved HGNC names + any original symbols that
    # didn't look like accessions (already HGNC or unknown format)
    unresolved_originals = gene_symbols - set(acc_to_sym.keys())
    lookup_symbols: set[str] = set(acc_to_sym.values()) | unresolved_originals

    tsv_map = _load_loeuf_from_tsv()

    # Build result keyed by original gene_symbol
    result: dict[str, Optional[float]] = {}
    missing_hgnc: list[str] = []

    for hgnc_sym in lookup_symbols:
        if hgnc_sym in tsv_map:
            orig_key = sym_to_acc.get(hgnc_sym, hgnc_sym)
            result[orig_key] = tsv_map[hgnc_sym]
        else:
            missing_hgnc.append(hgnc_sym)

    tsv_hit_pct = 100 * len(result) / len(gene_symbols) if gene_symbols else 0
    log.info(
        "gnomAD LOEUF: %d / %d symbols resolved from TSV (%.0f%%). "
        "%d will fall back to API.",
        len(result), len(gene_symbols), tsv_hit_pct, len(missing_hgnc),
    )

    if missing_hgnc:
        workers = int(get_tool_config().get("gnomad_workers", 20))
        log.info("  gnomAD API fallback: fetching %d symbols (%d workers)...",
                 len(missing_hgnc), workers)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_fetch_loeuf_api, gs): gs for gs in missing_hgnc}
            for future in as_completed(futures):
                hgnc_sym = futures[future]
                orig_key = sym_to_acc.get(hgnc_sym, hgnc_sym)
                result[orig_key] = future.result()

    # Ensure every original symbol has an entry (None = unknown)
    for gs in gene_symbols:
        result.setdefault(gs, None)

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

    Classification matrix:
    ┌─────────────────────────────────────────────────────────────────────────────┐
    │ ESM destab (LLR < -2) + AM high (> 0.6) + LOEUF tolerant (> 0.5)          │
    │   → loss_of_function  (gene can tolerate LoF; likely protective LoF)       │
    │ ESM destab + AM high + LOEUF intolerant (≤ 0.5)                            │
    │   → functional_shift  (LoF-intolerant gene changed anyway — adaptive)      │
    │ ESM destab + AM high + LOEUF unknown (None)                                │
    │   → functional_shift  (direction unknown but clearly functional change)    │
    │ AM high only (ESM neutral/absent)                                           │
    │   → gain_of_function  (AM says pathogenic but not structurally destab)     │
    │ All other                                                                   │
    │   → neutral                                                                 │
    └─────────────────────────────────────────────────────────────────────────────┘

    Correction vs previous version: when loeuf is None, we previously returned
    'neutral' even with strong ESM + AM signals. This caused the 90% neutral rate
    in step4d because ~40% of genes lack gnomAD LOEUF entries. We now return
    'functional_shift' — we know the variant is functional, direction uncertain.

    Note: These scores are derived from human-centric ML models. A variant
    flagged as 'functional_shift' may actually represent successful adaptation
    in a resilient species. Treat all labels as hypotheses, not conclusions.
    """
    am_high = (consequence_score is not None and consequence_score > _AM_PATHOGENIC_THRESHOLD)
    esm_destab = (esm1v_score is not None and esm1v_score < _ESM_DESTAB_THRESHOLD)

    if esm_destab and am_high:
        if loeuf is None:
            return "functional_shift"
        if loeuf > _LOEUF_INTOLERANT_THRESHOLD:
            return "loss_of_function"
        else:
            return "functional_shift"
    elif am_high and not esm_destab:
        return "gain_of_function"
    else:
        return "neutral"


# ---------------------------------------------------------------------------
# Fast SQL-based bulk update
# ---------------------------------------------------------------------------

def _bulk_update_via_sql(conn, loeuf_map: dict[str, Optional[float]]) -> int:
    """Classify all motifs and write results in a single SQL UPDATE.

    Strategy:
      1. COPY loeuf values into a temp table _loeuf(gene_symbol TEXT, loeuf DOUBLE).
      2. Run a single UPDATE divergent_motif SET motif_direction = CASE ... END
         joined against the temp table.

    This avoids Python-level row iteration for the write path completely.
    On 1.69M rows this takes ~30–60s on RDS vs ~3h for executemany.

    Returns the number of rows updated.
    """
    with conn.cursor() as cur:
        # Build the LOEUF temp table
        cur.execute("""
            CREATE TEMP TABLE IF NOT EXISTS _loeuf_tmp (
                gene_symbol TEXT PRIMARY KEY,
                loeuf       DOUBLE PRECISION
            ) ON COMMIT DROP
        """)
        buf = io.StringIO()
        for sym, val in loeuf_map.items():
            if val is None:
                buf.write(f"{sym}\t\\N\n")
            else:
                buf.write(f"{sym}\t{val}\n")
        buf.seek(0)
        cur.copy_from(buf, "_loeuf_tmp", columns=("gene_symbol", "loeuf"), null="\\N")
        log.info("LOEUF temp table populated: %d rows.", len(loeuf_map))

        # Single UPDATE using a CTE that LEFT JOINs loeuf onto each motif.
        # The CASE mirrors classify_motif_direction() exactly.
        cur.execute(f"""
            WITH motif_loeuf AS (
                SELECT dm.id AS motif_id,
                       dm.esm1v_score,
                       dm.consequence_score,
                       lt.loeuf
                FROM   divergent_motif dm
                JOIN   ortholog o        ON o.id          = dm.ortholog_id
                JOIN   gene     g        ON g.id          = o.gene_id
                LEFT JOIN _loeuf_tmp lt  ON lt.gene_symbol = g.gene_symbol
            )
            UPDATE divergent_motif AS dm
            SET motif_direction = CASE
                -- ESM destabilising + AM high → branch on LOEUF
                WHEN ml.esm1v_score < {_ESM_DESTAB_THRESHOLD}
                  AND ml.consequence_score > {_AM_PATHOGENIC_THRESHOLD}
                THEN
                    CASE
                        WHEN ml.loeuf IS NULL                        THEN 'functional_shift'
                        WHEN ml.loeuf > {_LOEUF_INTOLERANT_THRESHOLD} THEN 'loss_of_function'
                        ELSE                                              'functional_shift'
                    END
                -- AM high only (ESM neutral or missing) → gain_of_function
                WHEN ml.consequence_score > {_AM_PATHOGENIC_THRESHOLD}
                THEN 'gain_of_function'
                -- Everything else → neutral
                ELSE 'neutral'
            END
            FROM motif_loeuf ml
            WHERE dm.id = ml.motif_id
        """)
        updated = cur.rowcount
        conn.commit()
        log.info("SQL UPDATE complete: %d motif rows classified.", updated)
        return updated


def _get_direction_distribution(conn) -> dict[str, int]:
    with conn.cursor() as cur:
        cur.execute("""
            SELECT motif_direction, COUNT(*)
            FROM divergent_motif
            GROUP BY motif_direction
        """)
        return {row[0] or "null": row[1] for row in cur.fetchall()}


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def annotate_variant_directions(gene_ids: Optional[list[str]] = None) -> int:
    """Classify each DivergentMotif and store motif_direction.

    Uses a SQL-native bulk UPDATE via a LOEUF temp table for maximum throughput.
    When gene_ids is provided, only those genes are updated (scatter mode).

    Args:
        gene_ids: Optional list of gene IDs to limit processing.

    Returns:
        Number of motifs annotated.
    """
    # ------------------------------------------------------------------
    # 1. Determine which gene symbols we need LOEUF for.
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

    log.info("Step 4d: fetching LOEUF for %d gene symbols...", len(gene_symbols_needed))
    loeuf_map: dict[str, Optional[float]] = _fetch_loeuf_bulk(gene_symbols_needed)

    # ------------------------------------------------------------------
    # 2. Run fast SQL-based bulk classification.
    # ------------------------------------------------------------------
    conn = psycopg2.connect(get_db_url())
    try:
        if gene_ids:
            # Scatter mode: restrict the UPDATE to the specified gene subset.
            # Classify rows in Python then COPY results into a temp table and
            # do a targeted UPDATE by motif_id — keeps scatter tasks isolated.
            updated = _scatter_update(conn, loeuf_map, gene_ids)
        else:
            # Full-table fast SQL path.
            updated = _bulk_update_via_sql(conn, loeuf_map)

        dist = _get_direction_distribution(conn)
        log.info("Variant direction distribution: %s", dist)
    finally:
        conn.close()

    return updated


def _scatter_update(
    conn,
    loeuf_map: dict[str, Optional[float]],
    gene_ids: list[str],
) -> int:
    """Targeted UPDATE for a subset of genes (scatter mode).

    Fetches only the motifs for the given gene_ids, classifies them in Python,
    then writes results via COPY + temp table UPDATE. This keeps each scatter
    worker independent and avoids full-table locks.
    """
    with conn.cursor() as cur:
        cur.execute("""
            SELECT dm.id, dm.esm1v_score, dm.consequence_score, g.gene_symbol
            FROM   divergent_motif dm
            JOIN   ortholog o ON o.id = dm.ortholog_id
            JOIN   gene     g ON g.id = o.gene_id
            WHERE  g.id = ANY(%s)
        """, (gene_ids,))
        rows = cur.fetchall()

    log.info("Scatter: classifying %d motifs for %d genes...", len(rows), len(gene_ids))

    updates: list[tuple[str, str]] = []  # (direction, motif_id)
    dist: dict[str, int] = {}
    for motif_id, esm1v_score, consequence_score, gene_symbol in rows:
        loeuf = loeuf_map.get(gene_symbol) if gene_symbol else None
        direction = classify_motif_direction(esm1v_score, consequence_score, loeuf)
        updates.append((direction, str(motif_id)))
        dist[direction] = dist.get(direction, 0) + 1

    if not updates:
        return 0

    # Write via COPY + temp table UPDATE — no executemany
    with conn.cursor() as cur:
        cur.execute("""
            CREATE TEMP TABLE _dir_update (
                motif_id  TEXT,
                direction TEXT
            ) ON COMMIT DROP
        """)
        buf = io.StringIO()
        for direction, mid in updates:
            buf.write(f"{mid}\t{direction}\n")
        buf.seek(0)
        cur.copy_from(buf, "_dir_update", columns=("motif_id", "direction"))

        cur.execute("""
            UPDATE divergent_motif dm
            SET    motif_direction = u.direction
            FROM   _dir_update u
            WHERE  dm.id = u.motif_id
        """)
        n = cur.rowcount
        conn.commit()

    log.info("Scatter update: %d motifs written. Distribution: %s", n, dist)
    return n


def run_variant_direction_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    return annotate_variant_directions(gene_ids)
