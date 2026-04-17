"""One-time backfill: populate gene.loeuf from gnomAD v4.1 constraint TSV.

This script downloads gnomAD v4.1 LOEUF values and writes them directly into
the gene table. Once done, step4d reads loeuf from the DB instead of
downloading gnomAD at runtime inside a throwaway Batch container (which was
silently failing and causing loss_of_function = 0 for every gene).

Run once before the next step4d rerun:

    python -m scripts.backfill_loeuf

Requirements:
    - DATABASE_URL env var (or .env file)
    - Network access to gnomAD S3 (us-east-1) or GCS fallback
    - ~500 MB disk space for the TSV download (cached, gzip-compressed)

After running, verify:
    SELECT COUNT(*) FROM gene WHERE loeuf IS NOT NULL;
    -- Expected: ~11,000-13,000 out of 12,795 genes
    -- (~85-90% of human gene symbols resolve to gnomAD entries)
"""

import csv
import gzip
import io
import logging
import os
import sys
from pathlib import Path

import psycopg2
import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

_GNOMAD_TSV_URL_S3 = (
    "https://gnomad-public-us-east-1.s3.amazonaws.com"
    "/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
)
_GNOMAD_TSV_URL_GCS = (
    "https://storage.googleapis.com/gcp-public-data--gnomad"
    "/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
)

_CACHE_PATH = Path("/tmp/bioresilient/gnomad/gnomad.v4.1.constraint_metrics.tsv.gz")


def _get_db_url() -> str:
    url = os.environ.get("DATABASE_URL", "")
    url = url.replace("${RDS_PASSWORD}", os.environ.get("RDS_PASSWORD", ""))
    url = url.replace("${RDS_HOST}", os.environ.get("RDS_HOST", ""))
    if not url:
        raise RuntimeError("DATABASE_URL not set")
    return url


def _download_gnomad_tsv() -> Path:
    if _CACHE_PATH.exists():
        log.info("gnomAD TSV already cached at %s", _CACHE_PATH)
        return _CACHE_PATH

    _CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    tmp = _CACHE_PATH.with_suffix(".tmp")

    for url in (_GNOMAD_TSV_URL_S3, _GNOMAD_TSV_URL_GCS):
        log.info("Downloading gnomAD v4.1 constraint TSV from %s ...", url)
        try:
            with requests.get(url, stream=True, timeout=300) as r:
                r.raise_for_status()
                total_mb = int(r.headers.get("content-length", 0)) / 1_048_576
                downloaded = 0
                with gzip.open(tmp, "wt", encoding="utf-8") as fh:
                    for chunk in r.iter_content(chunk_size=4 << 20):
                        fh.write(chunk.decode("utf-8", errors="replace"))
                        downloaded += len(chunk)
                        if total_mb:
                            pct = 100 * downloaded / (total_mb * 1_048_576)
                            if pct % 10 < (100 * 4 << 20) / (total_mb * 1_048_576):
                                log.info("  %.0f%% (%.0f / %.0f MB)", pct, downloaded / 1_048_576, total_mb)
            tmp.rename(_CACHE_PATH)
            log.info("Saved to %s", _CACHE_PATH)
            return _CACHE_PATH
        except Exception as exc:
            log.warning("Download from %s failed: %s", url, exc)
            if tmp.exists():
                tmp.unlink()

    raise RuntimeError(
        "Could not download gnomAD TSV from S3 or GCS. "
        "Check network access or place the file manually at: " + str(_CACHE_PATH)
    )


def _load_loeuf_from_tsv(path: Path) -> dict[str, float]:
    """Parse gnomAD TSV → {gene_symbol: loeuf}. Canonical transcripts only."""
    loeuf_map: dict[str, float] = {}
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
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
    log.info("gnomAD TSV parsed: LOEUF available for %d gene symbols.", len(loeuf_map))
    return loeuf_map


def _strip_human_suffix(symbol: str) -> str:
    """Convert 'TP53_HUMAN' → 'TP53' for gnomAD lookup."""
    if symbol.endswith("_HUMAN"):
        return symbol[:-6]
    return symbol


def backfill(dry_run: bool = False) -> None:
    path = _download_gnomad_tsv()
    gnomad_map = _load_loeuf_from_tsv(path)

    conn = psycopg2.connect(_get_db_url())
    try:
        with conn.cursor() as cur:
            cur.execute("SELECT id, gene_symbol FROM gene WHERE gene_symbol IS NOT NULL")
            genes = cur.fetchall()  # [(gene_id, gene_symbol), ...]
        log.info("Gene table has %d rows to update.", len(genes))

        updates = []
        missing = []
        for gene_id, gene_symbol in genes:
            hgnc = _strip_human_suffix(gene_symbol)
            loeuf = gnomad_map.get(hgnc)
            if loeuf is not None:
                updates.append((loeuf, gene_id))
            else:
                missing.append(hgnc)

        log.info(
            "LOEUF resolved for %d / %d genes (%.0f%%). Missing: %d.",
            len(updates),
            len(genes),
            100 * len(updates) / max(len(genes), 1),
            len(missing),
        )
        if missing[:10]:
            log.info("Sample missing symbols: %s", missing[:10])

        if dry_run:
            log.info("DRY RUN — no changes written.")
            return

        # Bulk update via a COPY + temp table strategy (fast on RDS t4g.micro)
        buf = io.StringIO()
        for loeuf_val, gene_id in updates:
            buf.write(f"{gene_id}\t{loeuf_val}\n")
        buf.seek(0)

        with conn.cursor() as cur:
            cur.execute("""
                CREATE TEMP TABLE _loeuf_bulk (
                    gene_id TEXT,
                    loeuf   FLOAT
                ) ON COMMIT DROP
            """)
            cur.copy_from(buf, "_loeuf_bulk", columns=("gene_id", "loeuf"))
            cur.execute("""
                UPDATE gene g
                SET loeuf = lb.loeuf
                FROM _loeuf_bulk lb
                WHERE g.id = lb.gene_id
            """)
            updated_rows = cur.rowcount
        conn.commit()
        log.info("Updated gene.loeuf for %d rows.", updated_rows)

        # Verify
        with conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) FROM gene WHERE loeuf IS NOT NULL")
            count = cur.fetchone()[0]
        log.info("Verification: gene.loeuf populated for %d genes.", count)

    finally:
        conn.close()


if __name__ == "__main__":
    dry_run = "--dry-run" in sys.argv
    backfill(dry_run=dry_run)
    log.info("Done. Now rerun step4d to fix loss_of_function classification.")
