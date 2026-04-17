"""Step 14b-i — DepMap CRISPR essentiality scoring for drug target safety.

DepMap (Dependency Map) provides genome-wide CRISPR knockout fitness data
across hundreds of cancer cell lines. The Chronos score measures gene
essentiality: values near -1 mean the gene is broadly essential (knocking
it out kills most cell lines); values near 0 mean it's non-essential.

For drug target safety, we want genes that are:
  - NOT broadly essential across many tissue types (avoids toxicity)
  - May be essential in tumour cells but not normal cells (window for therapy)

We use DepMap's public summary file which provides the mean Chronos score
per gene across all screened cell lines.

Data source:
  DepMap Public Release (https://depmap.org/portal/download/all/)
  File: CRISPRGeneDependency.csv (gene-level mean scores)
  Cached at {storage_root}/depmap/CRISPRGeneDependency.csv

Score interpretation (probability of dependency, 0–1):
  > 0.7: very broadly essential — high safety risk
  0.5–0.7: broadly essential — moderate safety concern
  < 0.3: not broadly essential — safer target
"""

import csv
import io
import logging
import time
from pathlib import Path
from typing import Optional

import requests

from db.models import Gene, SafetyFlag
from db.session import get_session
from pipeline.config import get_local_storage_root

log = logging.getLogger(__name__)

DEPMAP_FILES_API = "https://depmap.org/portal/api/download/files"
DEPMAP_TARGET_FILE = "CRISPRGeneDependency.csv"
REQUEST_TIMEOUT = 60
DEPMAP_FALLBACK_URLS = [
    # Known working direct file URLs
    "https://storage.googleapis.com/depmap-external-downloads/downloads-by-canonical-id/26q1-public-3b44.1/CRISPRGeneDependency.csv",
    "https://ndownloader.figshare.com/files/51064631",  # 24Q4
    "https://ndownloader.figshare.com/files/46489021",  # 24Q2
]


def _depmap_cache_path() -> Path:
    root = Path(get_local_storage_root()) / "depmap"
    root.mkdir(parents=True, exist_ok=True)
    return root / "CRISPRGeneDependency.csv"


_DEPMAP_S3_KEY = "cache/depmap/CRISPRGeneDependency.csv"


def _try_s3_download(s3_key: str, local_path: Path) -> bool:
    """Try to restore *s3_key* from the pipeline S3 bucket to *local_path*.

    Returns True on success, False if the object doesn't exist or on error.
    Same-region downloads (ap-south-1 → ap-south-1) are effectively free
    and run at ~100 MB/s, making this much faster than hitting the DepMap
    portal every retry.
    """
    from pipeline.config import get_storage_root
    storage_root = get_storage_root()
    if not storage_root.startswith("s3://"):
        return False
    bucket = storage_root[len("s3://"):].split("/")[0]
    try:
        import boto3
        boto3.client("s3").download_file(bucket, s3_key, str(local_path))
        log.info("DepMap: loaded from S3 cache s3://%s/%s", bucket, s3_key)
        return True
    except Exception:
        return False


def _try_s3_upload(s3_key: str, local_path: Path) -> None:
    """Upload *local_path* to the pipeline S3 bucket for future retries."""
    from pipeline.config import get_storage_root
    storage_root = get_storage_root()
    if not storage_root.startswith("s3://"):
        return
    bucket = storage_root[len("s3://"):].split("/")[0]
    try:
        import boto3
        boto3.client("s3").upload_file(str(local_path), bucket, s3_key)
        log.info("DepMap: cached to s3://%s/%s for future retries.", bucket, s3_key)
    except Exception as exc:
        log.debug("DepMap S3 cache upload failed (non-fatal): %s", exc)


def _resolve_depmap_url() -> Optional[str]:
    """Resolve the current download URL for CRISPRGeneDependency.csv from DepMap.

    DepMap's /api/download/files returns a CSV manifest with columns:
      release, release_date, filename, url, md5_hash
    We prefer the latest release with a GCS URL; fall back to Figshare.
    """
    try:
        r = requests.get(
            DEPMAP_FILES_API,
            timeout=30,
            headers={"User-Agent": "BioResilient/1.0"},
        )
        if r.status_code != 200:
            log.warning("DepMap files API: HTTP %d", r.status_code)
            return DEPMAP_FALLBACK_URLS[0]
        reader = csv.DictReader(io.StringIO(r.text))
        candidates = []
        for row in reader:
            fname = row.get("filename", "") or row.get("name", "")
            url = row.get("url", "") or row.get("downloadUrl", "")
            release = row.get("release", "")
            if fname == DEPMAP_TARGET_FILE and url:
                candidates.append((release, url))
        if not candidates:
            log.warning("DepMap files API: %s not found in manifest.", DEPMAP_TARGET_FILE)
            return DEPMAP_FALLBACK_URLS[0]
        # Sort by release name descending (latest first) — prefer GCS over Figshare
        candidates.sort(key=lambda x: x[0], reverse=True)
        gcs = [(r, u) for r, u in candidates if "storage.googleapis.com" in u]
        chosen = (gcs or candidates)[0]
        log.info("DepMap: resolved %s → %s (%s)", DEPMAP_TARGET_FILE, chosen[0], chosen[1][:80])
        return chosen[1]
    except Exception as exc:
        log.warning("DepMap URL resolution failed: %s", exc)
        return DEPMAP_FALLBACK_URLS[0]


def download_depmap_scores(force: bool = False) -> Path:
    """Download DepMap CRISPR gene dependency scores.

    Lookup order (fastest first):
      1. Local disk (/tmp/bioresilient/depmap/) — present after first run in same container
      2. S3 bucket (cache/depmap/) — present after any previous successful run
      3. DepMap files API (dynamic URL from CSV manifest) → GCS or Figshare
    """
    path = _depmap_cache_path()
    if path.exists() and not force:
        log.info("DepMap scores already cached locally at %s", path)
        return path

    # S3 cache check — avoids large download on Spot retries
    if not force and _try_s3_download(_DEPMAP_S3_KEY, path):
        return path

    url = _resolve_depmap_url()
    if not url:
        log.warning("DepMap: cannot resolve download URL — skipping essentiality scoring.")
        return path

    log.info("Downloading DepMap CRISPR scores from %s ...", url[:80])
    tried = [url] + [u for u in DEPMAP_FALLBACK_URLS if u != url]
    for candidate in tried:
        try:
            with requests.get(
                candidate,
                stream=True,
                timeout=REQUEST_TIMEOUT,
                headers={"User-Agent": "BioResilient/1.0"},
            ) as r:
                r.raise_for_status()
                with open(path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        f.write(chunk)
            log.info("DepMap download complete: %s", path)
            _try_s3_upload(_DEPMAP_S3_KEY, path)
            return path
        except Exception as exc:
            log.warning("DepMap download failed from %s: %s", candidate[:80], exc)
            if path.exists():
                path.unlink()
            continue
    log.warning("DepMap: all download sources failed — skipping essentiality scoring.")
    return path


def load_depmap_index(csv_path: Path) -> dict[str, float]:
    """Parse DepMap CSV and build {gene_symbol: mean_chronos_score} index.

    The CSV has gene symbols in columns (format: "SYMBOL (ENTREZ_ID)").
    The rows are cell lines. We compute the mean across all cell lines per gene.
    """
    if not csv_path.exists():
        return {}

    log.info("Building DepMap index from %s...", csv_path)
    gene_scores: dict[str, list[float]] = {}

    try:
        with open(csv_path, newline="", errors="replace") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header is None:
                return {}

            # Parse gene symbols from header
            gene_symbols = []
            for h in header[1:]:  # skip cell line ID column
                # Format: "SYMBOL (ENTREZ_ID)" or just "SYMBOL"
                symbol = h.split(" (")[0].strip()
                gene_symbols.append(symbol)
                gene_scores[symbol] = []

            for row in reader:
                for i, val in enumerate(row[1:]):
                    if i >= len(gene_symbols):
                        break
                    symbol = gene_symbols[i]
                    try:
                        gene_scores[symbol].append(float(val))
                    except ValueError:
                        pass

    except Exception as exc:
        log.warning("DepMap index build failed: %s", exc)
        return {}

    # Compute mean per gene
    result = {}
    for symbol, scores in gene_scores.items():
        if scores:
            result[symbol] = round(sum(scores) / len(scores), 4)

    log.info("DepMap index built: %d genes.", len(result))
    return result


def _hgnc_symbol(raw: str) -> str:
    """Convert OMA entry name (AKT1_HUMAN) to bare HGNC symbol (AKT1)."""
    if raw and "_" in raw:
        return raw.split("_")[0]
    return raw or ""


def annotate_depmap(
    depmap_index: dict[str, float],
    gene_ids: Optional[list[str]] = None,
) -> int:
    """Annotate SafetyFlag rows with DepMap essentiality scores.

    Returns number of genes annotated.
    """
    with get_session() as session:
        if gene_ids is None:
            gene_ids = [g.id for g in session.query(Gene).all()]
        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()

    annotated = 0
    with get_session() as session:
        for gene in genes:
            symbol = _hgnc_symbol(gene.gene_symbol)
            if not symbol:
                continue
            if symbol not in depmap_index:
                log.info("DepMap: no score for %s (raw: %s)", symbol, gene.gene_symbol)
                continue

            score = depmap_index[symbol]
            sf = session.get(SafetyFlag, gene.id)
            if sf is None:
                sf = SafetyFlag(gene_id=gene.id)
                session.add(sf)

            sf.depmap_score = score
            log.info("DepMap %s: chronos=%.4f", symbol, score)
            annotated += 1

        session.commit()

    log.info("DepMap annotation complete: %d genes annotated.", annotated)
    return annotated


def run_depmap_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step14b."""
    path = download_depmap_scores()
    if not path.exists():
        return 0
    index = load_depmap_index(path)
    if not index:
        return 0
    return annotate_depmap(index, gene_ids)
