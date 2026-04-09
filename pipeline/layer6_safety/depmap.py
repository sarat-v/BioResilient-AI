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

Score interpretation:
  < -0.5: broadly essential — high safety risk
  -0.5 to -0.2: context-dependent essentiality — moderate risk
  > -0.2: not broadly essential — safer target
"""

import csv
import logging
import time
from pathlib import Path
from typing import Optional

import requests

from db.models import Gene, SafetyFlag
from db.session import get_session
from pipeline.config import get_local_storage_root

log = logging.getLogger(__name__)

DEPMAP_SUMMARY_URL = "https://depmap.org/portal/api/download/file?file_name=CRISPRGeneDependency.csv&release=DepMap+Public+24Q2"
REQUEST_TIMEOUT = 30


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


def download_depmap_scores(force: bool = False) -> Path:
    """Download DepMap CRISPR gene dependency scores.

    Lookup order (fastest first):
      1. Local disk (/tmp/bioresilient/depmap/) — present after first run in same container
      2. S3 bucket (cache/depmap/) — present after any previous successful run
      3. DepMap portal (~150 MB) — original source; result is uploaded to S3 for caching
    """
    path = _depmap_cache_path()
    if path.exists() and not force:
        log.info("DepMap scores already cached locally at %s", path)
        return path

    # S3 cache check — avoids 150 MB portal download on Spot retries
    if not force and _try_s3_download(_DEPMAP_S3_KEY, path):
        return path

    log.info("Downloading DepMap CRISPR scores from portal (~150 MB)...")
    try:
        with requests.get(DEPMAP_SUMMARY_URL, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(path, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    f.write(chunk)
        log.info("DepMap download complete: %s", path)
        # Cache to S3 so spot retries skip the portal download
        _try_s3_upload(_DEPMAP_S3_KEY, path)
    except Exception as exc:
        log.warning("DepMap download failed: %s — skipping essentiality scoring.", exc)
        if path.exists():
            path.unlink()
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


def annotate_depmap(
    depmap_index: dict[str, float],
    gene_ids: Optional[list[str]] = None,
) -> int:
    """Annotate SafetyFlag rows with DepMap essentiality scores.

    Returns number of genes annotated.
    """
    with get_session() as session:
        if gene_ids is None:
            from db.models import Gene
            gene_ids = [g.id for g in session.query(Gene).all()]
        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()

    annotated = 0
    with get_session() as session:
        for gene in genes:
            symbol = gene.gene_symbol
            if not symbol or symbol not in depmap_index:
                continue

            score = depmap_index[symbol]
            sf = session.get(SafetyFlag, gene.id)
            if sf is None:
                sf = SafetyFlag(gene_id=gene.id)
                session.add(sf)

            sf.depmap_score = score
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
