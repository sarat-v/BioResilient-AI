"""AlphaFold DB — download pre-computed structures for human proteins.

Caching strategy (three-tier, fastest first):
  1. Local disk (/tmp/bioresilient/structures/) — same container, same run
  2. S3 bucket  (cache/structures/{uniprot}.pdb) — survives spot interruptions
  3. EBI AlphaFold portal — original source; result cached to S3 immediately

On a 14-gene Tier1 set, skipping EBI downloads saves ~60-90s per retry and
reduces EBI API load. The PDB files are ~2-5 MB each (v4 model), so S3 I/O
is negligible (same-region, effectively free within the Batch task's VPC).
"""

import logging
from pathlib import Path
from typing import Optional

import requests

from db.models import Gene
from db.session import get_session
from pipeline.config import get_local_storage_root, get_storage_root

log = logging.getLogger(__name__)

ALPHAFOLD_FILES_BASE = "https://alphafold.ebi.ac.uk/files"
_S3_STRUCTURES_PREFIX = "cache/structures"


def _structures_dir() -> Path:
    d = Path(get_local_storage_root()) / "structures"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _s3_bucket() -> Optional[str]:
    root = get_storage_root()
    if root.startswith("s3://"):
        return root[len("s3://"):].split("/")[0]
    return None


def _s3_restore(uniprot_id: str, local_path: Path) -> bool:
    bucket = _s3_bucket()
    if not bucket:
        return False
    key = f"{_S3_STRUCTURES_PREFIX}/{uniprot_id}.pdb"
    try:
        import boto3
        boto3.client("s3").download_file(bucket, key, str(local_path))
        log.debug("AlphaFold S3 hit: %s", uniprot_id)
        return True
    except Exception:
        return False


def _s3_store(uniprot_id: str, local_path: Path) -> None:
    bucket = _s3_bucket()
    if not bucket or not local_path.exists():
        return
    key = f"{_S3_STRUCTURES_PREFIX}/{uniprot_id}.pdb"
    try:
        import boto3
        boto3.client("s3").upload_file(str(local_path), bucket, key)
        log.debug("AlphaFold cached to S3: %s", uniprot_id)
    except Exception as exc:
        log.debug("AlphaFold S3 store failed (non-fatal): %s", exc)


def download_alphafold_structure(uniprot_id: str) -> Optional[Path]:
    """Download AlphaFold PDB for a UniProt ID.

    Uses three-tier caching: local disk → S3 → EBI portal.
    Returns path to local file or None on failure.
    """
    if not uniprot_id or len(uniprot_id) < 2:
        return None
    out_path = _structures_dir() / f"{uniprot_id}.pdb"

    # Tier 1: local disk (already in container from this run)
    if out_path.exists() and out_path.stat().st_size > 500:
        return out_path

    # Tier 2: S3 cache (survives spot interruptions)
    if _s3_restore(uniprot_id, out_path):
        if out_path.exists() and out_path.stat().st_size > 500:
            return out_path

    # Tier 3: EBI AlphaFold portal (v4 model)
    url = f"{ALPHAFOLD_FILES_BASE}/AF-{uniprot_id}-F1-model_v4.pdb"
    try:
        r = requests.get(url, timeout=60, stream=True)
        if r.status_code != 200:
            return None
        out_path.write_bytes(r.content)
        # Populate S3 cache for future retries
        _s3_store(uniprot_id, out_path)
        return out_path
    except Exception as exc:
        log.debug("AlphaFold download %s: %s", uniprot_id, exc)
        return None


def ensure_structures_for_genes(gene_ids: list[str]) -> dict[str, Path]:
    """Download AlphaFold structures for genes (using Gene.human_protein as UniProt ID).

    Returns {gene_id: path}.
    """
    with get_session() as session:
        gene_list = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()
    result = {}
    for g in gene_list:
        uniprot = (g.human_protein or "").strip()
        if not uniprot:
            continue
        path = download_alphafold_structure(uniprot)
        if path:
            result[g.id] = path
    log.info("AlphaFold structures: %d / %d genes.", len(result), len(gene_ids))
    return result
