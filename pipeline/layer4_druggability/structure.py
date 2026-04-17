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
ALPHAFOLD_API_BASE  = "https://alphafold.ebi.ac.uk/api"
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


def _oma_to_entry_name(human_protein: str) -> str:
    """Strip OMA species prefix from human_protein field.

    OMA stores IDs as 'human|AKT1_HUMAN'. AlphaFold API expects 'AKT1_HUMAN'.
    """
    return human_protein.split("|", 1)[-1] if "|" in human_protein else human_protein


def _entry_name_to_accession(entry_name: str) -> Optional[str]:
    """Resolve a UniProt entry name (e.g. AKT1_HUMAN) to a UniProt accession (e.g. P31749).

    AlphaFold API requires accession, not entry name.
    """
    try:
        r = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={"query": f"id:{entry_name}", "fields": "accession", "format": "json"},
            timeout=15,
        )
        if r.status_code == 200:
            results = r.json().get("results", [])
            if results:
                return results[0].get("primaryAccession")
    except Exception as exc:
        log.debug("UniProt accession lookup %s: %s", entry_name, exc)
    return None


def _accession_to_pdb_url(accession: str) -> Optional[str]:
    """Get the latest AlphaFold PDB URL for a UniProt accession via the AlphaFold API."""
    try:
        r = requests.get(
            f"{ALPHAFOLD_API_BASE}/prediction/{accession}",
            timeout=15,
        )
        if r.status_code == 200:
            data = r.json()
            if data and data[0].get("pdbUrl"):
                return data[0]["pdbUrl"]
    except Exception as exc:
        log.debug("AlphaFold API lookup %s: %s", accession, exc)
    return None


def download_alphafold_structure(human_protein: str) -> Optional[Path]:
    """Download AlphaFold PDB for a gene's human_protein identifier.

    Accepts both OMA-style ('human|AKT1_HUMAN') and bare entry names ('AKT1_HUMAN').
    Uses three-tier caching: local disk → S3 → EBI portal.
    Returns path to local PDB file or None on failure.
    """
    if not human_protein:
        return None

    entry_name = _oma_to_entry_name(human_protein.strip())
    if not entry_name:
        return None

    out_path = _structures_dir() / f"{entry_name}.pdb"

    # Tier 1: local disk
    if out_path.exists() and out_path.stat().st_size > 500:
        return out_path

    # Tier 3: EBI AlphaFold portal — resolve entry name → accession → PDB URL
    accession = _entry_name_to_accession(entry_name)
    if not accession:
        log.debug("No UniProt accession found for %s", entry_name)
        return None

    # Also check S3 cache by accession
    acc_path = _structures_dir() / f"{accession}.pdb"
    if _s3_restore(accession, acc_path):
        if acc_path.exists() and acc_path.stat().st_size > 500:
            return acc_path

    pdb_url = _accession_to_pdb_url(accession)
    if not pdb_url:
        log.debug("No AlphaFold structure found for %s (%s)", entry_name, accession)
        return None
    try:
        r = requests.get(pdb_url, timeout=60, stream=True)
        if r.status_code != 200:
            log.debug("AlphaFold download %s (%s): HTTP %d", entry_name, accession, r.status_code)
            return None
        acc_path.write_bytes(r.content)
        _s3_store(accession, acc_path)
        return acc_path
    except Exception as exc:
        log.debug("AlphaFold download %s: %s", entry_name, exc)
        return None


def ensure_structures_for_genes(gene_ids: list[str]) -> dict[str, Path]:
    """Download AlphaFold structures for genes using Gene.human_protein.

    Returns {gene_id: local_pdb_path}.
    """
    with get_session() as session:
        gene_list = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()
    result = {}
    for g in gene_list:
        if not g.human_protein:
            continue
        path = download_alphafold_structure(g.human_protein)
        if path:
            result[g.id] = path
    log.info("AlphaFold structures: %d / %d genes.", len(result), len(gene_ids))
    return result
