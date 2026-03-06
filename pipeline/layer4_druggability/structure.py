"""AlphaFold DB — download pre-computed structures for human proteins."""

import logging
from pathlib import Path
from typing import Optional

import requests

from db.models import Gene
from db.session import get_session
from pipeline.config import get_storage_root

log = logging.getLogger(__name__)

ALPHAFOLD_FILES_BASE = "https://alphafold.ebi.ac.uk/files"


def _structures_dir() -> Path:
    root = get_storage_root()
    if root.startswith("s3://"):
        d = Path("/tmp/bioresilient/structures")
    else:
        d = Path(root) / "structures"
    d.mkdir(parents=True, exist_ok=True)
    return d


def download_alphafold_structure(uniprot_id: str) -> Optional[Path]:
    """Download AlphaFold PDB for a UniProt ID. Returns path to local file or None."""
    if not uniprot_id or len(uniprot_id) < 2:
        return None
    out_path = _structures_dir() / f"{uniprot_id}.pdb"
    if out_path.exists() and out_path.stat().st_size > 500:
        return out_path
    # EBI AlphaFold file naming: AF-{uniprot}-F1-model_v4.pdb
    url = f"{ALPHAFOLD_FILES_BASE}/AF-{uniprot_id}-F1-model_v4.pdb"
    try:
        r = requests.get(url, timeout=60, stream=True)
        if r.status_code != 200:
            return None
        out_path.write_bytes(r.content)
        return out_path
    except Exception as exc:
        log.debug("AlphaFold download %s: %s", uniprot_id, exc)
        return None


def ensure_structures_for_genes(gene_ids: list[str]) -> dict[str, Path]:
    """Download AlphaFold structures for genes (using Gene.human_protein as UniProt ID). Returns {gene_id: path}."""
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
