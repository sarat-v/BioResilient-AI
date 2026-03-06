"""Step 2 — Download proteomes from NCBI for all registered species.

For each species in species_registry.json:
  1. Query NCBI datasets API using taxid to locate the reference proteome FASTA.
  2. Download to ./data/proteomes/{species_id}.faa (local) or s3://…/proteomes/ (cloud).
  3. Validate: protein count, no truncated FASTA records.

Also downloads the human UniProt canonical proteome (~20k proteins) used as baseline.
"""

import logging
import os
import re
import time
from pathlib import Path
from typing import Optional
from urllib.parse import urlencode

import requests
from Bio import SeqIO

from pipeline.config import (
    get_ncbi_api_key,
    get_ncbi_email,
    get_storage_root,
    get_deployment,
)

log = logging.getLogger(__name__)

NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_DATASETS_DOWNLOAD = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
UNIPROT_HUMAN_FASTA = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed%3Atrue)%20AND%20(organism_id%3A9606)"

_RATE_LIMIT_SLEEP = 0.11   # 10 req/sec with API key; 0.34 without


def _ncbi_params(extra: dict) -> dict:
    params = {
        "api_key": get_ncbi_api_key(),
        "email": get_ncbi_email(),
        "retmode": "text",
    }
    params.update(extra)
    return {k: v for k, v in params.items() if v}


def _proteomes_dir() -> Path:
    root = get_storage_root()
    if root.startswith("s3://"):
        # For cloud, still cache locally under /tmp
        d = Path("/tmp/bioresillient/proteomes")
    else:
        d = Path(root) / "proteomes"
    d.mkdir(parents=True, exist_ok=True)
    return d


def download_proteome_ncbi(species_id: str, taxid: int, assembly: str) -> Optional[Path]:
    """Download reference proteome FASTA for a given assembly accession.

    Uses NCBI Datasets v2 API to fetch protein FASTA for the assembly.
    Falls back to Entrez protein search if the Datasets endpoint fails.
    """
    out_dir = _proteomes_dir()
    out_path = out_dir / f"{species_id}.faa"

    if out_path.exists() and out_path.stat().st_size > 10_000:
        log.info("  %s: already downloaded (%s bytes)", species_id, out_path.stat().st_size)
        return out_path

    log.info("  Downloading proteome for %s (taxid=%s, assembly=%s)...", species_id, taxid, assembly)

    # Primary: NCBI Datasets v2 — protein FASTA package
    try:
        url = NCBI_DATASETS_DOWNLOAD.format(accession=assembly)
        params = {
            "include_annotation_type": "PROT_FASTA",
            "api-key": get_ncbi_api_key(),
        }
        headers = {"Accept": "application/zip"}
        r = requests.get(url, params=params, headers=headers, stream=True, timeout=300)
        if r.status_code == 200:
            import io
            import zipfile

            zip_bytes = io.BytesIO(r.content)
            with zipfile.ZipFile(zip_bytes) as zf:
                faa_files = [n for n in zf.namelist() if n.endswith(".faa") or n.endswith(".faa.gz")]
                if faa_files:
                    with zf.open(faa_files[0]) as src, open(out_path, "wb") as dst:
                        if faa_files[0].endswith(".gz"):
                            import gzip
                            with gzip.open(src) as gz:
                                dst.write(gz.read())
                        else:
                            dst.write(src.read())
                    count = _count_sequences(out_path)
                    if count > 0:
                        log.info("    ✓ %s: %d proteins (Datasets API)", species_id, count)
                        return out_path
    except Exception as exc:
        log.warning("    Datasets API failed for %s: %s — trying Entrez fallback", species_id, exc)

    # Fallback: Entrez protein search by taxid
    time.sleep(_RATE_LIMIT_SLEEP)
    out_path = _download_via_entrez(species_id, taxid, out_path)
    return out_path


def _download_via_entrez(species_id: str, taxid: int, out_path: Path) -> Optional[Path]:
    """Entrez protein database search + fetch as FASTA fallback."""
    try:
        search_params = _ncbi_params({
            "db": "protein",
            "term": f"txid{taxid}[Organism:exp] AND refseq[Filter]",
            "retmax": 50000,
            "usehistory": "y",
        })
        r = requests.get(NCBI_ESEARCH, params=search_params, timeout=60)
        r.raise_for_status()
        time.sleep(_RATE_LIMIT_SLEEP)

        webenv = re.search(r"<WebEnv>(\S+)</WebEnv>", r.text)
        query_key = re.search(r"<QueryKey>(\d+)</QueryKey>", r.text)
        count = re.search(r"<Count>(\d+)</Count>", r.text)

        if not webenv or not query_key:
            log.error("    Entrez search returned no results for %s", species_id)
            return None

        n_proteins = int(count.group(1)) if count else 0
        log.info("    Entrez: %d proteins found for taxid %s", n_proteins, taxid)

        fetch_params = _ncbi_params({
            "db": "protein",
            "query_key": query_key.group(1),
            "WebEnv": webenv.group(1),
            "rettype": "fasta",
            "retmode": "text",
            "retmax": min(n_proteins, 50000),
        })
        r2 = requests.get(NCBI_EFETCH, params=fetch_params, timeout=600, stream=True)
        r2.raise_for_status()

        with open(out_path, "wb") as f:
            for chunk in r2.iter_content(chunk_size=65536):
                f.write(chunk)

        count = _count_sequences(out_path)
        log.info("    ✓ %s: %d proteins (Entrez)", species_id, count)
        return out_path

    except Exception as exc:
        log.error("    Entrez download failed for %s: %s", species_id, exc)
        return None


def download_human_proteome() -> Path:
    """Download the human UniProt Swiss-Prot reviewed proteome."""
    out_dir = _proteomes_dir()
    out_path = out_dir / "human.faa"

    if out_path.exists() and out_path.stat().st_size > 5_000_000:
        log.info("  human: already downloaded (%s bytes)", out_path.stat().st_size)
        return out_path

    log.info("  Downloading human UniProt reviewed proteome...")
    try:
        r = requests.get(UNIPROT_HUMAN_FASTA, timeout=600, stream=True)
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=65536):
                f.write(chunk)
        count = _count_sequences(out_path)
        log.info("  ✓ human: %d proteins (UniProt SwissProt)", count)
        return out_path
    except Exception as exc:
        log.error("  Failed to download human proteome: %s", exc)
        raise


def _count_sequences(path: Path) -> int:
    """Count FASTA records in a file. Returns 0 on parse error."""
    try:
        return sum(1 for _ in SeqIO.parse(str(path), "fasta"))
    except Exception:
        return 0


def reheader_fasta(path: Path, species_id: str, out_path: Optional[Path] = None) -> Path:
    """Rewrite FASTA headers to format: >{species_id}|{protein_id}

    OrthoFinder requires clean, unique headers. This ensures compatibility.
    """
    if out_path is None:
        out_path = path.with_suffix(".reheadered.faa")

    records = []
    for rec in SeqIO.parse(str(path), "fasta"):
        # Take only the first word of the current ID
        clean_id = rec.id.split()[0].split("|")[-1]
        rec.id = f"{species_id}|{clean_id}"
        rec.description = ""
        records.append(rec)

    SeqIO.write(records, str(out_path), "fasta")
    log.info("  Reheadered %d sequences for %s → %s", len(records), species_id, out_path)
    return out_path


def validate_proteome(path: Path, min_proteins: int = 100) -> bool:
    """Basic sanity check: file exists, non-empty, minimum protein count."""
    if not path.exists():
        return False
    count = _count_sequences(path)
    if count < min_proteins:
        log.warning("  Validation failed: %s has only %d proteins (min=%d)", path, count, min_proteins)
        return False
    return True


def run_downloads(species_list: list[dict]) -> dict[str, Path]:
    """Download proteomes for all species. Returns map of species_id → FASTA path."""
    paths: dict[str, Path] = {}

    for entry in species_list:
        sid = entry["id"]
        if sid == "human":
            # Human comes from UniProt
            continue
        try:
            raw = download_proteome_ncbi(sid, entry["taxid"], entry.get("genome_assembly", ""))
            if raw and validate_proteome(raw):
                reheadered = reheader_fasta(raw, sid)
                paths[sid] = reheadered
            else:
                log.error("  FAILED: %s proteome download or validation failed", sid)
        except Exception as exc:
            log.error("  FAILED: %s — %s", sid, exc)
        time.sleep(_RATE_LIMIT_SLEEP)

    # Human
    try:
        human_raw = download_human_proteome()
        if validate_proteome(human_raw, min_proteins=10000):
            human_reheadered = reheader_fasta(human_raw, "human")
            paths["human"] = human_reheadered
    except Exception as exc:
        log.error("  FAILED: human proteome — %s", exc)

    log.info("Download complete: %d / %d species succeeded.", len(paths), len(species_list))
    return paths
