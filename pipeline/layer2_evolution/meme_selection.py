"""Codon-guided MEME selection analysis.

Replaces the protein divergence proxy with statistically rigorous episodic
positive selection detection using HyPhy MEME (Mixed Effects Model of Evolution).

Pipeline for each orthogroup:
  1. Fetch CDS (nucleotide coding sequences) from NCBI for each protein accession
  2. Build a codon alignment guided by the existing protein (amino acid) alignment
     using PAL2NAL-style logic — ensures codons are in the correct reading frame
  3. Run HyPhy MEME on the codon alignment + species tree
  4. Parse per-site selection p-values and report:
       - fraction of sites under episodic positive selection
       - mean effect size (beta+) at significant sites
       - which branches show selection at each site

MEME detects *episodic* selection — sites that were under positive selection
at some point in some lineages, not necessarily in all. This is exactly the
right model for cancer resistance: selection happened once, in specific lineages.

Requirements:
  - HyPhy ≥ 2.5 installed (provides `meme` analysis via hyphy LIBPATH)
  - NCBI API key (for CDS fetching at 10 req/sec)
  - Entrez biopython for CDS lookup
"""

import json
import logging
import os
import pickle
import random
import re
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional
from urllib.parse import urlencode
from urllib.request import urlopen, Request

_rng = random.Random(42)

from Bio import Entrez, SeqIO

from db.models import Species
from db.session import get_session
from pipeline.config import get_ncbi_api_key, get_ncbi_email, get_local_storage_root, get_tool_config

log = logging.getLogger(__name__)

MEME_SIGNIFICANCE_THRESHOLD = 0.05   # per-site p-value threshold for MEME


def export_cds_cache_pkl(aligned_orthogroups: dict[str, dict[str, str]], out_path: str) -> str:
    """Pre-fetch ALL CDS for every accession and export as a single pickle.

    This runs once before the HyPhy scatter so batch tasks never hit NCBI.
    Returns the output path.
    """
    _prefetch_all_cds(aligned_orthogroups)

    cache_dir = Path(get_local_storage_root()) / "cds"
    cds_cache: dict[str, str] = {}
    for fna in cache_dir.glob("*.fna"):
        content = fna.read_text().strip()
        if content:
            cds_cache[fna.stem] = content

    log.info("CDS cache: %d accessions with sequences", len(cds_cache))
    with open(out_path, "wb") as f:
        pickle.dump(cds_cache, f, protocol=pickle.HIGHEST_PROTOCOL)
    log.info("CDS cache written to %s (%.1f MB)", out_path,
             Path(out_path).stat().st_size / 1e6)
    return out_path


_IN_MEMORY_CDS_CACHE: dict[str, str] = {}


def load_cds_cache_pkl(pkl_path: str) -> None:
    """Load a pre-fetched CDS cache pickle into an in-memory dict.

    Each batch task calls this once so fetch_cds_for_protein() can
    look up sequences instantly without disk I/O or NCBI calls.
    """
    global _IN_MEMORY_CDS_CACHE

    with open(pkl_path, "rb") as f:
        _IN_MEMORY_CDS_CACHE = pickle.load(f)

    log.info("CDS cache loaded in memory: %d entries", len(_IN_MEMORY_CDS_CACHE))


# ---------------------------------------------------------------------------
# CDS fetching
# ---------------------------------------------------------------------------

def _setup_entrez() -> None:
    Entrez.email = get_ncbi_email()
    api_key = get_ncbi_api_key()
    if api_key:
        Entrez.api_key = api_key


def _map_uniprot_to_refseq(uniprot_mnemonics: list[str], cache_dir: Path) -> dict[str, str]:
    """Map UniProt entry names (e.g. CASP8_HUMAN) to NCBI RefSeq protein IDs.

    Uses UniProt REST search API with xref_refseq field. Returns {mnemonic: refseq_id}.
    """
    if not uniprot_mnemonics:
        return {}

    mapping: dict[str, str] = {}
    map_cache = cache_dir / "_uniprot_refseq_map.json"
    if map_cache.exists():
        try:
            cached = json.loads(map_cache.read_text())
            mapping.update(cached)
        except Exception:
            pass

    to_map = [m for m in uniprot_mnemonics if m not in mapping]
    if not to_map:
        log.info("UniProt→RefSeq: all %d mappings cached.", len(mapping))
        return {k: v for k, v in mapping.items() if k in set(uniprot_mnemonics)}

    log.info("UniProt→RefSeq: mapping %d accessions via UniProt search API...", len(to_map))

    _BATCH = 100
    for i in range(0, len(to_map), _BATCH):
        batch = to_map[i:i + _BATCH]
        query = " OR ".join(f"id:{m}" for m in batch)
        params = urlencode({
            "query": f"({query})",
            "fields": "accession,id,xref_refseq",
            "format": "json",
            "size": str(len(batch)),
        })
        url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
        try:
            req = Request(url, headers={"Accept": "application/json"})
            with urlopen(req, timeout=120) as resp:
                data = json.loads(resp.read().decode())

            for entry in data.get("results", []):
                entry_name = entry.get("uniProtkbId", "")
                xrefs = entry.get("uniProtKBCrossReferences", [])
                for xref in xrefs:
                    if xref.get("database") == "RefSeq":
                        refseq_id = xref.get("id", "")
                        if refseq_id.startswith("NP_"):
                            mapping[entry_name] = refseq_id
                            break

            if (i + _BATCH) < len(to_map):
                time.sleep(0.5)

        except Exception as exc:
            log.warning("UniProt search batch %d failed: %s", i, exc)

    try:
        map_cache.write_text(json.dumps(mapping))
    except Exception:
        pass

    mapped = {k: v for k, v in mapping.items() if k in set(uniprot_mnemonics)}
    log.info("UniProt→RefSeq: %d/%d mapped successfully.", len(mapped), len(uniprot_mnemonics))
    return mapped


def _prefetch_all_cds(aligned_orthogroups: dict[str, dict[str, str]]) -> None:
    """Pre-fetch CDS for all unique protein accessions before spawning workers.

    Strategy:
    - Phase 0: Map UniProt mnemonics (human) → NCBI RefSeq protein IDs
    - Phase 1: efetch fasta_cds_na in parallel threads (8 workers, batches of 200)
      fasta_cds_na returns only the CDS FASTA — much smaller/faster than GenBank
    - Results cached to disk so workers never hit NCBI again
    """
    import threading

    _setup_entrez()
    api_key = get_ncbi_api_key()
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    cache_dir = Path(get_local_storage_root()) / "cds"
    cache_dir.mkdir(parents=True, exist_ok=True)

    _UNIPROT_PAT = re.compile(
        r"^[A-Z0-9]{1,11}_[A-Z]{3,5}$"
        r"|^[OPQ][0-9][A-Z0-9]{3}[0-9]$"
        r"|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
    )

    def _is_uniprot(acc: str) -> bool:
        if "_" in acc:
            after = acc.split("_", 1)[1].lstrip("0123456789.")
            if after:
                return bool(_UNIPROT_PAT.match(acc))
            return False
        return bool(_UNIPROT_PAT.match(acc))

    all_accessions: set[str] = set()
    for seqs in aligned_orthogroups.values():
        for label in seqs:
            acc = label.split("|")[-1] if "|" in label else label
            all_accessions.add(acc)

    uniprot_accs = sorted(a for a in all_accessions if _is_uniprot(a))
    uniprot_needing_map = []
    for a in uniprot_accs:
        cf = cache_dir / f"{a}.fna"
        if not cf.exists() or cf.stat().st_size == 0:
            uniprot_needing_map.append(a)
    uniprot_to_refseq: dict[str, str] = {}
    if uniprot_needing_map:
        uniprot_to_refseq = _map_uniprot_to_refseq(uniprot_needing_map, cache_dir)

    to_fetch = []
    for acc in sorted(all_accessions):
        cf = cache_dir / f"{acc}.fna"
        if _is_uniprot(acc):
            refseq = uniprot_to_refseq.get(acc)
            if refseq:
                to_fetch.append((acc, refseq))
            elif not cf.exists():
                cf.write_text("")
            continue
        if cf.exists() and cf.stat().st_size == 0:
            cf.unlink()
        if not cf.exists():
            to_fetch.append((acc, acc))

    if not to_fetch:
        log.info("CDS pre-fetch: all %d accessions already cached.", len(all_accessions))
        return

    ncbi_direct = sum(1 for _, r in to_fetch if not _is_uniprot(r))
    uniprot_mapped = sum(1 for o, r in to_fetch if o != r)
    log.info("CDS pre-fetch: %d to fetch (%d NCBI direct, %d UniProt→RefSeq mapped).",
             len(to_fetch), ncbi_direct, uniprot_mapped)

    # Build mapping: fetch_acc → [original_accs] so we can cache under original names
    fetch_to_originals: dict[str, list[str]] = {}
    for orig, fetch in to_fetch:
        fetch_to_originals.setdefault(fetch, []).append(orig)

    fetch_ids = list(fetch_to_originals.keys())

    _EFETCH_BATCH = 200
    _EFETCH_WORKERS = 8
    done_fetch = 0
    total_written = 0
    fetch_lock = threading.Lock()

    _pid_pat = re.compile(r"\[protein_id=([^\]]+)\]")

    def _flush_cds(fetch_acc: str, lines: list[str]) -> int:
        seq = "".join(lines)
        count = 0
        for orig_acc in fetch_to_originals.get(fetch_acc, [fetch_acc]):
            cache_path = cache_dir / f"{orig_acc}.fna"
            existing = cache_path.read_text() if cache_path.exists() else ""
            if len(seq) > len(existing):
                cache_path.write_text(seq if seq else "")
            elif not cache_path.exists():
                cache_path.write_text(seq if seq else "")
            if seq:
                count += 1
        return count

    def _fetch_cds_batch_direct(batch: list[str]) -> int:
        """Fetch fasta_cds_na for a batch of protein accessions with retry on 429."""
        _MAX_RETRIES = 5
        for attempt in range(_MAX_RETRIES):
            params = {
                "db": "protein",
                "id": ",".join(batch),
                "rettype": "fasta_cds_na",
                "retmode": "text",
            }
            if api_key:
                params["api_key"] = api_key
            try:
                req = Request(
                    base_url + "efetch.fcgi",
                    data=urlencode(params).encode(),
                    method="POST",
                )
                with urlopen(req, timeout=120) as resp:
                    text = resp.read().decode("utf-8", errors="replace")
                written = 0
                current_fetch_acc = None
                current_lines: list[str] = []

                for line in text.splitlines():
                    if line.startswith(">"):
                        if current_fetch_acc and current_lines:
                            written += _flush_cds(current_fetch_acc, current_lines)
                        current_lines = []
                        current_fetch_acc = None
                        m = _pid_pat.search(line)
                        if m:
                            pid = m.group(1)
                            if pid in batch:
                                current_fetch_acc = pid
                            else:
                                pid_base = pid.rsplit(".", 1)[0]
                                for b in batch:
                                    if b.rsplit(".", 1)[0] == pid_base:
                                        current_fetch_acc = b
                                        break
                    else:
                        if current_fetch_acc:
                            current_lines.append(line.strip())

                if current_fetch_acc and current_lines:
                    written += _flush_cds(current_fetch_acc, current_lines)

                for fetch_acc in batch:
                    for orig_acc in fetch_to_originals.get(fetch_acc, [fetch_acc]):
                        if not (cache_dir / f"{orig_acc}.fna").exists():
                            (cache_dir / f"{orig_acc}.fna").write_text("")
                return written
            except Exception as exc:
                is_retryable = any(s in str(exc) for s in ("429", "Network is unreachable",
                    "Connection reset", "timed out", "503", "500", "urlopen error"))
                if is_retryable and attempt < _MAX_RETRIES - 1:
                    wait = (2 ** attempt) * 5 + _rng.uniform(0, 5)
                    log.info("NCBI transient error, retry %d/%d in %.0fs: %s",
                             attempt + 1, _MAX_RETRIES, wait, exc)
                    time.sleep(wait)
                    continue
                log.warning("efetch batch failed (attempt %d/%d): %s",
                            attempt + 1, _MAX_RETRIES, exc)
                for fetch_acc in batch:
                    for orig_acc in fetch_to_originals.get(fetch_acc, [fetch_acc]):
                        if not (cache_dir / f"{orig_acc}.fna").exists():
                            (cache_dir / f"{orig_acc}.fna").write_text("")
                return 0
        log.warning("efetch batch exhausted all %d retries", _MAX_RETRIES)
        return 0

    log.info("  CDS pre-fetch: efetch fasta_cds_na from protein db, "
             "%d fetch-accessions, %d workers, batches of %d...",
             len(fetch_ids), _EFETCH_WORKERS, _EFETCH_BATCH)

    efetch_batches = [fetch_ids[i:i + _EFETCH_BATCH] for i in range(0, len(fetch_ids), _EFETCH_BATCH)]
    with ThreadPoolExecutor(max_workers=_EFETCH_WORKERS) as pool:
        futures2 = {pool.submit(_fetch_cds_batch_direct, b): b for b in efetch_batches}
        for fut in as_completed(futures2):
            written = fut.result()
            with fetch_lock:
                total_written += written
                done_fetch += len(futures2[fut])
                if done_fetch % 2000 == 0 or done_fetch >= len(fetch_ids):
                    log.info("  efetch: %d / %d done, %d CDS written.",
                             done_fetch, len(fetch_ids), total_written)
            time.sleep(0.02)

    log.info("CDS pre-fetch complete: %d written, %d failed/no-CDS.",
             total_written, len(fetch_ids) - total_written)


def fetch_cds_for_protein(protein_accession: str) -> Optional[str]:
    """Fetch the CDS nucleotide sequence for a protein accession from NCBI.

    Checks in-memory cache first (populated by load_cds_cache_pkl),
    then falls back to disk cache, then NCBI as last resort.
    """
    acc = protein_accession.split("|")[-1] if "|" in protein_accession else protein_accession

    if acc in _IN_MEMORY_CDS_CACHE:
        seq = _IN_MEMORY_CDS_CACHE[acc]
        return seq if seq else None

    _setup_entrez()

    cache_dir = Path(get_local_storage_root()) / "cds"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{acc}.fna"
    if cache_file.exists():
        cds = cache_file.read_text().strip()
        return cds if cds else None

    if re.match(r"^[A-Z0-9]{1,11}_[A-Z]{3,5}$", acc):
        return None
    if re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9])", acc):
        return None

    try:
        # Link protein → nucleotide
        handle = Entrez.elink(dbfrom="protein", db="nuccore", id=acc, linkname="protein_nuccore_mrna")
        link_record = Entrez.read(handle)
        handle.close()
        time.sleep(0.12)

        nuc_ids = []
        if link_record and link_record[0].get("LinkSetDb"):
            for link in link_record[0]["LinkSetDb"]:
                if link["LinkName"] in ("protein_nuccore_mrna", "protein_nuccore"):
                    nuc_ids = [l["Id"] for l in link["Link"]]
                    break

        if not nuc_ids:
            cache_file.write_text("")   # cache negative result
            return None

        # Fetch the first mRNA record
        handle = Entrez.efetch(db="nuccore", id=nuc_ids[0], rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        time.sleep(0.12)

        # Extract CDS feature
        for feature in record.features:
            if feature.type == "CDS":
                cds_seq = str(feature.extract(record.seq))
                # Validate: must be divisible by 3 and start with ATG
                if len(cds_seq) % 3 == 0 and cds_seq.upper().startswith("ATG"):
                    cache_file.write_text(cds_seq)   # persist for future runs
                    return cds_seq

        cache_file.write_text("")   # CDS not found — cache to skip next time
        return None

    except Exception as exc:
        log.debug("CDS fetch failed for %s: %s", acc, exc)
        return None


# ---------------------------------------------------------------------------
# Codon alignment (PAL2NAL-style)
# ---------------------------------------------------------------------------

def protein_to_codon_alignment(
    protein_alignment: dict[str, str],
    cds_seqs: dict[str, str],
    min_coverage: float = 0.60,
    min_seqs: int = 3,
) -> Optional[dict[str, str]]:
    """Build a codon alignment from a protein alignment and CDS sequences.

    Resilient: skips individual sequences whose CDS length doesn't match
    (different isoform, trimmed protein, etc.) and proceeds with the rest
    as long as coverage and count thresholds are met.

    Args:
        protein_alignment: {label: aligned_protein_seq} — contains gap characters
        cds_seqs: {label: cds_nucleotide_seq} — no gaps, raw CDS
        min_coverage: fraction of sequences that must succeed (default 0.60)
        min_seqs: minimum number of valid codon-aligned sequences (default 3)

    Returns:
        {label: codon_aligned_nucleotide_seq} or None if thresholds not met.
    """
    codon_aln: dict[str, str] = {}
    skipped: list[str] = []

    for label, prot_aln in protein_alignment.items():
        cds = cds_seqs.get(label)
        if not cds:
            skipped.append(label)
            continue

        non_gap_count = sum(1 for aa in prot_aln if aa != "-")
        expected_cds_len = non_gap_count * 3

        if len(cds) < expected_cds_len:
            log.info("Codon skip %s: CDS len %d < expected %d (protein has %d non-gap AAs)",
                     label, len(cds), expected_cds_len, non_gap_count)
            skipped.append(label)
            continue

        cds_pos = 0
        codon_seq = []
        for aa in prot_aln:
            if aa == "-":
                codon_seq.append("---")
            else:
                codon_seq.append(cds[cds_pos:cds_pos + 3])
                cds_pos += 3

        codon_aln[label] = "".join(codon_seq)

    if skipped:
        log.info("Codon alignment: %d/%d succeeded, %d skipped",
                 len(codon_aln), len(protein_alignment), len(skipped))

    _STOP_CODONS = {"TAA", "TAG", "TGA"}
    clean_aln: dict[str, str] = {}
    for label, seq in codon_aln.items():
        has_internal_stop = False
        last_codon_start = max(0, len(seq) - 3)
        for i in range(0, last_codon_start, 3):
            codon = seq[i:i + 3].upper().replace("-", "")
            if codon in _STOP_CODONS:
                has_internal_stop = True
                break
        if has_internal_stop:
            log.info("Codon alignment: dropping %s (internal stop codon)", label)
        else:
            clean_aln[label] = seq
    if len(clean_aln) < len(codon_aln):
        log.info("Dropped %d seqs with stop codons, %d remain",
                 len(codon_aln) - len(clean_aln), len(clean_aln))
    codon_aln = clean_aln

    coverage = len(codon_aln) / len(protein_alignment) if protein_alignment else 0
    if len(codon_aln) < min_seqs or coverage < min_coverage:
        log.info("Codon alignment below threshold: %d seqs (need %d), %.0f%% (need %.0f%%)",
                 len(codon_aln), min_seqs, coverage * 100, min_coverage * 100)
        return None

    lengths = {len(v) for v in codon_aln.values()}
    if len(lengths) != 1:
        log.info("Codon alignment has unequal lengths: %s — trimming to common length", lengths)
        min_len = min(lengths)
        codon_aln = {k: v[:min_len] for k, v in codon_aln.items()}

    # --- Sequence QC ---

    # Remove sequences with >50% gap characters
    gap_filtered = {}
    for label, seq in codon_aln.items():
        gap_frac = seq.count("-") / len(seq) if seq else 1.0
        if gap_frac <= 0.50:
            gap_filtered[label] = seq
        else:
            log.info("Codon QC: dropping %s (%.0f%% gaps)", label, gap_frac * 100)
    codon_aln = gap_filtered

    # Remove sequences with >5% ambiguous bases
    _VALID_BASES = set("ATGCatgc-")
    clean2 = {}
    for label, seq in codon_aln.items():
        ambig_count = sum(1 for c in seq if c not in _VALID_BASES)
        ambig_frac = ambig_count / len(seq) if seq else 1.0
        if ambig_frac <= 0.05:
            clean2[label] = seq
        else:
            log.info("Codon QC: dropping %s (%.1f%% ambiguous bases)", label, ambig_frac * 100)
    codon_aln = clean2

    # Remove duplicate sequences (identical codons = same evolutionary signal)
    seen_seqs: dict[str, str] = {}
    deduped: dict[str, str] = {}
    for label, seq in codon_aln.items():
        seq_upper = seq.upper()
        if seq_upper not in seen_seqs:
            seen_seqs[seq_upper] = label
            deduped[label] = seq
        else:
            log.info("Codon QC: dropping %s (identical to %s)", label, seen_seqs[seq_upper])
    codon_aln = deduped

    if len(codon_aln) < min_seqs:
        log.info("Codon QC: only %d sequences remain after QC (need %d)", len(codon_aln), min_seqs)
        return None

    return codon_aln


def infer_gene_tree(codon_aln: dict[str, str], og_id: str) -> Optional[str]:
    """Infer a per-gene ML tree using IQ-TREE2 on the codon alignment.

    Uses GTR+G (general time-reversible + gamma rate heterogeneity) with
    single-threaded execution for efficiency when parallelised across genes.

    Returns the Newick tree string, or None on failure.
    """
    import shutil as _shutil

    root = Path(get_local_storage_root())
    tree_dir = root / "gene_trees" / og_id
    tree_dir.mkdir(parents=True, exist_ok=True)

    aln_path = tree_dir / "codon_aln.fna"
    with open(aln_path, "w") as f:
        for label, seq in codon_aln.items():
            species = label.split("|")[0]
            f.write(f">{species}\n{seq}\n")

    prefix = str(tree_dir / "gene")
    treefile = Path(f"{prefix}.treefile")

    if treefile.exists() and treefile.stat().st_size > 0:
        return treefile.read_text().strip()

    iqtree_bin = "iqtree2" if _shutil.which("iqtree2") else "iqtree"
    if not _shutil.which(iqtree_bin):
        log.warning("IQ-TREE2 not found — cannot infer gene tree for %s", og_id)
        return None

    cmd = [
        iqtree_bin,
        "-s", str(aln_path),
        "-m", "GTR+G",
        "-T", "1",
        "--prefix", prefix,
        "--redo",
        "-quiet",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0 and treefile.exists():
            return treefile.read_text().strip()
        log.warning("IQ-TREE failed for %s (rc=%d): %s",
                    og_id, result.returncode, result.stderr[:200])
    except subprocess.TimeoutExpired:
        log.warning("IQ-TREE timed out for %s (>120s)", og_id)
    except Exception as exc:
        log.warning("IQ-TREE error for %s: %s", og_id, exc)

    return None


# ---------------------------------------------------------------------------
# HyPhy MEME runner
# ---------------------------------------------------------------------------

def run_meme(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
) -> Optional[dict]:
    """Run HyPhy MEME on a codon alignment and return the JSON results.

    MEME (Mixed Effects Model of Evolution) detects sites under episodic
    positive selection — the correct model for convergent adaptive evolution.

    Returns the parsed MEME JSON output, or None on failure.
    """
    tools = get_tool_config()
    hyphy_bin = tools.get("hyphy_bin", "hyphy")

    root = Path(get_local_storage_root())
    meme_dir = root / "meme" / og_id
    meme_dir.mkdir(parents=True, exist_ok=True)

    aln_path = meme_dir / "codon_aln.fna"
    tree_path = meme_dir / "species.treefile"
    out_path = meme_dir / "meme.json"

    species_map = {}
    with open(aln_path, "w") as f:
        for label, seq in codon_aln.items():
            species = label.split("|")[0]
            if species in species_map:
                log.warning("Duplicate species %s in OG %s, skipping %s", species, og_id, label)
                continue
            species_map[species] = label
            f.write(f">{species}\n{seq}\n")

    tree_path.write_text(tree_newick)

    cpus = int(os.environ.get("HYPHY_CPUS", 1))
    hyphy_env = {**os.environ, "CPU": str(cpus)}
    meme_timeout = int(os.environ.get("MEME_TIMEOUT", 180))

    try:
        cmd = [
            hyphy_bin, "meme",
            "--alignment", str(aln_path),
            "--tree", str(tree_path),
            "--output", str(out_path),
            "--branches", "All",
        ]
        log.info("Running MEME for %s (CPU=%d, timeout=%ds): %s", og_id, cpus, meme_timeout, " ".join(cmd[:3]))
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=meme_timeout,
            env=hyphy_env,
        )
        if result.returncode != 0:
            log.warning("MEME failed for %s (rc=%d): %s",
                        og_id, result.returncode, result.stderr[:500])
            return None

        if not out_path.exists():
            log.warning("MEME output missing for %s", og_id)
            return None

        log.info("MEME succeeded for %s", og_id)
        return json.loads(out_path.read_text())

    except subprocess.TimeoutExpired:
        log.warning("MEME timed out for %s (>%ds)", og_id, meme_timeout)
        return None
    except Exception as exc:
        log.warning("MEME error for %s: %s", og_id, exc)
        return None


# ---------------------------------------------------------------------------
# MEME result parsing
# ---------------------------------------------------------------------------

def _extract_meme_resistant_branches(meme_json: dict) -> list[str]:
    """Return species IDs of test (resistant/target) branches with elevated dN/dS in MEME.

    Reads the `branch attributes` section of the MEME JSON — HyPhy always writes
    per-branch MG94 omega estimates there.  We intersect the elevated-omega branches
    with the registry-derived test-species set so the filter is phenotype-agnostic.

    This is the raw MEME output preserved for reference; the caller decides how to
    use this list (e.g. as a phenotype-specificity weight in scoring).
    """
    test_set, _ref, _out = _get_relax_branch_sets()
    if not test_set:
        return []

    elevated: list[str] = []
    branch_attrs = meme_json.get("branch attributes", {})
    # Branch attributes may be nested under a model key "0"
    if isinstance(branch_attrs, dict) and "0" in branch_attrs:
        branch_attrs = branch_attrs["0"]

    for branch_name, attrs in branch_attrs.items():
        if not isinstance(attrs, dict):
            continue
        # Try various key names HyPhy uses for per-branch omega across versions
        omega = (
            attrs.get("MG94xREV with separate rates for branch sets")
            or attrs.get("omega")
            or attrs.get("dN/dS")
            or attrs.get("Omega")
        )
        if omega is None:
            continue
        try:
            if float(omega) > 1.0:
                # Strip any HyPhy branch annotations like "{Test}" or "{Reference}"
                clean = branch_name.split("{")[0].strip()
                if clean in test_set:
                    elevated.append(clean)
        except (ValueError, TypeError):
            continue

    return elevated


def parse_meme_results(meme_json: dict, og_id: str) -> dict:
    """Extract selection statistics from HyPhy MEME output.

    MEME outputs per-site results with:
      - alpha: synonymous rate
      - beta-: purifying selection intensity
      - beta+: positive selection intensity
      - p-value: episodic diversification p-value per site

    Also extracts which test (resistant/target) species branches show elevated
    omega (dN/dS > 1) from the `branch attributes` section, stored in
    `branches_under_selection`.  Raw MEME output is preserved in full; the
    scoring layer uses branch specificity as a phenotype-relevance multiplier.

    Returns a dict compatible with load_selection_scores().
    """
    try:
        mle_content = meme_json.get("MLE", {}).get("content", {}).get("0", [])
        if not mle_content:
            return _meme_null_result(og_id)

        n_sites = len(mle_content)
        selected_sites = []
        beta_plus_values = []

        for site_data in mle_content:
            # MEME MLE columns: [alpha, beta-, omega-, beta+, omega+, p-value, ...]
            if len(site_data) < 6:
                continue
            pvalue = site_data[5]
            beta_plus = site_data[3]
            if pvalue is not None and pvalue < MEME_SIGNIFICANCE_THRESHOLD:
                selected_sites.append(pvalue)
                if beta_plus is not None:
                    beta_plus_values.append(beta_plus)

        if n_sites == 0:
            return _meme_null_result(og_id)

        fraction_selected = len(selected_sites) / n_sites
        mean_beta_plus = sum(beta_plus_values) / len(beta_plus_values) if beta_plus_values else 0.0

        pseudo_dnds = fraction_selected * min(mean_beta_plus / 5.0, 2.0)
        if selected_sites:
            import math
            log_sum = sum(math.log(p) for p in selected_sites if p > 0)
            pseudo_pvalue = max(math.exp(log_sum / len(selected_sites)), 1e-10)
        else:
            pseudo_pvalue = 1.0

        # Identify which resistant/target branches carry the elevated signal.
        # Kept separate from per-site filtering so raw MEME output is untouched.
        resistant_branches = _extract_meme_resistant_branches(meme_json)

        log.debug("  MEME %s: %d/%d sites selected, mean_beta+=%.2f, pseudo_dnds=%.3f, "
                  "resistant_branches=%s",
                  og_id, len(selected_sites), n_sites, mean_beta_plus, pseudo_dnds,
                  resistant_branches)

        return {
            "dnds_ratio": round(pseudo_dnds, 4),
            "dnds_pvalue": round(pseudo_pvalue, 6),
            "selection_model": "MEME_episodic",
            "branches_under_selection": resistant_branches,
            "fraction_sites_selected": round(fraction_selected, 4),
            "n_sites_selected": len(selected_sites),
            "mean_beta_plus": round(mean_beta_plus, 4),
        }

    except Exception as exc:
        log.debug("MEME parse error for %s: %s", og_id, exc)
        return _meme_null_result(og_id)


def _meme_null_result(og_id: str) -> dict:
    return {
        "dnds_ratio": 0.0,
        "dnds_pvalue": 1.0,
        "selection_model": "MEME_no_signal",
        "branches_under_selection": [],
        "fraction_sites_selected": 0.0,
        "n_sites_selected": 0,
        "mean_beta_plus": 0.0,
    }


# ---------------------------------------------------------------------------
# FUBAR — Fast Unconstrained Bayesian AppRoximation (pervasive selection pre-screen)
# ---------------------------------------------------------------------------

def run_fubar(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
) -> Optional[dict]:
    """Run HyPhy FUBAR on a codon alignment.

    FUBAR is ~10-20x faster than MEME and detects *pervasive* positive selection.
    Used as a pre-screen: only OGs where FUBAR finds signal are sent to MEME.
    """
    tools = get_tool_config()
    hyphy_bin = tools.get("hyphy_bin", "hyphy")

    root = Path(get_local_storage_root())
    fubar_dir = root / "fubar" / og_id
    fubar_dir.mkdir(parents=True, exist_ok=True)

    aln_path = fubar_dir / "codon_aln.fna"
    tree_path = fubar_dir / "species.treefile"
    out_path = fubar_dir / "fubar.json"

    if out_path.exists() and out_path.stat().st_size > 0:
        try:
            return json.loads(out_path.read_text())
        except Exception:
            pass

    with open(aln_path, "w") as f:
        seen_species = set()
        for label, seq in codon_aln.items():
            species = label.split("|")[0]
            if species in seen_species:
                continue
            seen_species.add(species)
            f.write(f">{species}\n{seq}\n")
    tree_path.write_text(tree_newick)

    cpus = int(os.environ.get("HYPHY_CPUS", 1))
    hyphy_env = {**os.environ, "CPU": str(cpus)}

    try:
        cmd = [
            hyphy_bin, "fubar",
            "--alignment", str(aln_path),
            "--tree", str(tree_path),
            "--output", str(out_path),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120, env=hyphy_env)
        if result.returncode != 0:
            log.debug("FUBAR failed for %s (rc=%d)", og_id, result.returncode)
            return None
        if not out_path.exists():
            return None
        return json.loads(out_path.read_text())
    except subprocess.TimeoutExpired:
        log.debug("FUBAR timed out for %s", og_id)
        return None
    except Exception as exc:
        log.debug("FUBAR error for %s: %s", og_id, exc)
        return None


def parse_fubar_results(fubar_json: dict) -> dict:
    """Count sites under pervasive positive selection from HyPhy FUBAR output.

    Returns {"fubar_pos_sites": int, "fubar_neg_sites": int}.
    Positive selection: posterior probability > 0.9 that beta > alpha.
    """
    try:
        mle = fubar_json.get("MLE", {})
        mle_content = mle.get("content", {}).get("0", [])
        if not mle_content:
            return {"fubar_pos_sites": 0, "fubar_neg_sites": 0}

        pos_sites = 0
        neg_sites = 0
        for site_data in mle_content:
            if len(site_data) < 5:
                continue
            # FUBAR columns: [alpha, beta, beta-alpha, Prob(alpha<beta), Prob(alpha>beta)]
            prob_pos = site_data[3] if len(site_data) > 3 else 0
            prob_neg = site_data[4] if len(site_data) > 4 else 0
            if prob_pos is not None and prob_pos > 0.9:
                pos_sites += 1
            if prob_neg is not None and prob_neg > 0.9:
                neg_sites += 1

        return {"fubar_pos_sites": pos_sites, "fubar_neg_sites": neg_sites}
    except Exception as exc:
        log.debug("FUBAR parse error: %s", exc)
        return {"fubar_pos_sites": 0, "fubar_neg_sites": 0}


def _codon_aln_passes_qc(codon_aln: dict[str, str], og_id: str) -> bool:
    """Pre-flight QC: reject alignments that will waste time in MEME."""
    if len(codon_aln) < 3:
        log.debug("Skipping %s: only %d sequences", og_id, len(codon_aln))
        return False

    aln_len = len(next(iter(codon_aln.values()), ""))
    n_codons = aln_len // 3
    if n_codons < 100:
        log.debug("Skipping %s: only %d codons (need ≥100)", og_id, n_codons)
        return False

    total_chars = sum(len(s) for s in codon_aln.values())
    total_gaps = sum(s.count("-") for s in codon_aln.values())
    gap_frac = total_gaps / total_chars if total_chars else 1.0
    if gap_frac > 0.40:
        log.debug("Skipping %s: %.0f%% gaps in alignment", og_id, gap_frac * 100)
        return False

    return True


# ---------------------------------------------------------------------------
# FEL — Fixed Effects Likelihood (pervasive selection)
# ---------------------------------------------------------------------------

def run_fel(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
) -> Optional[dict]:
    """Run HyPhy FEL on a codon alignment.

    FEL (Fixed Effects Likelihood) detects sites under *pervasive* (constant)
    positive or purifying selection — complementary to MEME's episodic detection.

    Returns the parsed FEL JSON output, or None on failure.
    """
    tools = get_tool_config()
    hyphy_bin = tools.get("hyphy_bin", "hyphy")

    root = Path(get_local_storage_root())
    fel_dir = root / "fel" / og_id
    fel_dir.mkdir(parents=True, exist_ok=True)

    aln_path = fel_dir / "codon_aln.fna"
    tree_path = fel_dir / "species.treefile"
    out_path = fel_dir / "fel.json"

    with open(aln_path, "w") as f:
        seen_species = set()
        for label, seq in codon_aln.items():
            species = label.split("|")[0]
            if species in seen_species:
                continue
            seen_species.add(species)
            f.write(f">{species}\n{seq}\n")
    tree_path.write_text(tree_newick)

    cpus = int(os.environ.get("HYPHY_CPUS", 1))
    hyphy_env = {**os.environ, "CPU": str(cpus)}

    try:
        cmd = [
            hyphy_bin, "fel",
            "--alignment", str(aln_path),
            "--tree", str(tree_path),
            "--output", str(out_path),
            "--full-model", "No",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=hyphy_env)
        if result.returncode != 0:
            log.warning("FEL failed for %s (rc=%d): %s", og_id, result.returncode, result.stderr[:300])
            return None
        if not out_path.exists():
            log.warning("FEL output missing for %s", og_id)
            return None
        return json.loads(out_path.read_text())
    except subprocess.TimeoutExpired:
        log.debug("FEL timed out for %s", og_id)
        return None
    except Exception as exc:
        log.debug("FEL error for %s: %s", og_id, exc)
        return None


def parse_fel_results(fel_json: dict) -> dict:
    """Count sites under pervasive positive selection from HyPhy FEL output.

    Dynamically detects column positions from FEL JSON headers to handle
    format differences across HyPhy versions (2.5.33+ added an extra column).

    Returns {"fel_sites": int} — count of sites with dN > dS at p < 0.05.
    """
    try:
        mle = fel_json.get("MLE", {})
        mle_content = mle.get("content", {}).get("0", [])
        if not mle_content:
            return {"fel_sites": 0}

        headers = mle.get("headers", [])
        header_names = [h[0].lower() if isinstance(h, list) else str(h).lower() for h in headers]

        alpha_idx = 0
        beta_idx = 1
        pvalue_idx = 4

        for i, name in enumerate(header_names):
            if "alpha" in name:
                alpha_idx = i
            elif "beta" in name and "alpha" not in name:
                beta_idx = i
            elif "p-value" in name or name == "p-value" or "p_value" in name:
                pvalue_idx = i

        pos_sites = 0
        for site_data in mle_content:
            if len(site_data) <= max(alpha_idx, beta_idx, pvalue_idx):
                continue
            alpha = site_data[alpha_idx]
            beta = site_data[beta_idx]
            pvalue = site_data[pvalue_idx]
            if pvalue is not None and pvalue < 0.05 and beta is not None and alpha is not None:
                if beta > alpha:
                    pos_sites += 1

        return {"fel_sites": pos_sites}

    except Exception as exc:
        log.debug("FEL parse error: %s", exc)
        return {"fel_sites": 0}


# ---------------------------------------------------------------------------
# BUSTED — Branch-Site Unrestricted Test for Episodic Diversification
# ---------------------------------------------------------------------------

def run_busted(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
) -> Optional[dict]:
    """Run HyPhy BUSTED on a codon alignment.

    BUSTED tests whether the gene *as a whole* has experienced episodic positive
    selection anywhere in the tree, providing a single gene-level p-value.
    Complementary to MEME (which is site-level).

    Returns the parsed BUSTED JSON output, or None on failure.
    """
    tools = get_tool_config()
    hyphy_bin = tools.get("hyphy_bin", "hyphy")

    root = Path(get_local_storage_root())
    busted_dir = root / "busted" / og_id
    busted_dir.mkdir(parents=True, exist_ok=True)

    aln_path = busted_dir / "codon_aln.fna"
    tree_path = busted_dir / "species.treefile"
    out_path = busted_dir / "busted.json"

    with open(aln_path, "w") as f:
        seen_species = set()
        for label, seq in codon_aln.items():
            species = label.split("|")[0]
            if species in seen_species:
                continue
            seen_species.add(species)
            f.write(f">{species}\n{seq}\n")
    tree_path.write_text(tree_newick)

    cpus = int(os.environ.get("HYPHY_CPUS", 1))
    hyphy_env = {**os.environ, "CPU": str(cpus)}

    try:
        cmd = [
            hyphy_bin, "busted",
            "--alignment", str(aln_path),
            "--tree", str(tree_path),
            "--output", str(out_path),
            # --srv Yes omitted: SRV is 2-4x slower with marginal sensitivity gain at pipeline scale
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=hyphy_env)
        if result.returncode != 0:
            log.warning("BUSTED failed for %s (rc=%d): %s", og_id, result.returncode, result.stderr[:300])
            return None
        if not out_path.exists():
            return None
        return json.loads(out_path.read_text())
    except subprocess.TimeoutExpired:
        log.debug("BUSTED timed out for %s", og_id)
        return None
    except Exception as exc:
        log.debug("BUSTED error for %s: %s", og_id, exc)
        return None


def parse_busted_results(busted_json: dict) -> dict:
    """Extract gene-level p-value from HyPhy BUSTED output.

    Returns {"busted_pvalue": float} — p-value for gene-wide episodic selection.
    """
    try:
        # BUSTED stores test p-value at json["test results"]["p-value"]
        pvalue = (
            busted_json
            .get("test results", {})
            .get("p-value")
        )
        if pvalue is None:
            return {"busted_pvalue": 1.0}
        return {"busted_pvalue": round(float(pvalue), 6)}
    except Exception as exc:
        log.debug("BUSTED parse error: %s", exc)
        return {"busted_pvalue": 1.0}


# ---------------------------------------------------------------------------
# BUSTED-PH — BUSTED restricted to foreground (phenotype-bearing) branches
# ---------------------------------------------------------------------------

def run_busted_ph(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
) -> Optional[dict]:
    """Run HyPhy BUSTED restricted to test (phenotype-bearing) foreground branches.

    This is the phenotype-specific equivalent of MEME — ~50-200x faster because
    it computes a single gene-level LRT instead of per-site tests.  It answers:
    "Is there evidence of positive selection specifically in the species that carry
    the phenotype of interest?"

    Implemented as standard HyPhy BUSTED with `--branches Test`, where test-species
    branches are annotated `{Test}` in the Newick tree.  The LRT compares:
      H0: the foreground rate distribution = background
      H1: foreground is allowed a positive ω component (ω > 1)

    Returns the parsed BUSTED JSON, or None on failure.
    """
    tools = get_tool_config()
    hyphy_bin = tools.get("hyphy_bin", "hyphy")

    root = Path(get_local_storage_root())
    bph_dir = root / "busted_ph" / og_id
    bph_dir.mkdir(parents=True, exist_ok=True)

    aln_path = bph_dir / "codon_aln.fna"
    tree_path = bph_dir / "species.treefile"
    out_path = bph_dir / "busted_ph.json"

    if out_path.exists() and out_path.stat().st_size > 0:
        try:
            return json.loads(out_path.read_text())
        except Exception:
            pass

    with open(aln_path, "w") as f:
        seen_species = set()
        for label, seq in codon_aln.items():
            species = label.split("|")[0]
            if species in seen_species:
                continue
            seen_species.add(species)
            f.write(f">{species}\n{seq}\n")

    # Label test-phenotype branches with {Test} so BUSTED can restrict to them.
    # Reference and outgroup branches are unlabeled (treated as background).
    test_set, _ref, _out = _get_relax_branch_sets()
    all_aln_species = list({label.split("|")[0] for label in codon_aln})
    test_branches = [sp for sp in all_aln_species if sp in test_set]

    if not test_branches:
        log.info("BUSTED-PH %s: no test branches in alignment — skipping", og_id)
        return None

    labeled_tree = tree_newick
    for branch in sorted(test_branches, key=len, reverse=True):
        labeled_tree = labeled_tree.replace(f"{branch}:", f"{branch}{{Test}}:")
    tree_path.write_text(labeled_tree)

    cpus = int(os.environ.get("HYPHY_CPUS", 1))
    hyphy_env = {**os.environ, "CPU": str(cpus)}

    try:
        cmd = [
            hyphy_bin, "busted",
            "--alignment", str(aln_path),
            "--tree", str(tree_path),
            "--output", str(out_path),
            "--branches", "Test",
            # --srv Yes omitted: SRV is 2-4x slower with marginal sensitivity gain at pipeline scale
        ]
        log.info("Running BUSTED-PH for %s (CPU=%d)", og_id, cpus)
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=hyphy_env)
        if result.returncode != 0:
            log.warning("BUSTED-PH failed for %s (rc=%d): %s",
                        og_id, result.returncode, result.stderr[:300])
            return None
        if not out_path.exists():
            log.warning("BUSTED-PH output missing for %s", og_id)
            return None
        return json.loads(out_path.read_text())
    except subprocess.TimeoutExpired:
        log.warning("BUSTED-PH timed out for %s (>600s)", og_id)
        return None
    except Exception as exc:
        log.warning("BUSTED-PH error for %s: %s", og_id, exc)
        return None


def parse_busted_ph_results(busted_ph_json: dict) -> dict:
    """Extract phenotype-specific gene-level p-value from BUSTED-PH output.

    Returns {"busted_ph_pvalue": float} — p-value for positive selection
    specifically in the foreground (phenotype-bearing) branches.
    """
    try:
        pvalue = (
            busted_ph_json
            .get("test results", {})
            .get("p-value")
        )
        if pvalue is None:
            return {"busted_ph_pvalue": 1.0}
        return {"busted_ph_pvalue": round(float(pvalue), 6)}
    except Exception as exc:
        log.debug("BUSTED-PH parse error: %s", exc)
        return {"busted_ph_pvalue": 1.0}


# ---------------------------------------------------------------------------
# FEL + BUSTED pipeline entry point
# ---------------------------------------------------------------------------

def _fel_busted_worker(args: tuple) -> tuple[str, dict]:
    """Top-level worker: run FEL + BUSTED for one orthogroup."""
    og_id, aligned_seqs, species_treefile_str = args
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species
    from pathlib import Path as _Path
    from Bio import SeqIO as _SeqIO

    root = get_local_storage_root()

    codon_aln_path = _Path(root) / "meme" / og_id / "codon_aln.fna"
    codon_aln = None
    if codon_aln_path.exists():
        try:
            codon_aln = {r.id: str(r.seq) for r in _SeqIO.parse(str(codon_aln_path), "fasta")}
        except Exception:
            codon_aln = None

    if codon_aln is None:
        cds_seqs: dict[str, str] = {}
        for label in aligned_seqs:
            parts = label.split("|")
            accession = parts[-1] if parts else label
            cds = fetch_cds_for_protein(accession)
            if cds:
                cds_seqs[label] = cds
        if len(cds_seqs) == len(aligned_seqs):
            codon_aln = protein_to_codon_alignment(aligned_seqs, cds_seqs)

    if codon_aln is None:
        return og_id, {}

    species_ids = [label.split("|")[0] for label in codon_aln]
    tree = prune_tree_to_species(_Path(species_treefile_str), species_ids)

    fel_json = run_fel(codon_aln, tree, og_id)
    busted_json = run_busted(codon_aln, tree, og_id)
    fel_result = parse_fel_results(fel_json) if fel_json else {"fel_sites": 0}
    busted_result = parse_busted_results(busted_json) if busted_json else {"busted_pvalue": 1.0}
    return og_id, {**fel_result, **busted_result}


def run_fel_busted_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
    motifs_by_og: dict[str, list],
    species_treefile: Path,
    flush_callback=None,
) -> dict[str, dict]:
    """Run FEL and BUSTED in parallel for all candidate orthogroups.

    Args:
        flush_callback: optional callable(batch: dict) called every 500 OGs
            to write results to DB incrementally (enables spot-safe resume).

    Returns {og_id: {"fel_sites": int, "busted_pvalue": float}}.
    """
    from pipeline.layer2_evolution.selection import _should_run_hyphy

    candidates = [og for og in aligned_orthogroups if _should_run_hyphy(og, motifs_by_og)]

    # Skip OGs already scored — enables resume after spot termination
    if flush_callback:
        already_scored = _get_already_scored_fel_busted_ogs(candidates)
        if already_scored:
            log.info("  FEL+BUSTED: skipping %d already-scored OGs (resuming).", len(already_scored))
            candidates = [og for og in candidates if og not in already_scored]

    log.info("Running FEL + BUSTED on %d candidate orthogroups...", len(candidates))
    if not candidates:
        log.info("  FEL+BUSTED: all OGs already scored.")
        return {}

    n_cpu = os.cpu_count() or 8
    n_workers = max(1, n_cpu - 2)
    total = len(candidates)
    all_results: dict[str, dict] = {}
    pending_flush: dict[str, dict] = {}
    _FLUSH_INTERVAL = 500
    done = 0

    work_items = [
        (og_id, aligned_orthogroups[og_id], str(species_treefile))
        for og_id in candidates
    ]

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_fel_busted_worker, item): item[0] for item in work_items}
        for future in as_completed(futures):
            og_id, result = future.result()
            done += 1
            if result:
                all_results[og_id] = result
                pending_flush[og_id] = result
            if done % 100 == 0 or done == total:
                log.info("  FEL+BUSTED: %d / %d orthogroups (%.0f%%).", done, total, 100 * done / total)
            if flush_callback and len(pending_flush) >= _FLUSH_INTERVAL:
                try:
                    flush_callback(pending_flush)
                    pending_flush.clear()
                except Exception as exc:
                    log.warning("FEL+BUSTED flush failed: %s", exc)

    if flush_callback and pending_flush:
        try:
            flush_callback(pending_flush)
            pending_flush.clear()
        except Exception as exc:
            log.warning("FEL+BUSTED final flush failed: %s", exc)

    log.info("FEL + BUSTED complete: %d / %d orthogroups processed.", len(all_results), len(candidates))
    return {} if flush_callback else all_results


def _get_already_scored_fel_busted_ogs(og_ids: list[str]) -> set[str]:
    """Return set of og_ids that already have fel_sites scored in the DB."""
    try:
        from db.session import get_session
        from db.models import EvolutionScore, Ortholog
        from sqlalchemy import and_
        with get_session() as session:
            scored_genes = {
                row.gene_id for row in session.query(EvolutionScore.gene_id).filter(
                    EvolutionScore.fel_sites.isnot(None)
                )
            }
            og_map = {
                row.orthofinder_og: row.gene_id
                for row in session.query(Ortholog.orthofinder_og, Ortholog.gene_id).filter(
                    Ortholog.orthofinder_og.in_(og_ids)
                )
            }
        return {og for og, gid in og_map.items() if gid in scored_genes}
    except Exception as exc:
        log.debug("Could not check FEL+BUSTED scored OGs: %s", exc)
        return set()



# ---------------------------------------------------------------------------
# RELAX — Branch-specific rate acceleration / relaxation
# ---------------------------------------------------------------------------

_relax_branch_cache: Optional[tuple[set[str], set[str], set[str]]] = None


def _get_relax_branch_sets() -> tuple[set[str], set[str], set[str]]:
    """Return (test_set, reference_set, outgroup_set) derived from the Species DB table.

    Classification is driven entirely by the species registry — no hardcoded
    species names — so this works correctly when the phenotype or species panel
    changes.

    Test      = species with the target phenotype (not baseline / control / outgroup)
    Reference = species with is_control=True OR 'baseline' in phenotypes
    Outgroup  = species with 'outgroup' in phenotypes (folded into Reference for RELAX,
                which requires every branch to be labeled)
    """
    global _relax_branch_cache
    if _relax_branch_cache is not None:
        return _relax_branch_cache

    test: set[str] = set()
    reference: set[str] = set()
    outgroup: set[str] = set()

    try:
        with get_session() as session:
            rows = session.query(Species.id, Species.phenotypes, Species.is_control).all()
        for sp_id, phenotypes, is_ctrl in rows:
            phenotypes = phenotypes or []
            if "outgroup" in phenotypes:
                outgroup.add(sp_id)
            elif is_ctrl or "baseline" in phenotypes:
                reference.add(sp_id)
            else:
                test.add(sp_id)
    except Exception as exc:
        log.warning("Could not load RELAX branch sets from DB: %s", exc)

    _relax_branch_cache = (test, reference, outgroup)
    return _relax_branch_cache


def run_relax(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
    test_branches: Optional[list[str]] = None,
    reference_branches: Optional[list[str]] = None,
) -> Optional[dict]:
    """Run HyPhy RELAX on a codon alignment to detect branch-specific rate shifts.

    RELAX compares the *intensity* of selection (the k parameter) between
    Test branches (trait-selected species) and Reference branches
    (baseline + control + outgroup species).
    k > 1 → selection intensified; k < 1 → relaxation.

    Branch classification is registry-driven: Test = species carrying the target
    phenotype; Reference = baseline / control / outgroup species.  HyPhy requires
    every branch to be labeled, so outgroup species are folded into Reference.

    Args:
        codon_aln: {label: codon_sequence}
        tree_newick: newick tree string
        og_id: orthogroup ID (for output directory)
        test_branches: override Test labels (default: registry-derived)
        reference_branches: override Reference labels (default: registry-derived)

    Returns:
        Raw RELAX JSON dict or None on failure.
    """
    tools = get_tool_config()
    hyphy_bin = tools.get("hyphy_bin", "hyphy")

    root = Path(get_local_storage_root())
    relax_dir = root / "relax" / og_id
    relax_dir.mkdir(parents=True, exist_ok=True)

    aln_path = relax_dir / "codon_aln.fna"
    tree_path = relax_dir / "species.treefile"
    out_path = relax_dir / "relax.json"

    if out_path.exists() and out_path.stat().st_size > 0:
        try:
            with open(out_path) as f:
                return json.load(f)
        except Exception:
            pass

    species_labels = []
    with open(aln_path, "w") as f:
        seen_species = set()
        for label, seq in codon_aln.items():
            species = label.split("|")[0]
            if species in seen_species:
                continue
            seen_species.add(species)
            species_labels.append(species)
            f.write(f">{species}\n{seq}\n")

    test_set, reference_set, outgroup_set = _get_relax_branch_sets()
    if test_branches is None:
        test_branches = [sp for sp in species_labels if sp in test_set]
    if reference_branches is None:
        # Outgroups fold into Reference: HyPhy errors on unlabeled branches.
        reference_branches = [sp for sp in species_labels if sp in reference_set | outgroup_set]
        # Any remaining species not in either set default to Reference (safe fallback).
        labeled = set(test_branches) | set(reference_branches)
        unlabeled = [sp for sp in species_labels if sp not in labeled]
        if unlabeled:
            log.debug("RELAX %s: %d unlabeled species defaulting to Reference: %s",
                      og_id, len(unlabeled), unlabeled)
            reference_branches.extend(unlabeled)

    if len(test_branches) < 1 or len(reference_branches) < 1:
        log.info("RELAX %s: need ≥1 test and ≥1 reference branch, got %d/%d",
                 og_id, len(test_branches), len(reference_branches))
        return None

    labeled_tree = tree_newick
    for branch in sorted(test_branches, key=len, reverse=True):
        labeled_tree = labeled_tree.replace(f"{branch}:", f"{branch}{{Test}}:")
    for branch in sorted(reference_branches, key=len, reverse=True):
        labeled_tree = labeled_tree.replace(f"{branch}:", f"{branch}{{Reference}}:")

    tree_path.write_text(labeled_tree)

    cpus = int(os.environ.get("HYPHY_CPUS", 1))
    hyphy_env = {**os.environ, "CPU": str(cpus)}

    cmd = [
        hyphy_bin, "relax",
        "--alignment", str(aln_path),
        "--tree", str(tree_path),
        "--output", str(out_path),
        "--test", "Test",
        "--reference", "Reference",   # required: without this HyPhy prompts interactively → rc=1
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, check=False, env=hyphy_env)
        if result.returncode != 0 or not out_path.exists():
            log.warning("RELAX failed for %s (rc=%d): stdout=%s stderr=%s",
                        og_id, result.returncode if hasattr(result, 'returncode') else -1,
                        result.stdout[:300], result.stderr[:300])
            return None
        with open(out_path) as f:
            return json.load(f)
    except Exception as exc:
        log.warning("RELAX error for %s: %s", og_id, exc)
        return None


def parse_relax_results(relax_json: dict) -> dict:
    """Extract k (rate intensity) and p-value from RELAX output.

    Returns {"relax_k": float, "relax_pvalue": float}.
    k > 1 means test branches have intensified selection (acceleration).
    """
    try:
        test_results = relax_json.get("test results", {})
        pvalue = test_results.get("p-value")
        k = test_results.get("relaxation or intensification parameter")
        return {
            "relax_k": round(float(k), 4) if k is not None else None,
            "relax_pvalue": round(float(pvalue), 6) if pvalue is not None else None,
        }
    except Exception as exc:
        log.debug("RELAX parse error: %s", exc)
        return {"relax_k": None, "relax_pvalue": None}


def _relax_worker(args: tuple) -> tuple[str, Optional[dict]]:
    """Top-level worker: run RELAX for one orthogroup."""
    og_id, aligned_seqs, species_treefile_str, storage_root = args
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species
    from pathlib import Path as _Path
    from Bio import SeqIO as _SeqIO

    codon_aln_path = _Path(storage_root) / "meme" / og_id / "codon_aln.fna"
    if not codon_aln_path.exists():
        return og_id, None
    try:
        codon_aln = {rec.id: str(rec.seq) for rec in _SeqIO.parse(str(codon_aln_path), "fasta")}
    except Exception:
        return og_id, None
    if not codon_aln:
        return og_id, None

    species_ids = [label.split("|")[0] for label in codon_aln]
    tree = prune_tree_to_species(_Path(species_treefile_str), species_ids)

    relax_json = run_relax(codon_aln, tree, og_id)
    if relax_json is not None:
        return og_id, parse_relax_results(relax_json)
    return og_id, None


def run_relax_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
    motifs_by_og: dict[str, list],
    species_treefile: Path,
    flush_callback=None,
) -> dict[str, dict]:
    """Run RELAX in parallel for all candidate orthogroups that have codon alignments.

    Args:
        flush_callback: optional callable(batch: dict) called every 500 OGs
            to write results to DB incrementally (enables spot-safe resume).

    Returns {og_id: {"relax_k": float, "relax_pvalue": float}}.
    """
    from pipeline.layer2_evolution.selection import _should_run_hyphy

    candidates = [og for og in aligned_orthogroups if _should_run_hyphy(og, motifs_by_og)]

    # Skip OGs already scored — enables resume after spot termination
    if flush_callback:
        already_scored = _get_already_scored_relax_ogs(candidates)
        if already_scored:
            log.info("  RELAX: skipping %d already-scored OGs (resuming).", len(already_scored))
            candidates = [og for og in candidates if og not in already_scored]

    log.info("Running RELAX on %d candidate orthogroups...", len(candidates))
    if not candidates:
        log.info("  RELAX: all OGs already scored.")
        return {}

    n_cpu = os.cpu_count() or 8
    n_workers = max(1, n_cpu - 2)
    total = len(candidates)
    all_results: dict[str, dict] = {}
    pending_flush: dict[str, dict] = {}
    _FLUSH_INTERVAL = 500
    done = 0
    storage_root = str(get_local_storage_root())

    work_items = [
        (og_id, aligned_orthogroups[og_id], str(species_treefile), storage_root)
        for og_id in candidates
    ]

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_relax_worker, item): item[0] for item in work_items}
        for future in as_completed(futures):
            og_id, result = future.result()
            done += 1
            if result is not None:
                all_results[og_id] = result
                pending_flush[og_id] = result
            if done % 100 == 0 or done == total:
                log.info("  RELAX: %d / %d orthogroups (%.0f%%).", done, total, 100 * done / total)
            if flush_callback and len(pending_flush) >= _FLUSH_INTERVAL:
                try:
                    flush_callback(pending_flush)
                    pending_flush.clear()
                except Exception as exc:
                    log.warning("RELAX flush failed: %s", exc)

    if flush_callback and pending_flush:
        try:
            flush_callback(pending_flush)
            pending_flush.clear()
        except Exception as exc:
            log.warning("RELAX final flush failed: %s", exc)

    log.info("RELAX complete: %d / %d orthogroups", len(all_results), len(candidates))
    return {} if flush_callback else all_results


def _get_already_scored_relax_ogs(og_ids: list[str]) -> set[str]:
    """Return set of og_ids that already have relax_k scored in the DB."""
    try:
        from db.session import get_session
        from db.models import EvolutionScore, Ortholog
        with get_session() as session:
            scored_genes = {
                row.gene_id for row in session.query(EvolutionScore.gene_id).filter(
                    EvolutionScore.relax_k.isnot(None)
                )
            }
            og_map = {
                row.orthofinder_og: row.gene_id
                for row in session.query(Ortholog.orthofinder_og, Ortholog.gene_id).filter(
                    Ortholog.orthofinder_og.in_(og_ids)
                )
            }
        return {og for og, gid in og_map.items() if gid in scored_genes}
    except Exception as exc:
        log.debug("Could not check RELAX scored OGs: %s", exc)
        return set()


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def _meme_worker(args: tuple) -> tuple[str, Optional[dict]]:
    """Top-level worker: FUBAR pre-screen → MEME (or aBSREL proxy).

    Two-tier strategy:
      1. Run FUBAR (~10-30s) to detect pervasive selection signal.
      2. Only run MEME (~2-10min) if FUBAR finds ≥1 positively selected site,
         OR if the alignment passes QC but FUBAR is inconclusive.
      3. Fall back to protein aBSREL proxy if CDS unavailable.

    Uses pruned species tree — HyPhy re-optimizes branch lengths internally,
    so the species tree topology (well-resolved with 1000 bootstraps) is
    sufficient and avoids noisy per-gene tree inference on small alignments.
    """
    og_id, aligned_seqs, species_treefile_str, motifs_by_og = args
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species
    from pipeline.layer2_evolution.selection import (
        _should_run_hyphy, parse_absrel_results, run_absrel, write_hyphy_input,
    )

    cds_seqs: dict[str, str] = {}
    for label in aligned_seqs:
        parts = label.split("|")
        accession = parts[-1] if parts else label
        cds = fetch_cds_for_protein(accession)
        if cds:
            cds_seqs[label] = cds

    cds_coverage = len(cds_seqs) / max(len(aligned_seqs), 1)
    if len(cds_seqs) < len(aligned_seqs):
        missing = [lbl for lbl in aligned_seqs if lbl not in cds_seqs]
        log.debug(
            "%s: CDS fetched for %d/%d species (%.0f%%); missing: %s — "
            "building partial codon alignment",
            og_id, len(cds_seqs), len(aligned_seqs), 100 * cds_coverage,
            [lbl.split("|")[0] for lbl in missing],
        )

    # Build codon alignment from whatever species have CDS.
    # protein_to_codon_alignment enforces min_coverage=0.60 and min_seqs=3 internally,
    # so partial inputs are safe — it returns None if thresholds aren't met.
    # Previously required ALL species to have CDS (strict equality), which silently
    # skipped the entire orthogroup whenever any reference species CDS failed.
    codon_aln = None
    if cds_seqs:
        codon_aln = protein_to_codon_alignment(aligned_seqs, cds_seqs)

    if codon_aln is not None:
        if not _codon_aln_passes_qc(codon_aln, og_id):
            return og_id, None

        species_ids = [label.split("|")[0] for label in codon_aln]
        gene_tree = prune_tree_to_species(Path(species_treefile_str), species_ids)

        # FUBAR pre-screen: fast Bayesian scan for any selection signal
        fubar_json = run_fubar(codon_aln, gene_tree, og_id)
        fubar_result = parse_fubar_results(fubar_json) if fubar_json else {"fubar_pos_sites": 0}
        run_full_meme = fubar_result.get("fubar_pos_sites", 0) > 0

        if run_full_meme:
            meme_json = run_meme(codon_aln, gene_tree, og_id)
            if meme_json is not None:
                return og_id, ("meme", parse_meme_results(meme_json, og_id))
        else:
            log.debug("FUBAR found no signal for %s — skipping MEME", og_id)

        # If FUBAR found no signal but we have good codon alignment, still record
        # the FUBAR result as a weak negative
        if not run_full_meme:
            return og_id, ("fubar_only", _meme_null_result(og_id))

    # Fallback: protein-based aBSREL proxy
    species_ids = [label.split("|")[0] for label in aligned_seqs]
    pruned_tree = prune_tree_to_species(Path(species_treefile_str), species_ids)
    aln_path, tree_path = write_hyphy_input(og_id, aligned_seqs, pruned_tree)
    raw_result = run_absrel(aln_path, tree_path, og_id)
    if raw_result is not None:
        return og_id, ("proxy", parse_absrel_results(raw_result))
    return og_id, None


def run_meme_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
    motifs_by_og: dict[str, list],
    species_treefile: Path,
    flush_callback=None,
) -> dict[str, dict]:
    """Run MEME in parallel for all candidate orthogroups.

    Falls back to protein divergence proxy for any orthogroup where CDS
    fetching fails (e.g. UniProt accessions, missing NCBI links).

    Args:
        flush_callback: optional callable(batch: dict) called every 500 OGs
            to write results to DB incrementally.

    Returns {og_id: parsed_selection_result} for any remaining unflushed results.
    """
    from pipeline.layer2_evolution.selection import _should_run_hyphy
    candidates = [og for og in aligned_orthogroups if _should_run_hyphy(og, motifs_by_og)]
    log.info("Running MEME on %d candidate orthogroups (falls back to proxy if CDS unavailable)...",
             len(candidates))

    # Pre-fetch all CDS sequences to disk cache before spawning workers.
    # This avoids 14 workers competing for NCBI's rate limit simultaneously.
    _prefetch_all_cds(aligned_orthogroups)

    n_cpu = os.cpu_count() or 8
    n_workers = max(1, n_cpu - 2)
    total = len(candidates)
    all_results: dict[str, dict] = {}
    meme_success = 0
    fubar_only = 0
    proxy_fallback = 0
    done = 0
    pending_flush: dict[str, dict] = {}
    _FLUSH_INTERVAL = 500

    work_items = [
        (og_id, aligned_orthogroups[og_id], str(species_treefile), motifs_by_og)
        for og_id in candidates
    ]

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_meme_worker, item): item[0] for item in work_items}
        for future in as_completed(futures):
            og_id, result = future.result()
            done += 1
            if result is not None:
                kind, parsed = result
                pending_flush[og_id] = parsed
                all_results[og_id] = parsed
                if kind == "meme":
                    meme_success += 1
                elif kind == "fubar_only":
                    fubar_only += 1
                else:
                    proxy_fallback += 1
            if done % 100 == 0 or done == total:
                log.info("  MEME: %d / %d (%.0f%%) — %d MEME, %d FUBAR-only, %d proxy, %d failed.",
                         done, total, 100 * done / total,
                         meme_success, fubar_only, proxy_fallback,
                         done - meme_success - fubar_only - proxy_fallback)
            if flush_callback and len(pending_flush) >= _FLUSH_INTERVAL:
                try:
                    flush_callback(pending_flush)
                    pending_flush.clear()
                except Exception as exc:
                    log.warning("MEME flush failed: %s", exc)

    # Flush any remaining
    if flush_callback and pending_flush:
        try:
            flush_callback(pending_flush)
            pending_flush.clear()
        except Exception as exc:
            log.warning("MEME final flush failed: %s", exc)

    log.info("MEME complete: %d MEME, %d FUBAR-only, %d proxy, %d failed.",
             meme_success, fubar_only, proxy_fallback,
             total - meme_success - fubar_only - proxy_fallback)

    # Apply Benjamini-Hochberg FDR correction across all site-level p-values
    from pipeline.stats import apply_bh_to_meme_results
    all_results = apply_bh_to_meme_results(all_results)

    # Return only unflushed results (empty if flush_callback handled everything)
    return {} if flush_callback else all_results
