"""Step 4b-i — Pfam/InterPro domain annotation for divergent motifs.

For each gene's human protein (UniProt accession in Gene.human_protein):
  1. Fetch domain annotations (positions, names) via UniProt REST API (batch) or
     InterPro REST API (per-protein fallback).
  2. For each DivergentMotif, check if its [start_pos, end_pos] overlaps any domain.
  3. Set domain_name and in_functional_domain on the motif.

UniProt strategy (default, _STRATEGY = "uniprot"):
  - Sends batches of up to _UNIPROT_BATCH_SIZE accessions per HTTP request.
  - ~12 requests for 2 400 genes vs. 2 400 requests with InterPro — ~200× faster.
  - Returns UniProt curated Domain + Repeat features with exact residue positions.
  - Results are cached to disk ({storage_root}/pfam/{accession}.json).

InterPro strategy (fallback, _STRATEGY = "interpro"):
  - One request per protein, 50 concurrent threads.
  - Slower but returns richer cross-database annotations (Pfam, PRINTS, SMART …).

Domain annotation is purely additive — it does not change divergence scores,
only labels motifs to allow downstream filtering for functional significance.
"""

import json
import logging
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import requests

from db.models import DivergentMotif, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_local_storage_root, get_tool_config

log = logging.getLogger(__name__)

# --------------------------------------------------------------------------- #
# Strategy selection
# --------------------------------------------------------------------------- #
# "uniprot"  — fast batch REST queries (~200 accessions / request, ~12 total)
# "interpro" — original per-protein InterPro queries (50 threads)
_STRATEGY: str = "uniprot"

# --------------------------------------------------------------------------- #
# Shared constants
# --------------------------------------------------------------------------- #
REQUEST_TIMEOUT = 30
RETRY_DELAY = 2.0

# --------------------------------------------------------------------------- #
# UniProt strategy constants
# --------------------------------------------------------------------------- #
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
# UniProt accepts up to 200 accessions per request; this is a protocol limit,
# not a tuning parameter — do not raise above 200.
_UNIPROT_BATCH_SIZE = 200
_UNIPROT_FIELDS = "accession,ft_domain,ft_repeat"
# Feature types from UniProt that correspond to functional domains
_UNIPROT_DOMAIN_TYPES = {"domain", "repeat"}


def _uniprot_batch_workers() -> int:
    """Return number of parallel UniProt batch workers from config (default 4)."""
    return int(get_tool_config().get("uniprot_batch_workers", 4))


# --------------------------------------------------------------------------- #
# InterPro strategy constants
# --------------------------------------------------------------------------- #
INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"


def _interpro_fetch_workers() -> int:
    """Return number of parallel InterPro workers from config (default 50)."""
    return int(get_tool_config().get("interpro_fetch_workers", 50))


# InterPro entry types treated as functional regions
_INTERPRO_DOMAIN_TYPES = {"domain", "family", "homologous_superfamily", "repeat"}


# --------------------------------------------------------------------------- #
# Disk cache helpers (shared by both strategies)
# --------------------------------------------------------------------------- #

def _cache_path(accession: str) -> Path:
    root = Path(get_local_storage_root()) / "pfam"
    root.mkdir(parents=True, exist_ok=True)
    return root / f"{accession}.json"


def _load_cache(accession: str) -> Optional[list[dict]]:
    path = _cache_path(accession)
    if path.exists():
        try:
            return json.loads(path.read_text())
        except Exception:
            pass
    return None


def _write_cache(accession: str, domains: list[dict]) -> None:
    try:
        _cache_path(accession).write_text(json.dumps(domains))
    except Exception as exc:
        log.debug("Cache write failed for %s: %s", accession, exc)


# --------------------------------------------------------------------------- #
# UniProt batch strategy
# --------------------------------------------------------------------------- #

def _parse_uniprot_features(entry: dict) -> list[dict]:
    """Extract domain/repeat features with positions from a UniProt entry dict."""
    domains: list[dict] = []
    for feature in entry.get("features", []):
        ftype = feature.get("type", "").lower()
        if ftype not in _UNIPROT_DOMAIN_TYPES:
            continue
        loc = feature.get("location", {})
        start_info = loc.get("start", {})
        end_info = loc.get("end", {})
        start = start_info.get("value")
        end = end_info.get("value")
        if start is None or end is None:
            continue
        name = feature.get("description", "") or ftype.capitalize()
        domains.append({
            "name": name,
            "start": int(start),
            "end": int(end),
            "database": "uniprot",
        })
    return domains


def _fetch_uniprot_batch(accessions: list[str]) -> dict[str, list[dict]]:
    """Fetch domains for a batch of UniProt accessions in a single POST request.

    Uses POST to avoid URL-length limits that cause HTTP 400 with large batches.
    Returns {accession: [domain, ...]} for every accession in the batch.
    Accessions not found in UniProt map to an empty list.
    """
    query = " OR ".join(f"accession:{acc}" for acc in accessions)
    data = {
        "query": query,
        "fields": _UNIPROT_FIELDS,
        "format": "json",
        "size": len(accessions),
    }

    for attempt in range(4):
        try:
            r = requests.post(UNIPROT_API, data=data, timeout=REQUEST_TIMEOUT)
            if r.status_code == 200:
                break
            if r.status_code == 429:
                wait = RETRY_DELAY * (attempt + 2) * 2
                log.debug("UniProt rate limit (attempt %d), sleeping %.1fs", attempt + 1, wait)
                time.sleep(wait)
            elif r.status_code >= 500:
                time.sleep(RETRY_DELAY * (attempt + 1))
            else:
                log.warning("UniProt returned HTTP %d for batch of %d — response: %s",
                            r.status_code, len(accessions), r.text[:200])
                break
        except requests.RequestException as exc:
            log.debug("UniProt request error (attempt %d): %s", attempt + 1, exc)
            time.sleep(RETRY_DELAY * (attempt + 1))
    else:
        return {acc: [] for acc in accessions}

    result: dict[str, list[dict]] = {acc: [] for acc in accessions}
    try:
        for entry in r.json().get("results", []):
            acc = entry.get("primaryAccession", "").upper()
            if acc:
                result[acc] = _parse_uniprot_features(entry)
    except Exception as exc:
        log.warning("UniProt batch parse error: %s", exc)

    return result


def fetch_domains_uniprot(accession_list: list[str]) -> dict[str, list[dict]]:
    """Fetch domains for all accessions using UniProt batch API.

    Splits the list into chunks of _UNIPROT_BATCH_SIZE, sends one HTTP
    request per chunk, and merges results.  Cache hits are returned
    immediately without a network call.

    Returns {accession: [domain, ...]} for every input accession.
    """
    domain_map: dict[str, list[dict]] = {}
    to_fetch: list[str] = []

    for acc in accession_list:
        cached = _load_cache(acc)
        if cached is not None:
            domain_map[acc] = cached
        else:
            to_fetch.append(acc)

    if not to_fetch:
        return domain_map

    total = len(to_fetch)
    batches = [to_fetch[i:i + _UNIPROT_BATCH_SIZE] for i in range(0, total, _UNIPROT_BATCH_SIZE)]
    log.info("  UniProt: fetching %d accessions in %d batches of up to %d ...",
             total, len(batches), _UNIPROT_BATCH_SIZE)

    fetched = 0
    for batch in batches:
        batch_upper = [a.upper() for a in batch]
        batch_result = _fetch_uniprot_batch(batch_upper)
        for acc, domains in batch_result.items():
            domain_map[acc] = domains
            _write_cache(acc, domains)
        fetched += len(batch)
        if fetched % 1000 == 0 or fetched == total:
            log.info("  UniProt fetched %d / %d accessions (%.0f%%).",
                     fetched, total, 100 * fetched / total)

    return domain_map


# --------------------------------------------------------------------------- #
# InterPro per-protein strategy (original implementation, kept as fallback)
# --------------------------------------------------------------------------- #

def fetch_domains_for_protein(accession: str) -> list[dict]:
    """Fetch Pfam + InterPro domain annotations for a single UniProt accession.

    Returns a list of domain dicts:
        [{"name": str, "start": int, "end": int, "database": str}, ...]

    Start/end are 1-indexed residue positions (inclusive).
    Results are cached to disk.
    """
    cached = _load_cache(accession)
    if cached is not None:
        return cached

    domains: list[dict] = []
    url = f"{INTERPRO_API}/entry/all/protein/UniProt/{accession}/?page_size=200"

    r = None
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=REQUEST_TIMEOUT)
            if r.status_code == 404:
                _write_cache(accession, [])
                return []
            if r.status_code == 200:
                break
            if r.status_code == 429:
                time.sleep(RETRY_DELAY * (attempt + 2) * 2)
            else:
                time.sleep(RETRY_DELAY * (attempt + 1))
        except requests.RequestException as exc:
            log.debug("InterPro request failed (attempt %d): %s", attempt + 1, exc)
            time.sleep(RETRY_DELAY)
    else:
        return domains

    if r is None or r.status_code != 200:
        return domains

    try:
        data = r.json()
        for entry in data.get("results", []):
            entry_meta = entry.get("metadata", {})
            db = entry_meta.get("source_database", "")
            name = entry_meta.get("name", "") or entry_meta.get("accession", "")
            entry_type = entry_meta.get("type", "")

            if entry_type.lower() not in _INTERPRO_DOMAIN_TYPES:
                continue

            for prot in entry.get("proteins", []):
                if prot.get("accession", "").upper() != accession.upper():
                    continue
                for location in prot.get("entry_protein_locations", []):
                    for fragment in location.get("fragments", []):
                        start = fragment.get("start")
                        end_ = fragment.get("end")
                        if start is not None and end_ is not None:
                            domains.append({
                                "name": name,
                                "start": int(start),
                                "end": int(end_),
                                "database": db,
                            })
    except Exception as exc:
        log.debug("InterPro parse error for %s: %s", accession, exc)

    seen: set[tuple] = set()
    unique: list[dict] = []
    for d in domains:
        key = (d["name"], d["start"], d["end"])
        if key not in seen:
            seen.add(key)
            unique.append(d)
    unique.sort(key=lambda d: d["start"])

    _write_cache(accession, unique)
    log.debug("  %s: %d domains from InterPro", accession, len(unique))
    return unique


def _fetch_worker(args: tuple) -> tuple[str, str, list[dict]]:
    """Fetch domains for one gene via InterPro (thread-pool worker)."""
    gene_id, accession = args
    domains = fetch_domains_for_protein(accession)
    return gene_id, accession, domains


def fetch_domains_interpro(fetch_args: list[tuple[str, str]]) -> dict[str, list[dict]]:
    """Fetch domains for all (gene_id, accession) pairs using InterPro, one call per protein.

    Returns {gene_id: [domain, ...]} for every input pair.
    """
    domain_map: dict[str, list[dict]] = {}
    done = 0
    total = len(fetch_args)

    with ThreadPoolExecutor(max_workers=_interpro_fetch_workers()) as pool:
        futures = {pool.submit(_fetch_worker, args): args for args in fetch_args}
        for future in as_completed(futures):
            try:
                gene_id, _accession, domains = future.result()
                domain_map[gene_id] = domains
            except Exception as exc:
                gene_id = futures[future][0]
                log.warning("Domain fetch failed for gene %s: %s", gene_id, exc)
                domain_map[gene_id] = []
            done += 1
            if done % 500 == 0 or done == total:
                log.info("  InterPro fetched %d / %d gene domain annotations (%.0f%%).",
                         done, total, 100 * done / total)

    return domain_map


# --------------------------------------------------------------------------- #
# Shared overlap check
# --------------------------------------------------------------------------- #

def _motif_overlaps_domain(motif_start: int, motif_end: int, domain: dict) -> bool:
    """Return True if the motif position range overlaps the domain range.

    Motif positions are 0-indexed (from MAFFT alignment position).
    Domain positions are 1-indexed residue positions (from InterPro/UniProt).
    """
    m_start = motif_start + 1  # convert to 1-indexed
    m_end = motif_end + 1
    return not (m_end < domain["start"] or m_start > domain["end"])


# --------------------------------------------------------------------------- #
# Main annotation entry point
# --------------------------------------------------------------------------- #

def annotate_motif_domains(gene_ids: Optional[list[str]] = None) -> int:
    """Annotate DivergentMotif rows with domain information.

    1. Load genes that have divergent motifs (skip the rest).
    2. Fetch domain positions using the configured _STRATEGY.
    3. Apply overlaps and write to DB in a single session commit.

    Returns number of motifs annotated as in_functional_domain=True.
    """
    with get_session() as session:
        q = session.query(Gene)
        if gene_ids:
            q = q.filter(Gene.id.in_(gene_ids))
        all_genes = q.all()

        genes_with_motifs = {
            row[0] for row in
            session.query(Ortholog.gene_id)
            .join(DivergentMotif, DivergentMotif.ortholog_id == Ortholog.id)
            .distinct()
            .all()
        }
        genes = [g for g in all_genes if g.id in genes_with_motifs]

    log.info("Annotating domains for %d genes with motifs (strategy=%s)...",
             len(genes), _STRATEGY)

    # Build accession list — extract the bare UniProt accession from whatever
    # format gene.human_protein was stored in (reheadered FASTA gives either
    # "human|P12345" or "human|BRCA1_HUMAN" depending on the source FASTA).
    def _clean_accession(raw: str) -> str:
        parts = raw.strip().split("|")
        # "sp|P12345|BRCA1_HUMAN" → index 1 is the accession
        if len(parts) >= 3:
            candidate = parts[1].split("-")[0].strip()
            if re.match(r"^[A-Z][A-Z0-9]{4,9}$", candidate) and "_" not in candidate:
                return candidate.upper()
        # Try each part: first one with no underscore and 5-10 chars is the accession
        for part in parts:
            candidate = part.split("-")[0].strip()
            if re.match(r"^[A-Z][A-Z0-9]{4,9}$", candidate) and "_" not in candidate:
                return candidate.upper()
        # Last resort: strip isoform suffix from the whole value
        return parts[-1].split("-")[0].strip().upper()

    fetch_args: list[tuple[str, str]] = []
    for gene in genes:
        if not gene.human_protein:
            continue
        acc = _clean_accession(gene.human_protein)
        fetch_args.append((gene.id, acc))

    # ------------------------------------------------------------------ #
    # Fetch domains
    # ------------------------------------------------------------------ #
    domain_map: dict[str, list[dict]]   # keyed by gene_id

    if _STRATEGY == "uniprot":
        # Build accession → [gene_ids] index so we can map results back
        acc_to_gene_ids: dict[str, list[str]] = {}
        for gene_id, acc in fetch_args:
            acc_to_gene_ids.setdefault(acc, []).append(gene_id)

        unique_accessions = list(acc_to_gene_ids.keys())
        acc_domains = fetch_domains_uniprot(unique_accessions)  # {acc: [domain,...]}

        # Map back to gene_ids (multiple genes can share an accession)
        domain_map = {}
        for acc, doms in acc_domains.items():
            for gid in acc_to_gene_ids.get(acc, []):
                domain_map[gid] = doms

    else:
        # InterPro: domain_map is already keyed by gene_id
        domain_map = fetch_domains_interpro(fetch_args)

    log.info("All domain annotations fetched. Applying to motifs and writing to DB...")

    # ------------------------------------------------------------------ #
    # Apply overlaps in a single DB session
    # ------------------------------------------------------------------ #
    functional_count = 0
    with get_session() as session:
        for gene in genes:
            domains = domain_map.get(gene.id, [])
            if not domains:
                continue

            motifs = (
                session.query(DivergentMotif)
                .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
                .filter(Ortholog.gene_id == gene.id)
                .all()
            )

            for motif in motifs:
                best_domain = None
                for domain in domains:
                    if _motif_overlaps_domain(motif.start_pos, motif.end_pos, domain):
                        # Prefer UniProt curated > Pfam > other databases
                        if best_domain is None:
                            best_domain = domain
                        elif domain["database"].lower() in ("uniprot", "pfam"):
                            best_domain = domain

                if best_domain:
                    motif.domain_name = best_domain["name"]
                    motif.in_functional_domain = True
                    functional_count += 1
                else:
                    motif.in_functional_domain = False

        session.commit()

    log.info("Domain annotation complete: %d motifs in functional domains.", functional_count)
    return functional_count


def run_pfam_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step4b."""
    return annotate_motif_domains(gene_ids)
