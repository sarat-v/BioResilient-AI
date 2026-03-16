"""Step 4b-i — Pfam/InterPro domain annotation for divergent motifs.

For each gene's human protein (UniProt accession in Gene.human_protein):
  1. Fetch Pfam + InterPro domain annotations from InterPro REST API.
  2. For each DivergentMotif, check if its [start_pos, end_pos] overlaps any domain.
  3. Set domain_name and in_functional_domain on the motif.

Results are cached to disk ({storage_root}/pfam/{accession}.json) to avoid
re-fetching on re-runs.

Domain annotation is purely additive — it does not change divergence scores,
only labels motifs to allow downstream filtering for functional significance.
"""

import json
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import requests

from db.models import DivergentMotif, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_storage_root

log = logging.getLogger(__name__)

INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"
REQUEST_TIMEOUT = 20
RETRY_DELAY = 2.0
# InterPro tolerates ~50 concurrent connections before rate limiting
_FETCH_WORKERS = 50


def _cache_path(accession: str) -> Path:
    root = Path(get_storage_root()) / "pfam"
    root.mkdir(parents=True, exist_ok=True)
    return root / f"{accession}.json"


def fetch_domains_for_protein(accession: str) -> list[dict]:
    """Fetch Pfam + InterPro domain annotations for a UniProt accession.

    Returns a list of domain dicts:
        [{"name": str, "start": int, "end": int, "database": str}, ...]

    Start/end are 1-indexed residue positions (inclusive).
    Results are cached to disk.
    """
    cache = _cache_path(accession)
    if cache.exists():
        try:
            return json.loads(cache.read_text())
        except Exception:
            pass

    domains: list[dict] = []
    url = f"{INTERPRO_API}/entry/all/protein/UniProt/{accession}/?page_size=200"

    for attempt in range(3):
        try:
            r = requests.get(url, timeout=REQUEST_TIMEOUT)
            if r.status_code == 404:
                cache.write_text(json.dumps([]))
                return []
            if r.status_code == 200:
                break
            if r.status_code == 429:
                # Rate limited — back off longer
                time.sleep(RETRY_DELAY * (attempt + 2) * 2)
            else:
                time.sleep(RETRY_DELAY * (attempt + 1))
        except requests.RequestException as exc:
            log.debug("InterPro request failed (attempt %d): %s", attempt + 1, exc)
            time.sleep(RETRY_DELAY)
    else:
        return domains

    try:
        data = r.json()
        results = data.get("results", [])
        for entry in results:
            entry_meta = entry.get("metadata", {})
            db = entry_meta.get("source_database", "")
            name = entry_meta.get("name", "") or entry_meta.get("accession", "")
            entry_type = entry_meta.get("type", "")

            if entry_type.lower() not in ("domain", "family", "homologous_superfamily", "repeat"):
                continue

            proteins = entry.get("proteins", [])
            for prot in proteins:
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

    # Deduplicate and sort
    seen = set()
    unique = []
    for d in domains:
        key = (d["name"], d["start"], d["end"])
        if key not in seen:
            seen.add(key)
            unique.append(d)
    unique.sort(key=lambda d: d["start"])

    cache.write_text(json.dumps(unique))
    log.debug("  %s: %d domains from InterPro", accession, len(unique))
    return unique


def _motif_overlaps_domain(motif_start: int, motif_end: int, domain: dict) -> bool:
    """Return True if the motif position range overlaps the domain range.

    Motif positions are 0-indexed (from MAFFT alignment position).
    Domain positions are 1-indexed residue positions (from InterPro).
    Convert motif to 1-indexed for comparison.
    """
    m_start = motif_start + 1  # convert to 1-indexed
    m_end = motif_end + 1
    d_start = domain["start"]
    d_end = domain["end"]
    return not (m_end < d_start or m_start > d_end)


def _fetch_worker(args: tuple) -> tuple[str, str, list[dict]]:
    """Fetch domains for one gene in a thread. Returns (gene_id, accession, domains)."""
    gene_id, accession = args
    domains = fetch_domains_for_protein(accession)
    return gene_id, accession, domains


def annotate_motif_domains(gene_ids: Optional[list[str]] = None) -> int:
    """Annotate DivergentMotif rows with Pfam/InterPro domain information.

    Fetches InterPro annotations in parallel (thread pool), then applies
    domain overlaps and writes to DB in a single session commit.

    Returns number of motifs annotated as in_functional_domain=True.
    """
    with get_session() as session:
        q = session.query(Gene)
        if gene_ids:
            q = q.filter(Gene.id.in_(gene_ids))
        genes = q.all()

        # Only process genes that actually have divergent motifs — skip the rest
        genes_with_motifs = {
            row[0] for row in
            session.query(Ortholog.gene_id)
            .join(DivergentMotif, DivergentMotif.ortholog_id == Ortholog.id)
            .distinct()
            .all()
        }
        genes = [g for g in genes if g.id in genes_with_motifs]

    log.info("Annotating Pfam domains for %d genes with motifs (%d workers)...",
             len(genes), _FETCH_WORKERS)

    # Build list of (gene_id, accession) to fetch
    fetch_args = []
    for gene in genes:
        accession = gene.human_protein
        if not accession:
            continue
        if "|" in accession:
            accession = accession.split("|")[-1]
        fetch_args.append((gene.id, accession))

    # Parallel fetch — all I/O bound, safe to thread
    domain_map: dict[str, list[dict]] = {}
    done = 0
    total = len(fetch_args)

    with ThreadPoolExecutor(max_workers=_FETCH_WORKERS) as pool:
        futures = {pool.submit(_fetch_worker, args): args for args in fetch_args}
        for future in as_completed(futures):
            try:
                gene_id, accession, domains = future.result()
                domain_map[gene_id] = domains
            except Exception as exc:
                gene_id = futures[future][0]
                log.warning("Domain fetch failed for gene %s: %s", gene_id, exc)
                domain_map[gene_id] = []
            done += 1
            if done % 500 == 0 or done == total:
                log.info("  Fetched %d / %d gene domain annotations (%.0f%%).",
                         done, total, 100 * done / total)

    log.info("All domain annotations fetched. Applying to motifs and writing to DB...")

    # Single session for all DB writes
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
                        if best_domain is None or domain["database"].lower() == "pfam":
                            best_domain = domain

                if best_domain:
                    motif.domain_name = best_domain["name"]
                    motif.in_functional_domain = True
                    functional_count += 1
                else:
                    motif.in_functional_domain = False

        session.commit()

    log.info("Pfam domain annotation complete: %d motifs in functional domains.", functional_count)
    return functional_count


def run_pfam_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step4b."""
    return annotate_motif_domains(gene_ids)
