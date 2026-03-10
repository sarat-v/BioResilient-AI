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
    # Query all entries (Pfam, PANTHER, TIGRFAM, etc.) for this protein
    url = f"{INTERPRO_API}/entry/all/protein/UniProt/{accession}/?page_size=200"

    for attempt in range(3):
        try:
            r = requests.get(url, timeout=REQUEST_TIMEOUT)
            if r.status_code == 404:
                # No domains found — cache empty list
                cache.write_text(json.dumps([]))
                return []
            if r.status_code == 200:
                break
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

            # Only keep domain-type entries (skip sites, homologous superfamilies, etc.)
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
                        if start and end_:
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
    # Overlap: not (m_end < d_start or m_start > d_end)
    return not (m_end < d_start or m_start > d_end)


def annotate_motif_domains(gene_ids: Optional[list[str]] = None) -> int:
    """Annotate DivergentMotif rows with Pfam/InterPro domain information.

    For each gene (optionally restricted to gene_ids), fetches domain
    annotations for the human protein accession, then marks each DivergentMotif
    with domain_name and in_functional_domain.

    Returns number of motifs annotated as in_functional_domain=True.
    """
    functional_count = 0

    with get_session() as session:
        q = session.query(Gene)
        if gene_ids:
            q = q.filter(Gene.id.in_(gene_ids))
        genes = q.all()

        log.info("Annotating Pfam domains for %d genes...", len(genes))

        for gene in genes:
            accession = gene.human_protein
            if not accession:
                continue

            # Strip species prefix if present
            if "|" in accession:
                accession = accession.split("|")[-1]

            domains = fetch_domains_for_protein(accession)
            if not domains:
                continue

            # Get all motifs for this gene (via ortholog join)
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
                        # Prefer Pfam domain over others; take first match
                        if best_domain is None or domain["database"].lower() == "pfam":
                            best_domain = domain

                if best_domain:
                    motif.domain_name = best_domain["name"]
                    motif.in_functional_domain = True
                    functional_count += 1
                else:
                    motif.in_functional_domain = False

        session.commit()

    log.info(
        "Pfam domain annotation complete: %d motifs in functional domains.",
        functional_count,
    )
    return functional_count


def run_pfam_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step4b."""
    return annotate_motif_domains(gene_ids)
