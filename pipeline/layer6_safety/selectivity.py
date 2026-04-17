"""Selectivity — protein family size and essentiality (gnomAD pLI) for safety.

Protein family size is determined via UniProt: we first look up the gene to get
its protein family classification, then count how many reviewed human proteins
share that family. A larger family = more off-target paralogue risk.
"""

import logging
import re
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene, SafetyFlag
from db.session import get_session

log = logging.getLogger(__name__)

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
ESSENTIAL_PLI_THRESHOLD = 0.9


def _hgnc_symbol(raw: str) -> str:
    """Convert OMA entry name (AKT1_HUMAN) or 'human|AKT1_HUMAN' to bare HGNC symbol (AKT1)."""
    if not raw:
        return ""
    # Strip 'human|' prefix if present
    entry = raw.split("|")[-1].strip()
    # Strip '_HUMAN' suffix if present
    if "_" in entry:
        return entry.split("_")[0]
    return entry


def fetch_family_size(gene_symbol: str) -> Optional[int]:
    """Fetch number of reviewed human proteins in the same protein family via UniProt.

    Two-step:
      1. Look up the gene entry to get its protein_families classification text.
      2. Count reviewed human proteins sharing that family name.

    Returns 0 if gene has no family classification, None on API error.
    """
    try:
        # Step 1: get the UniProt entry for this gene
        r = requests.get(
            UNIPROT_SEARCH,
            params={
                "query": f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
                "fields": "accession,protein_families",
                "format": "json",
                "size": 1,
            },
            timeout=15,
        )
        if r.status_code != 200:
            log.warning("UniProt entry lookup %s: HTTP %d", gene_symbol, r.status_code)
            return None

        results = r.json().get("results", [])
        if not results:
            log.info("UniProt: no reviewed entry found for %s", gene_symbol)
            return 0

        # Extract protein family text (may be list of strings or absent)
        entry = results[0]
        families = entry.get("proteinDescription", {}).get("proteinFamilies", [])
        if not families:
            # Try alternate path
            comments = entry.get("comments", [])
            families = [
                c.get("texts", [{}])[0].get("value", "")
                for c in comments
                if c.get("commentType") == "SIMILARITY"
            ]

        if not families:
            log.info("UniProt %s: no protein family classification found.", gene_symbol)
            return 0

        family_text = families[0] if isinstance(families[0], str) else str(families[0])
        if not family_text:
            return 0

        # Extract the MOST SPECIFIC subfamily from the hierarchical SIMILARITY text.
        # UniProt SIMILARITY comments are period-delimited from broad to specific, e.g.:
        #   "Belongs to the protein kinase superfamily. AGC Ser/Thr protein kinase family. RAC subfamily."
        # Querying the broadest level ("protein kinase superfamily") yields hundreds of
        # members and overestimates off-target risk. We want the tightest classification
        # ("RAC subfamily" → AKT1, AKT2, AKT3) to give a realistic paralog count.
        segments = [s.strip() for s in family_text.split(".") if s.strip()]
        # Strip "Belongs to (the)..." prefix from the leading segment
        if segments:
            segments[0] = re.sub(
                r"^Belongs to (the )?", "", segments[0], flags=re.IGNORECASE
            ).strip()
        segments = [s for s in segments if s]

        if segments:
            # Use the most specific (last) non-empty segment
            cleaned = segments[-1].replace("\"", "").strip()
        else:
            cleaned = ""

        if not cleaned:
            cleaned = family_text.replace("\"", "")[:80]

        # Step 2: count reviewed human proteins in the same family.
        # UniProt query field for family classification is `family:`.
        r2 = requests.get(
            UNIPROT_SEARCH,
            params={
                "query": f'family:"{cleaned}" AND organism_id:9606 AND reviewed:true',
                "fields": "accession",
                "format": "json",
                "size": 500,
            },
            timeout=15,
        )
        if r2.status_code != 200:
            # Fallback: if family query is not accepted, do not fail the whole step.
            log.warning("UniProt family search %s: HTTP %d (query=%s)", gene_symbol, r2.status_code, cleaned)
            return 1

        members = r2.json().get("results", [])
        count = len(members)
        if count == 0:
            # Keep a conservative non-null value when family is known but query is sparse.
            count = 1
        log.info("UniProt %s: family='%s'  members=%d", gene_symbol, cleaned[:60], count)
        return count

    except Exception as exc:
        log.warning("UniProt family %s: %s", gene_symbol, exc)
        return None


def annotate_genes_selectivity(gene_ids: list[str]) -> int:
    """Populate SafetyFlag.family_size and is_essential (from DiseaseAnnotation.gnomad_pli)."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            symbol = _hgnc_symbol(gene.gene_symbol)
            log.info("Selectivity querying: %s (raw: %s)", symbol, gene.gene_symbol)

            family_size = fetch_family_size(symbol)
            ann = session.get(DiseaseAnnotation, gid)
            pli = ann.gnomad_pli if ann else None
            is_essential = pli is not None and pli >= ESSENTIAL_PLI_THRESHOLD

            sf = session.get(SafetyFlag, gid)
            if sf is None:
                sf = SafetyFlag(gene_id=gid)
                session.add(sf)
            if family_size is not None:
                sf.family_size = family_size
            sf.is_essential = is_essential
            log.info("Selectivity %s: family_size=%s  pLI=%s  is_essential=%s",
                     symbol, family_size, pli, is_essential)
            updated += 1
        session.commit()
    log.info("Selectivity: updated %d genes.", updated)
    return updated
