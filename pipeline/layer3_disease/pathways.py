"""Fetch GO terms and pathway IDs from UniProt for genes."""

import logging
import time
from typing import Optional

import requests

from db.models import Gene
from db.session import get_session

log = logging.getLogger(__name__)

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"


def fetch_go_and_pathways(uniprot_id: str) -> tuple[list[str], list[str]]:
    """Fetch GO terms and pathway IDs for a UniProt accession.

    Returns (go_terms, pathway_ids). Uses UniProt REST API.
    """
    go_terms = []
    pathway_ids = []
    url = f"{UNIPROT_API}/{uniprot_id}.json"
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        log.warning("UniProt fetch failed for %s: %s", uniprot_id, e)
        return go_terms, pathway_ids

    # Cross-references: GO (id is GO:xxxxx) and Reactome
    for xref in data.get("uniProtKBCrossReferences", []):
        db = xref.get("database", "")
        xref_id = xref.get("id", "")
        if db == "GO" and xref_id:
            go_terms.append(xref_id)
        elif db == "Reactome" and xref_id:
            pathway_ids.append(xref_id)

    # Deduplicate and filter empty
    go_terms = list(dict.fromkeys(g for g in go_terms if g))
    pathway_ids = list(dict.fromkeys(p for p in pathway_ids if p))
    return go_terms, pathway_ids


def annotate_genes_pathways(gene_ids: list[str]) -> int:
    """Fetch GO/pathway data for each gene with human_protein and update Gene rows."""
    updated = 0
    with get_session() as session:
        for gene_id in gene_ids:
            gene = session.get(Gene, gene_id)
            if not gene or not gene.human_protein:
                continue
            go_terms, pathway_ids = fetch_go_and_pathways(gene.human_protein)
            gene.go_terms = go_terms
            gene.pathway_ids = pathway_ids
            updated += 1
            time.sleep(0.2)
    return updated
