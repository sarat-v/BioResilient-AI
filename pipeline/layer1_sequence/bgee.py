"""Step 8b — Bgee cross-species expression supplement.

Bgee (https://www.bgee.org) provides curated, normalised gene expression data
across species and tissues, enabling reliable cross-species comparisons.

This module supplements the existing GEO/DESeq2 expression pipeline
(Step 8) with Bgee calls, especially for species with sparse GEO data.

For each Tier1/Tier2 gene's human ortholog:
  1. Resolve the human gene symbol to a Bgee gene ID.
  2. For each resilient species in the panel, look up expression calls in
     tissues relevant to the trait (e.g. liver, kidney for cancer resistance).
  3. If a "present" expression call is found for ≥2 resilient species, store
     an ExpressionResult row with geo_accession="bgee:{bgee_gene_id}".

Bgee REST API docs: https://www.bgee.org/support/api

Note: Bgee has inconsistent coverage for non-model organisms. Results are
stored alongside GEO rows in ExpressionResult and interpreted the same way.
Bgee data is tagged with a "bgee:" prefix on geo_accession for traceability.
"""

import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

import requests

from db.models import CandidateScore, ExpressionResult, Gene, Ortholog, Species
from db.session import get_session

log = logging.getLogger(__name__)

BGEE_API = "https://www.bgee.org/api"
REQUEST_TIMEOUT = 20
_BGEE_WORKERS = 10   # Bgee rate-limits at ~10–15 concurrent connections

# Tissues relevant for each trait (mapped to Uberon or Bgee tissue names)
TRAIT_TISSUES: dict[str, list[str]] = {
    "cancer_resistance": ["liver", "kidney", "skin", "lung", "colon"],
    "longevity":         ["brain", "heart", "liver", "muscle"],
    "viral_resistance":  ["lung", "spleen", "lymph node"],
    "regeneration":      ["limb", "heart", "spinal cord"],
    "default":           ["liver", "brain", "heart", "kidney"],
}

# Bgee species IDs for our panel (NCBI taxon IDs work directly in Bgee)
BGEE_SUPPORTED_TAXIDS = {
    9606:   "human",
    10181:  "naked_mole_rat",
    9785:   "african_elephant",
    59463:  "little_brown_bat",
    13146:  "budgerigar",
    8296:   "axolotl",
    30608:  "mouse_lemur",
}


def _bgee_gene_search(gene_symbol: str, species_taxid: int = 9606) -> Optional[str]:
    """Search Bgee for a gene by symbol and taxon. Returns bgee gene ID or None."""
    try:
        r = requests.get(
            f"{BGEE_API}/gene",
            params={"gene_id": gene_symbol, "species_id": species_taxid, "display_type": "json"},
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        genes = data.get("data", {}).get("genes", [])
        if genes:
            return genes[0].get("geneId")
    except Exception as exc:
        log.debug("Bgee gene search failed for %s: %s", gene_symbol, exc)
    return None


def _bgee_expression_calls(bgee_gene_id: str, species_taxid: int) -> list[dict]:
    """Fetch expression calls for a Bgee gene ID in a species.

    Returns a list of {tissue, expression_level, call_type, quality}.
    call_type is one of: "present", "absent".
    """
    try:
        r = requests.get(
            f"{BGEE_API}/expression",
            params={
                "gene_id": bgee_gene_id,
                "species_id": species_taxid,
                "display_type": "json",
                "data_type": "RNA-Seq",
                "condition_parameters": "anat_entity",
            },
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            return []
        data = r.json()
        calls = data.get("data", {}).get("expressionCalls", [])
        return calls
    except Exception as exc:
        log.debug("Bgee expression fetch failed for %s taxid %d: %s", bgee_gene_id, species_taxid, exc)
        return []


def _is_relevant_tissue(tissue_name: str, trait_tissues: list[str]) -> bool:
    """Check if a Bgee tissue name is relevant to the analysis trait."""
    tissue_lower = tissue_name.lower()
    return any(t in tissue_lower for t in trait_tissues)


def _query_bgee_for_gene(
    gene: Gene,
    trait_tissues: list[str],
) -> tuple[str, list[dict]]:
    """Query Bgee for all supported species for one gene.

    Returns (gene_id, [ExpressionResult-kwargs, ...]) — does NOT write to DB.
    Each dict is passed directly to ExpressionResult(**kwargs).
    """
    if not gene.gene_symbol:
        return gene.id, []

    bgee_human_id = _bgee_gene_search(gene.gene_symbol, species_taxid=9606)
    if not bgee_human_id:
        log.debug("  Bgee: no result for %s", gene.gene_symbol)
        return gene.id, []

    time.sleep(0.1)  # small per-gene courtesy pause

    expression_rows: list[dict] = []
    present_species: list[str] = []

    for taxid, species_id in BGEE_SUPPORTED_TAXIDS.items():
        if species_id == "human":
            continue

        bgee_species_id = _bgee_gene_search(gene.gene_symbol, species_taxid=taxid)
        if not bgee_species_id:
            continue

        calls = _bgee_expression_calls(bgee_species_id, taxid)
        time.sleep(0.05)

        for call in calls:
            tissue = ""
            anat = call.get("anatEntity") or {}
            if isinstance(anat, dict):
                tissue = anat.get("anatEntityName", "")
            elif isinstance(anat, str):
                tissue = anat

            call_type = (call.get("expressionState") or "").lower()
            if call_type in ("expressed", "present"):
                if not trait_tissues or _is_relevant_tissue(tissue, trait_tissues):
                    present_species.append(species_id)
                    expression_rows.append({
                        "gene_id": gene.id,
                        "geo_accession": f"bgee:{bgee_species_id}",
                        "comparison": species_id,
                        "log2fc": None,
                        "padj": None,
                        "n_samples": None,
                    })
                    break

    if len(present_species) >= 2:
        log.debug(
            "  Bgee: %s expressed in %d species (%s)",
            gene.gene_symbol, len(present_species), ", ".join(present_species[:4]),
        )

    return gene.id, expression_rows


def run_bgee_pipeline(
    gene_ids: Optional[list[str]] = None,
    trait_id: str = "",
) -> int:
    """Query Bgee for cross-species expression data for candidate genes.

    Fetches all genes in parallel (thread pool), then writes results to DB
    in a single session commit.

    Args:
        gene_ids: Optional list of gene IDs to process (default: Tier1+Tier2).
        trait_id: Trait preset id for tissue relevance mapping.

    Returns:
        Number of ExpressionResult rows written.
    """
    trait_tissues = TRAIT_TISSUES.get(trait_id, TRAIT_TISSUES["default"])

    with get_session() as session:
        if gene_ids is None:
            rows = (
                session.query(CandidateScore.gene_id)
                .filter(
                    CandidateScore.trait_id == (trait_id or ""),
                    CandidateScore.tier.in_(["Tier1", "Tier2", "Validated"]),
                )
                .all()
            )
            gene_ids = [r[0] for r in rows]

        if not gene_ids:
            log.info("  Bgee: no Tier1/Tier2 genes found — skipping.")
            return 0

        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()

    log.info("Querying Bgee for %d genes (trait_id=%r, %d workers)...",
             len(genes), trait_id, _BGEE_WORKERS)

    all_rows: list[dict] = []
    done = 0
    total = len(genes)

    with ThreadPoolExecutor(max_workers=_BGEE_WORKERS) as pool:
        futures = {pool.submit(_query_bgee_for_gene, gene, trait_tissues): gene.id
                   for gene in genes}
        for future in as_completed(futures):
            _gid, rows = future.result()
            all_rows.extend(rows)
            done += 1
            if done % 50 == 0 or done == total:
                log.info("  Bgee: %d / %d genes processed (%.0f%%).",
                         done, total, 100 * done / total)

    # Single session write for all results
    written = 0
    if all_rows:
        with get_session() as session:
            for kwargs in all_rows:
                session.add(ExpressionResult(**kwargs))
            session.commit()
            written = len(all_rows)

    log.info("Bgee supplement complete: %d expression rows written.", written)
    return written
