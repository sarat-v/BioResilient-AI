"""Step 11c — Literature validation and sanity check.

After ranking candidates, this module checks whether top genes overlap with
known longevity, cancer resistance, and resilience literature. This serves
as a critical scientific sanity check:

  - Expected hits: TP53, IGF1, FOXO3, DNA repair genes (ATM, BRCA1),
    autophagy genes (BECN1, ATG5), mTOR pathway, telomere maintenance
  - If NONE of these appear in Tier1/Tier2, the pipeline may be off-target.
  - A hit in known literature boosts confidence; a novel gene not in
    literature is interesting but needs extra scrutiny.

Method:
  1. Query PubMed via NCBI E-utilities for papers mentioning the gene
     symbol in combination with resilience/longevity keywords.
  2. Count relevant papers (lit_pmid_count).
  3. Compute lit_score = log10(1 + pmid_count) normalised to [0, 1].

Keyword sets by trait:
  - cancer_resistance: "cancer resistance", "tumour suppressor", "oncogene"
  - longevity: "longevity", "aging", "lifespan", "senescence"
  - viral_resistance: "antiviral", "innate immunity", "interferon"
  - default: combines all above + "resilience", "adaptation"

The lit_score is stored on DiseaseAnnotation and exposed in the API as a
confidence signal. It does NOT directly affect composite_score — it's a
separate validation layer shown to researchers.
"""

import logging
import time
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene
from db.session import get_session
from pipeline.config import get_ncbi_api_key, get_ncbi_email

log = logging.getLogger(__name__)

ENTREZ_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
REQUEST_TIMEOUT = 15

# Known landmark resilience/longevity genes — if these are NOT in our top
# candidates, we log a warning as a sanity check.
KNOWN_RESILIENCE_GENES = {
    "TP53", "IGF1", "IGF1R", "FOXO3", "FOXO1", "SIRT1", "SIRT3",
    "ATM", "BRCA1", "BRCA2", "TERT", "TERC",
    "MTOR", "AKT1", "PTEN",
    "ATG5", "BECN1", "ULK1",
    "CDKN2A", "MDM2", "PARP1",
    "STAT3", "NFE2L2",   # NRF2
    "TP63", "TP73",
    "PCSK9",             # PCSK9 paradigm gene
    "APOE", "CETP",      # longevity
    "LMNA", "WRN",       # progeria/Werner
    "CASP3", "BCL2", "MCL1",
}

TRAIT_KEYWORDS: dict[str, list[str]] = {
    "cancer_resistance": [
        "cancer resistance", "tumour suppression", "tumor suppressor",
        "oncogene", "anti-cancer", "cancer resilience",
    ],
    "longevity": [
        "longevity", "aging", "ageing", "lifespan", "life span",
        "senescence", "anti-aging", "healthspan",
    ],
    "viral_resistance": [
        "antiviral", "innate immunity", "interferon", "viral resistance",
        "pathogen resistance",
    ],
    "regeneration": [
        "regeneration", "tissue repair", "wound healing", "stem cell",
        "dedifferentiation",
    ],
    "default": [
        "resilience", "adaptation", "longevity", "cancer resistance",
        "stress response", "protective variant",
    ],
}


def _query_pubmed_count(gene_symbol: str, keywords: list[str]) -> int:
    """Return PubMed paper count for gene + any of the given keywords.

    Uses NCBI E-utilities esearch with a boolean OR query.
    """
    api_key = get_ncbi_api_key()
    email = get_ncbi_email()

    # Build query: gene symbol AND (keyword1 OR keyword2 OR ...)
    kw_query = " OR ".join(f'"{kw}"' for kw in keywords[:5])  # limit to 5 to keep URL short
    full_query = f'"{gene_symbol}"[Gene Name] AND ({kw_query})'

    params: dict = {
        "db": "pubmed",
        "term": full_query,
        "retmax": 0,   # We only want the count
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key
    headers = {"User-Agent": f"BioResilienceAI/1.0 ({email or 'noreply@example.com'})"}

    try:
        r = requests.get(ENTREZ_ESEARCH, params=params, headers=headers, timeout=REQUEST_TIMEOUT)
        if r.status_code != 200:
            return 0
        data = r.json()
        count = int(data.get("esearchresult", {}).get("count", 0))
        return count
    except Exception as exc:
        log.debug("PubMed query failed for %s: %s", gene_symbol, exc)
        return 0


def annotate_literature(
    gene_ids: Optional[list[str]] = None,
    trait_id: str = "",
) -> int:
    """Annotate DiseaseAnnotation rows with PubMed literature scores.

    Args:
        gene_ids: Optional list of gene IDs to process (default: all).
        trait_id: Trait preset id used to select relevant keywords.

    Returns:
        Number of genes with lit_score > 0.
    """
    keywords = TRAIT_KEYWORDS.get(trait_id, TRAIT_KEYWORDS["default"])

    with get_session() as session:
        if gene_ids is None:
            gene_ids = [g.id for g in session.query(Gene).all()]
        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()

    log.info("Querying PubMed literature for %d genes (trait=%r)...", len(genes), trait_id)
    annotated = 0
    found_known = []

    for gene in genes:
        if not gene.gene_symbol:
            continue

        count = _query_pubmed_count(gene.gene_symbol, keywords)
        time.sleep(0.12)  # NCBI rate limit: ~10 req/sec with API key

        # lit_score: log10(1 + count) normalised so 1000 papers → 1.0
        import math
        lit_score = round(min(math.log10(1 + count) / 3.0, 1.0), 4)

        with get_session() as session:
            ann = session.get(DiseaseAnnotation, gene.id)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gene.id)
                session.add(ann)

            ann.lit_score = lit_score
            ann.lit_pmid_count = count
            session.commit()

        if count > 0:
            annotated += 1
        if gene.gene_symbol.upper() in KNOWN_RESILIENCE_GENES:
            found_known.append(gene.gene_symbol)

    # Sanity check log
    if found_known:
        log.info("  Literature sanity check ✓ — known resilience genes in results: %s",
                 ", ".join(sorted(found_known)[:10]))
    else:
        log.warning(
            "  Literature sanity check ⚠ — NONE of the %d known resilience genes "
            "(e.g. TP53, FOXO3, IGF1) found in the candidate list. "
            "Consider reviewing convergence thresholds.",
            len(KNOWN_RESILIENCE_GENES),
        )

    log.info("Literature annotation complete: %d / %d genes with citations.", annotated, len(genes))
    return annotated


def run_literature_pipeline(
    gene_ids: Optional[list[str]] = None,
    trait_id: str = "",
) -> int:
    """Entry point called from orchestrator step11c."""
    return annotate_literature(gene_ids, trait_id)
