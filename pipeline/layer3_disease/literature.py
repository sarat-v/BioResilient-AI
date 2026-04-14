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
ENTREZ_EFETCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
REQUEST_TIMEOUT = 20
MAX_ABSTRACTS_PER_GENE = 20   # fetch up to 20 abstracts for semantic scoring

_SPECIES_SUFFIXES = ("_HUMAN", "_MOUSE", "_RAT", "_DOG", "_CAT", "_WHALE", "_BAT", "_SHARK")


def _hgnc_symbol(full_symbol: str) -> str:
    """Strip UniProt species suffix to get clean HGNC gene symbol."""
    for suffix in _SPECIES_SUFFIXES:
        if full_symbol.upper().endswith(suffix):
            return full_symbol[: -len(suffix)]
    return full_symbol

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
    # "viral_resistance" kept as alias; canonical key is "viral_immunity"
    # (matches functional_evidence_config.json).
    "viral_resistance": [
        "antiviral", "innate immunity", "interferon", "viral resistance",
        "pathogen resistance",
    ],
    "viral_immunity": [
        "antiviral", "innate immunity", "interferon", "viral immunity",
        "viral resistance", "pathogen resistance", "immune evasion",
        "toll-like receptor", "type I interferon", "ISG",
    ],
    "hypoxia_tolerance": [
        "hypoxia", "hypoxia tolerance", "oxygen deprivation", "HIF",
        "angiogenesis", "low oxygen", "ischemia", "altitude adaptation",
        "hypoxia-inducible factor",
    ],
    "dna_repair": [
        "DNA repair", "double-strand break", "base excision repair",
        "nucleotide excision repair", "mismatch repair", "DNA damage response",
        "genomic stability", "mutagenesis", "synthetic lethality",
    ],
    "regeneration": [
        "regeneration", "tissue repair", "wound healing", "stem cell",
        "dedifferentiation",
    ],
    "default": [
        "resilience", "adaptation", "stress response", "protective variant",
        "evolutionary conservation",
    ],
}


def _pubmed_headers(email: Optional[str] = None) -> dict:
    return {"User-Agent": f"BioResilienceAI/1.0 ({email or 'noreply@example.com'})"}


def _query_pubmed(gene_symbol: str, keywords: list[str]) -> tuple[int, list[str]]:
    """Search PubMed and return (count, list_of_pmids).

    Fetches up to MAX_ABSTRACTS_PER_GENE PMIDs for semantic scoring.
    """
    api_key = get_ncbi_api_key()
    email   = get_ncbi_email()

    kw_query  = " OR ".join(f'"{kw}"' for kw in keywords[:5])
    full_query = f'"{gene_symbol}"[Gene Name] AND ({kw_query})'

    params: dict = {
        "db": "pubmed",
        "term": full_query,
        "retmax": MAX_ABSTRACTS_PER_GENE,
        "retmode": "json",
        "usehistory": "n",
    }
    if api_key:
        params["api_key"] = api_key

    try:
        r = requests.get(
            ENTREZ_ESEARCH, params=params,
            headers=_pubmed_headers(email), timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            return 0, []
        data   = r.json().get("esearchresult", {})
        count  = int(data.get("count", 0))
        pmids  = data.get("idlist", [])
        return count, pmids
    except Exception as exc:
        log.debug("PubMed search failed for %s: %s", gene_symbol, exc)
        return 0, []


def _fetch_abstracts(pmids: list[str]) -> list[str]:
    """Fetch abstract text for a list of PMIDs via efetch."""
    if not pmids:
        return []
    api_key = get_ncbi_api_key()
    email   = get_ncbi_email()
    params: dict = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "rettype": "abstract",
        "retmode": "text",
    }
    if api_key:
        params["api_key"] = api_key
    try:
        r = requests.get(
            ENTREZ_EFETCH, params=params,
            headers=_pubmed_headers(email), timeout=REQUEST_TIMEOUT * 2,
        )
        if r.status_code != 200:
            return []
        # Split by paper delimiter — each abstract block starts with a number line
        text  = r.text
        blocks = [b.strip() for b in text.split("\n\n\n") if len(b.strip()) > 50]
        return blocks
    except Exception as exc:
        log.debug("PubMed efetch failed: %s", exc)
        return []


def _semantic_relevance_score(
    gene_symbol: str,
    abstracts: list[str],
    keywords: list[str],
) -> float:
    """Compute a 0-1 semantic relevance score from abstract text.

    Algorithm (deterministic, no LLM):
      1. Count keyword co-occurrences with gene_symbol across all abstracts.
      2. Score = matched_docs / total_docs (proportion of abstracts mentioning
         both the gene and at least one phenotype keyword).
      3. Clamp to [0, 1].

    This outperforms a raw paper count because a gene with 5 highly relevant
    papers scores higher than a gene with 200 tangentially related papers.
    """
    if not abstracts:
        return 0.0

    symbol_lower  = gene_symbol.lower()
    kws_lower     = [kw.lower() for kw in keywords]

    # Phenotype-neutral expansion: only add terms that apply to any resilience trait.
    # Do NOT include cancer-centric terms (tumor suppressor, oncogene) here — they
    # would bias abstract scoring for non-cancer phenotypes. The caller passes
    # phenotype-specific keywords via TRAIT_KEYWORDS; that is the right place.
    mechanism_terms = [
        "protective", "resistance", "resilience", "adaptation",
        "evolutionary conservation", "positive selection",
    ]
    all_terms = kws_lower + mechanism_terms

    matched = 0
    for abstract in abstracts:
        ab_lower = abstract.lower()
        if symbol_lower in ab_lower:
            if any(term in ab_lower for term in all_terms):
                matched += 1

    return round(matched / len(abstracts), 4)


def _query_pubmed_count(gene_symbol: str, keywords: list[str]) -> int:
    """Return PubMed paper count (backward-compat wrapper)."""
    count, _ = _query_pubmed(gene_symbol, keywords)
    return count


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

    import math
    for gene in genes:
        if not gene.gene_symbol:
            continue

        clean_symbol = _hgnc_symbol(gene.gene_symbol)

        # Phase 1: get count + top PMIDs
        count, pmids = _query_pubmed(clean_symbol, keywords)
        time.sleep(0.12)  # NCBI rate limit: ~10 req/sec with API key

        # Phase 2: fetch abstracts and compute semantic relevance score
        lit_score: float
        if pmids:
            abstracts = _fetch_abstracts(pmids[:MAX_ABSTRACTS_PER_GENE])
            time.sleep(0.12)
            if abstracts:
                lit_score = _semantic_relevance_score(clean_symbol, abstracts, keywords)
            else:
                # Fallback to log-count if efetch fails
                lit_score = round(min(math.log10(1 + count) / 3.0, 1.0), 4)
        else:
            lit_score = 0.0

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
        if clean_symbol.upper() in KNOWN_RESILIENCE_GENES:
            found_known.append(clean_symbol)

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
