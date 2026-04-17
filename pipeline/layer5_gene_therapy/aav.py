"""AAV compatibility — CDS size estimate and packaging limit, tissue tropism.

AAV packaging limit is ~4.7 kb for the payload (promoter + CDS + polyA).
We estimate CDS size as: UniProt protein_length × 3 bp/codon + 300 bp UTR buffer.
This is a better proxy than genomic span (which includes introns and can be >100 kb
even for small, AAV-compatible genes like AKT1 ~26 kb genomic / ~1.4 kb CDS).

Genomic span from NCBI is also stored for reference but NOT used for AAV compatibility.

Tissue tropism is inferred from GTEx expression data (top-3 expressed tissues mapped
to AAV serotypes). Stored as NULL when GTEx data is unavailable rather than hardcoding
a generic serotype that would be biologically incorrect for most genes.
"""

import logging
import time
from typing import Optional

import requests

from db.models import Gene, GeneTherapyScore
from db.session import get_session
from pipeline.config import get_ncbi_api_key, get_ncbi_email

log = logging.getLogger(__name__)

AAV_PACKAGING_LIMIT_BP = 4700   # ~4.7 kb payload limit for AAV
UTR_BUFFER_BP = 300             # estimated promoter + UTR overhead
UNIPROT_API = "https://rest.uniprot.org/uniprotkb"

# Maps GTEx tissueSiteDetailId prefix → recommended AAV serotype.
# Based on published tropism studies; references:
#   AAV8  → liver (Wang et al. 2010 Mol Ther)
#   AAV9  → heart, CNS, skeletal muscle (Foust et al. 2009 Nat Biotechnol)
#   AAV5  → lung, retina, CNS
#   AAV1  → skeletal muscle
#   AAV6  → heart, lung
#   AAV2  → CNS, retina
_TISSUE_SEROTYPE = [
    ("Liver",            "AAV8"),
    ("Pancreas",         "AAV8"),
    ("Heart",            "AAV9"),
    ("Brain",            "AAV9"),
    ("Spinal_cord",      "AAV9"),
    ("Nerve_tibial",     "AAV9"),
    ("Muscle",           "AAV1"),
    ("Lung",             "AAV5"),
    ("Retina",           "AAV5"),
    ("Eye",              "AAV2"),
    ("Kidney",           "AAV9"),
    ("Adrenal_Gland",    "AAV9"),
]

NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_ESEARCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"


def _hgnc_symbol(raw: str) -> str:
    """Convert OMA entry name (AKT1_HUMAN) to bare HGNC symbol (AKT1)."""
    if raw and "_" in raw:
        return raw.split("_")[0]
    return raw or ""


def _infer_tissue_tropism(gene_symbol: str) -> Optional[list[str]]:
    """Recommend AAV serotypes based on GTEx median expression across tissues.

    Queries the GTEx Portal API for the gene's tissue expression profile, identifies
    the top-expressed tissues, and maps them to known AAV serotype preferences.
    Returns a list of serotype strings (e.g. ['AAV8', 'AAV9']) in priority order,
    or None if GTEx data is unavailable.

    This is gene-specific: a liver-expressed gene gets AAV8, a cardiac gene gets
    AAV9, etc.  Returning None (rather than a generic fallback) is intentional —
    it signals that tropism was not determined from data.
    """
    try:
        from pipeline.layer6_safety.gtex import fetch_gtex_expression
    except ImportError:
        log.debug("GTEx import unavailable — skipping tissue tropism inference.")
        return None

    expr = fetch_gtex_expression(gene_symbol)
    if not expr:
        log.info("AAV %s: no GTEx expression data — tissue_tropism stored as NULL.", gene_symbol)
        return None

    # Rank tissues by median TPM; take top 5 for matching
    top_tissues = sorted(expr.items(), key=lambda x: x[1], reverse=True)[:5]
    log.info("AAV %s: top GTEx tissues: %s",
             gene_symbol, [(t, round(v, 1)) for t, v in top_tissues])

    serotypes: list[str] = []
    for tissue_id, _tpm in top_tissues:
        for fragment, serotype in _TISSUE_SEROTYPE:
            if fragment.lower() in tissue_id.lower() and serotype not in serotypes:
                serotypes.append(serotype)
                break

    if not serotypes:
        # No tissue match found in top-5; store NULL rather than a generic guess.
        log.info("AAV %s: no tissue-to-serotype match found — tissue_tropism stored as NULL.", gene_symbol)
        return None

    return serotypes[:3]


def _uniprot_accession(human_protein: str) -> Optional[str]:
    """Extract bare UniProt accession from OMA-style 'human|AKT1_HUMAN' or accession directly."""
    if not human_protein:
        return None
    # Try to get accession from human_protein field (could be entry name or accession)
    entry = human_protein.split("|")[-1].strip()  # e.g. AKT1_HUMAN
    try:
        r = requests.get(
            f"{UNIPROT_API}/search",
            params={"query": f"gene:{_hgnc_symbol(entry)} AND organism_id:9606 AND reviewed:true",
                    "fields": "accession,sequence",
                    "format": "json",
                    "size": 1},
            timeout=15,
        )
        if r.status_code == 200:
            results = r.json().get("results", [])
            if results:
                return results[0].get("primaryAccession")
    except Exception as exc:
        log.debug("UniProt lookup %s: %s", entry, exc)
    return None


def _uniprot_protein_length(human_protein: str) -> Optional[int]:
    """Return protein length (aa) from UniProt for a given entry."""
    if not human_protein:
        return None
    entry = human_protein.split("|")[-1].strip()
    symbol = _hgnc_symbol(entry)
    try:
        r = requests.get(
            f"{UNIPROT_API}/search",
            params={"query": f"gene:{symbol} AND organism_id:9606 AND reviewed:true",
                    "fields": "accession,length",
                    "format": "json",
                    "size": 1},
            timeout=15,
        )
        if r.status_code == 200:
            results = r.json().get("results", [])
            if results:
                seq_info = results[0].get("sequence", {})
                return seq_info.get("length")
    except Exception as exc:
        log.debug("UniProt length %s: %s", entry, exc)
    return None


def _ncbi_genomic_span(gene_symbol: str) -> Optional[int]:
    """Get genomic span (bp) from NCBI Gene — stored for reference, NOT used for AAV gating."""
    try:
        r = requests.get(
            NCBI_ESEARCH,
            params={
                "db": "gene",
                "term": f"{gene_symbol}[Gene Name] AND 9606[Taxonomy ID]",
                "retmax": 1,
                "retmode": "json",
                "api_key": get_ncbi_api_key(),
                "email": get_ncbi_email(),
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        id_list = r.json().get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return None
        time.sleep(0.12)
        r2 = requests.get(
            NCBI_ESUMMARY,
            params={
                "db": "gene",
                "id": id_list[0],
                "retmode": "json",
                "api_key": get_ncbi_api_key(),
                "email": get_ncbi_email(),
            },
            timeout=15,
        )
        if r2.status_code != 200:
            return None
        result = r2.json().get("result", {}).get(id_list[0], {})
        genomicinfo = result.get("genomicinfo", [{}])
        if not genomicinfo:
            return None
        start = genomicinfo[0].get("chrstart")
        end = genomicinfo[0].get("chrstop")
        if start is not None and end is not None:
            return abs(int(end) - int(start))  # abs() handles minus-strand genes
    except Exception as exc:
        log.debug("NCBI genomic span %s: %s", gene_symbol, exc)
    return None


def annotate_genes_aav(gene_ids: list[str]) -> int:
    """Populate GeneTherapyScore.gene_size_bp, aav_compatible, tissue_tropism.

    gene_size_bp    = estimated CDS size (protein_aa × 3 + UTR_BUFFER_BP).
    aav_compatible  = gene_size_bp <= AAV_PACKAGING_LIMIT_BP (4700 bp).
    tissue_tropism  = list of recommended AAV serotypes inferred from GTEx
                      expression profile, or NULL if data unavailable.
    """
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            symbol = _hgnc_symbol(gene.gene_symbol)

            # Primary: CDS estimate from UniProt protein length
            protein_aa = _uniprot_protein_length(gene.human_protein)
            if protein_aa:
                cds_estimate = protein_aa * 3 + UTR_BUFFER_BP
            else:
                # Fallback: NCBI genomic span (conservative — large genes look bigger than CDS)
                genomic = _ncbi_genomic_span(symbol)
                cds_estimate = genomic  # may overestimate; aav_compatible may be False

            # Infer tissue tropism from GTEx expression data.
            # This queries the GTEx API — a short sleep before is not needed here
            # because _infer_tissue_tropism already contains the request.
            tropism = _infer_tissue_tropism(symbol)

            gs = session.get(GeneTherapyScore, gid)
            if gs is None:
                gs = GeneTherapyScore(gene_id=gid)
                session.add(gs)

            gs.gene_size_bp = cds_estimate
            gs.aav_compatible = (cds_estimate is not None and cds_estimate <= AAV_PACKAGING_LIMIT_BP)
            gs.tissue_tropism = tropism  # None (NULL) if GTEx data not available
            log.info("AAV %s: protein=%s aa  cds_est=%s bp  aav_compatible=%s  tropism=%s",
                     symbol, protein_aa, cds_estimate, gs.aav_compatible, tropism)
            updated += 1
        session.commit()
    log.info("AAV: updated %d genes.", updated)
    return updated
