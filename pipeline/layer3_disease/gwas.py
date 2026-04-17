"""GWAS Catalog REST API — strongest association p-value per gene."""

import logging
import time
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene
from db.session import get_session

log = logging.getLogger(__name__)

GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"

_SPECIES_SUFFIXES = ("_HUMAN", "_MOUSE", "_RAT", "_BOVIN", "_DANRE", "_YEAST", "_CAEEL", "_DROME")


def _hgnc_symbol(gene_symbol: str) -> str:
    """Strip UniProt species suffix to get an HGNC-style gene symbol (AKT1_HUMAN → AKT1)."""
    upper = gene_symbol.upper()
    for suffix in _SPECIES_SUFFIXES:
        if upper.endswith(suffix):
            return gene_symbol[: -len(suffix)]
    return gene_symbol

# EBI asks for "reasonable use". 0.25 s between calls = ~4 req/s per IP.
# At 60 genes × 2 calls each that's 120 requests in ~30 s — well within limits.
_RATE_SLEEP = 0.25


def fetch_gwas_pvalue(gene_symbol: str) -> Optional[float]:
    """Fetch the strongest (minimum) p-value from GWAS Catalog for the given gene symbol.

    Uses SNP-based lookup since the /genes endpoint was removed from the GWAS Catalog REST API.
    Steps:
      1. Get up to 20 SNPs mapped to the gene via findByGene.
      2. For each SNP (up to 5), fetch associations and extract p-values from
         pvalueMantissa × 10^pvalueExponent.
      3. Return the minimum p-value across all associations found.
    """
    try:
        r = requests.get(
            f"{GWAS_API}/singleNucleotidePolymorphisms/search/findByGene",
            params={"geneName": gene_symbol, "size": 20},
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        snps = data.get("_embedded", {}).get("singleNucleotidePolymorphisms", [])
        if not snps:
            return None
    except Exception as exc:
        log.debug("GWAS SNP lookup for %s: %s", gene_symbol, exc)
        return None

    pvalues: list[float] = []
    for snp in snps[:5]:
        rs_id = snp.get("rsId")
        if not rs_id:
            continue
        try:
            r2 = requests.get(
                f"{GWAS_API}/singleNucleotidePolymorphisms/{rs_id}/associations",
                params={"size": 20},
                timeout=15,
            )
            if r2.status_code != 200:
                continue
            assoc_data = r2.json()
            for a in assoc_data.get("_embedded", {}).get("associations", []):
                mantissa = a.get("pvalueMantissa")
                exponent = a.get("pvalueExponent")
                if mantissa is not None and exponent is not None:
                    try:
                        m = float(mantissa)
                        # Skip zero-mantissa entries — these are GWAS Catalog data-quality
                        # artifacts (mantissa stored as 0 instead of a true p-value).
                        # min(pvalues) would otherwise pick up 0.0 as the "best" hit.
                        if m <= 0:
                            continue
                        pvalues.append(m * (10 ** int(exponent)))
                    except (TypeError, ValueError):
                        pass
            time.sleep(_RATE_SLEEP)
        except Exception as exc:
            log.debug("GWAS assoc for %s/%s: %s", gene_symbol, rs_id, exc)
            continue

    return min(pvalues) if pvalues else None


def annotate_genes_gwas(gene_ids: list[str]) -> int:
    """Populate DiseaseAnnotation.gwas_pvalue for the given genes.

    Idempotent: skips genes where gwas_pvalue is already populated.
    Bulk write: all DB writes in one session after all fetches complete.
    """
    with get_session() as session:
        gene_map = {g.id: g for g in session.query(Gene).filter(Gene.id.in_(gene_ids)).all()}
        already_done = {
            r.gene_id
            for r in session.query(DiseaseAnnotation.gene_id)
            .filter(
                DiseaseAnnotation.gene_id.in_(gene_ids),
                DiseaseAnnotation.gwas_pvalue.isnot(None),
            )
            .all()
        }

    todo = [gid for gid in gene_ids if gid not in already_done and gid in gene_map]
    if not todo:
        log.info("GWAS Catalog: all %d genes already annotated, skipping.", len(already_done))
        return 0

    # Fetch all p-values outside DB session (rate-limited to avoid EBI 429)
    fetch_results: dict[str, float] = {}
    for i, gid in enumerate(todo):
        pval = fetch_gwas_pvalue(_hgnc_symbol(gene_map[gid].gene_symbol))
        if pval is not None:
            fetch_results[gid] = pval
        if i < len(todo) - 1:
            time.sleep(_RATE_SLEEP)

    if not fetch_results:
        return 0

    # Bulk write in one session
    updated = 0
    with get_session() as session:
        ann_map = {
            r.gene_id: r
            for r in session.query(DiseaseAnnotation)
            .filter(DiseaseAnnotation.gene_id.in_(list(fetch_results)))
            .all()
        }
        for gid, pval in fetch_results.items():
            ann = ann_map.get(gid)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gid)
                session.add(ann)
            ann.gwas_pvalue = pval
            updated += 1

    log.info("GWAS Catalog: updated %d genes.", updated)
    return updated
