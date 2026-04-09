"""PheWAS / GWAS Catalog phenome-wide associations — safety signal.

LD-filter applied: only associations where the lead SNP location falls within
the gene body (chrStart–chrEnd) ± 1 kb are kept. This eliminates a common
source of false-positive safety signals driven by linkage disequilibrium with
neighbouring genes rather than the target gene itself.

Gene genomic coordinates are fetched from the GWAS Catalog gene endpoint
(chromosomalLocation field). If coordinates cannot be resolved the filter is
skipped conservatively (all hits retained with a warning flag).
"""

import logging
from typing import Any, Optional

import requests

from db.models import Gene, SafetyFlag
from db.session import get_session

log = logging.getLogger(__name__)

GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"

# Flanking window around gene body for SNP inclusion (bp)
_GENE_FLANK_BP = 1_000


def _fetch_gene_coordinates(gene_symbol: str) -> Optional[tuple[str, int, int]]:
    """Return (chromosome, start_bp, end_bp) for a gene from GWAS Catalog.

    Returns None if the gene is not found or has no location data.
    """
    try:
        r = requests.get(
            f"{GWAS_API}/genes",
            params={"geneName": gene_symbol, "size": 1},
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        genes = data.get("_embedded", {}).get("genes", [])
        if not genes:
            return None
        loc = genes[0].get("chromosomalLocation", {})
        chrom = loc.get("chromosomeName")
        start = loc.get("chromosomeStart") or loc.get("geneStart")
        end = loc.get("chromosomeEnd") or loc.get("geneStop")
        if chrom and start is not None and end is not None:
            return str(chrom), int(start), int(end)
    except Exception as exc:
        log.debug("GWAS gene coordinates %s: %s", gene_symbol, exc)
    return None


def _snp_location(association: dict) -> Optional[tuple[str, int]]:
    """Extract (chromosome, position_bp) from a GWAS Catalog association dict.

    Tries multiple paths in the association JSON structure.
    """
    # Path 1: loci → strongestRiskAlleles → (no position here usually)
    # Path 2: loci → genomicContexts → gene (but we need SNP position)
    # Path 3: directly on association SNPs
    snps = association.get("snps") or []
    for snp in snps:
        locs = snp.get("locations") or []
        for loc in locs:
            chrom = loc.get("chromosomeName")
            pos = loc.get("chromosomePosition")
            if chrom and pos is not None:
                return str(chrom), int(pos)
    # Fallback: loci list
    loci = association.get("loci") or []
    for locus in loci:
        for allele in locus.get("strongestRiskAlleles") or []:
            # Some endpoints embed location directly on the risk allele
            snp_detail = allele.get("snp") or {}
            for loc in snp_detail.get("locations") or []:
                chrom = loc.get("chromosomeName")
                pos = loc.get("chromosomePosition")
                if chrom and pos is not None:
                    return str(chrom), int(pos)
    return None


def fetch_phewas_hits(gene_symbol: str) -> Optional[dict[str, Any]]:
    """Fetch phenotype–pvalue map from GWAS Catalog for the gene (phenome-wide style).

    Returns {trait: pvalue} for associations where the lead SNP falls within
    the gene body ± _GENE_FLANK_BP. If coordinates cannot be resolved, all
    associations are returned (conservative fallback).

    The key improvement over the original implementation: LD-driven spurious
    safety signals from neighbouring genes are removed by the coordinate filter.
    """
    try:
        # Step 1: resolve gene
        r = requests.get(
            f"{GWAS_API}/genes",
            params={"geneName": gene_symbol, "size": 1},
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        embeds = data.get("_embedded", {}).get("genes", [])
        if not embeds:
            return None
        self_link = embeds[0].get("_links", {}).get("self", {}).get("href", "")
        if not self_link.startswith("/"):
            return None

        # Step 2: fetch associations
        r2 = requests.get(
            f"https://www.ebi.ac.uk{self_link}/associations",
            params={"size": 500},
            timeout=15,
        )
        if r2.status_code != 200:
            return None
        assoc = r2.json()
        assoc_list = assoc.get("_embedded", {}).get("associations", [])
        if not assoc_list:
            return None

        # Step 3: get gene genomic window for LD filter
        gene_coords = _fetch_gene_coordinates(gene_symbol)
        if gene_coords:
            g_chrom, g_start, g_end = gene_coords
            window_start = g_start - _GENE_FLANK_BP
            window_end = g_end + _GENE_FLANK_BP
            ld_filter_active = True
        else:
            # Cannot resolve coordinates — keep all hits but log a warning
            log.debug("PheWAS %s: could not resolve gene coordinates; LD filter skipped.", gene_symbol)
            g_chrom = None
            window_start = window_end = 0
            ld_filter_active = False

        hits: dict[str, float] = {}
        n_filtered = 0
        for a in assoc_list:
            # Apply LD coordinate filter
            if ld_filter_active:
                snp_loc = _snp_location(a)
                if snp_loc is not None:
                    snp_chrom, snp_pos = snp_loc
                    # Normalise chromosome names (remove "chr" prefix if present)
                    snp_chrom_norm = snp_chrom.replace("chr", "")
                    g_chrom_norm = g_chrom.replace("chr", "")  # type: ignore[union-attr]
                    if snp_chrom_norm != g_chrom_norm or not (window_start <= snp_pos <= window_end):
                        n_filtered += 1
                        continue
                # If SNP location is not resolvable, keep the association (conservative)

            efo = a.get("efoTraits") or []
            pval = a.get("pvalue")
            for trait in efo:
                name = trait.get("trait", "") if isinstance(trait, dict) else str(trait)
                if name and pval is not None:
                    try:
                        p = float(pval)
                        if name not in hits or p < hits[name]:
                            hits[name] = p
                    except (TypeError, ValueError):
                        pass

        if n_filtered > 0:
            log.debug("PheWAS %s: removed %d LD-proximal associations outside gene window.", gene_symbol, n_filtered)

        return hits if hits else None
    except Exception as exc:
        log.debug("PheWAS %s: %s", gene_symbol, exc)
        return None


def annotate_genes_phewas(gene_ids: list[str]) -> int:
    """Populate SafetyFlag.phewas_hits for the given genes (with LD filter)."""
    updated = 0
    with get_session() as session:
        for gid in gene_ids:
            gene = session.get(Gene, gid)
            if not gene:
                continue
            hits = fetch_phewas_hits(gene.gene_symbol)
            if hits is None:
                continue
            sf = session.get(SafetyFlag, gid)
            if sf is None:
                sf = SafetyFlag(gene_id=gid)
                session.add(sf)
            sf.phewas_hits = hits
            updated += 1
    log.info("PheWAS: updated %d genes (LD-filtered).", updated)
    return updated
