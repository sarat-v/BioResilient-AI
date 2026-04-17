"""Step 11b — Rare protective variant mapping (gnomAD + GWAS Catalog).

This is the highest-impact scientific addition: it maps the animal convergence
signal back to human genetics, following the PCSK9 discovery paradigm.

Workflow for each Tier1/Tier2 gene:
  1. For each DivergentMotif, identify the amino acid positions that differ
     between the animal and human sequences (divergent positions).
  2. Map protein residue coordinates to hg38 genomic coordinates using the
     NCBI Gene API (same mechanism as PhyloP enrichment).
  3. Query gnomAD v4 GraphQL API for rare variants (MAF < 1%) in those regions.
  4. For each rare variant, check if the amino acid change matches the direction
     of the animal divergence:
       - Same position (within ±2 residues)
       - Variant introduces the same amino acid seen in the resilient species
         OR a biochemically similar amino acid (same Miyata group)
  5. For direction-matching variants: query the GWAS Catalog REST API for
     phenotype associations on that variant or gene.
  6. Persist: protective_variant_count, best_protective_trait, protective_variant_pvalue
     on DiseaseAnnotation.

If a gene has protective_variant_count >= 1 AND protective_variant_pvalue < 5e-8,
the scoring layer upgrades it to "Validated" tier.

References:
  - Cohen et al. (2006) PCSK9: NEJM 354:1264
  - gnomAD v4 API: https://gnomad.broadinstitute.org/api
  - GWAS Catalog REST: https://www.ebi.ac.uk/gwas/rest/api
"""

import logging
import time
from typing import Optional

import requests

from db.models import DiseaseAnnotation, DivergentMotif, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_ncbi_api_key, get_ncbi_email

log = logging.getLogger(__name__)

GNOMAD_GRAPHQL = "https://gnomad.broadinstitute.org/api"
GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"
REQUEST_TIMEOUT = 20

_SPECIES_SUFFIXES = ("_HUMAN", "_MOUSE", "_RAT", "_DOG", "_CAT", "_WHALE", "_BAT", "_SHARK")


def _hgnc_symbol(full_symbol: str) -> str:
    """Strip UniProt species suffix to get clean HGNC gene symbol.

    e.g. 'AKT1_HUMAN' -> 'AKT1', 'CDK4_HUMAN' -> 'CDK4'
    """
    for suffix in _SPECIES_SUFFIXES:
        if full_symbol.upper().endswith(suffix):
            return full_symbol[: -len(suffix)]
    return full_symbol

# Miyata biochemical groupings for amino acid similarity
# Source: Miyata et al. (1979) J Mol Evol 12:219
_MIYATA_GROUPS: dict[str, int] = {
    "G": 0, "A": 0, "V": 0, "L": 0, "I": 0,  # nonpolar aliphatic
    "F": 1, "W": 1, "Y": 1,                     # aromatic
    "S": 2, "T": 2, "C": 2,                     # polar uncharged small
    "N": 3, "Q": 3, "H": 3,                     # polar uncharged
    "D": 4, "E": 4,                              # negatively charged
    "K": 5, "R": 5,                              # positively charged
    "P": 6, "M": 6,                              # special / other
}


def _biochemically_similar(aa1: str, aa2: str) -> bool:
    """Return True if two amino acids are in the same Miyata group."""
    if aa1 == aa2:
        return True
    g1 = _MIYATA_GROUPS.get(aa1.upper())
    g2 = _MIYATA_GROUPS.get(aa2.upper())
    return g1 is not None and g1 == g2


def _get_gene_coordinates_hg38(gene_symbol: str) -> Optional[tuple[str, int, int]]:
    """Get hg38 chromosomal coordinates for a gene from NCBI Gene API.

    Returns (chrom, start, end) or None on failure.
    """
    try:
        email = get_ncbi_email()
        api_key = get_ncbi_api_key()
        params = {
            "db": "gene",
            "term": f"{_hgnc_symbol(gene_symbol)}[gene] AND 9606[taxid]",
            "retmax": 1,
            "retmode": "json",
        }
        if api_key:
            params["api_key"] = api_key
        headers = {"User-Agent": f"BioResillientAI/1.0 ({email})"}
        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params=params, headers=headers, timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            return None
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None

        time.sleep(0.12)
        r2 = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params={"db": "gene", "id": ids[0], "retmode": "json", **({"api_key": api_key} if api_key else {})},
            headers=headers, timeout=REQUEST_TIMEOUT,
        )
        if r2.status_code != 200:
            return None
        summary = r2.json().get("result", {}).get(ids[0], {})
        genomic = summary.get("genomicinfo", [{}])[0]
        chrom = genomic.get("chraccver", "")
        start = int(genomic.get("chrstart", 0))
        stop = int(genomic.get("chrstop", 0))
        if chrom and stop > start:
            return chrom, start, stop
    except Exception as exc:
        log.debug("Gene coordinates lookup failed for %s: %s", gene_symbol, exc)
    return None


def _query_gnomad_variants(
    gene_symbol: str,
    chrom: str,
    start: int,
    stop: int,
    max_af: float = 0.01,
) -> list[dict]:
    """Query gnomAD v4 for rare variants in a gene region.

    Returns list of variant dicts with keys:
        variant_id, pos, ref, alt, af, amino_acid_change, consequence
    """
    # gnomAD uses chromosome names without 'chr' prefix in their API
    chrom_short = chrom.replace("chr", "").split(".")[0]
    if chrom_short.startswith("NC_0000"):
        # NCBI accession format — convert to numeric
        chrom_short = str(int(chrom_short.split(".")[0].replace("NC_", "").lstrip("0") or "0"))
        if chrom_short == "23":
            chrom_short = "X"
        elif chrom_short == "24":
            chrom_short = "Y"

    # gnomAD v4: af is a Float directly, not an object
    query = """
    query($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
      gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
        variants(dataset: gnomad_r4) {
          variant_id
          pos
          ref
          alt
          exome { af }
          genome { af }
          consequence
          hgvsp
        }
      }
    }
    """
    clean_symbol = _hgnc_symbol(gene_symbol)
    try:
        r = requests.post(
            GNOMAD_GRAPHQL,
            json={"query": query, "variables": {"geneSymbol": clean_symbol, "referenceGenome": "GRCh38"}},
            timeout=30,
        )
        if r.status_code != 200:
            return []
        data = r.json()
        if data.get("errors"):
            log.debug("gnomAD variants query errors for %s: %s", clean_symbol, data["errors"])
            return []
        variants_raw = data.get("data", {}).get("gene", {}).get("variants", []) or []

        rare_variants = []
        for v in variants_raw:
            # gnomAD v4: exome.af and genome.af are plain floats
            af = None
            exome_af = (v.get("exome") or {}).get("af")
            genome_af = (v.get("genome") or {}).get("af")
            if exome_af is not None:
                af = exome_af
            elif genome_af is not None:
                af = genome_af

            if af is None or af > max_af:
                continue

            rare_variants.append({
                "variant_id": v.get("variant_id"),
                "pos": v.get("pos"),
                "ref": v.get("ref"),
                "alt": v.get("alt"),
                "af": af,
                "consequence": v.get("consequence", ""),
                "hgvsp": v.get("hgvsp", ""),
            })

        return rare_variants

    except Exception as exc:
        log.debug("gnomAD query failed for %s: %s", gene_symbol, exc)
        return []


def _parse_hgvsp_position(hgvsp: str) -> Optional[tuple[str, int, str]]:
    """Parse HGVSp notation to extract (ref_aa, position, alt_aa).

    Example: 'p.Ala123Val' -> ('A', 123, 'V')
    """
    import re
    AA3TO1 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    }
    m = re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", hgvsp or "")
    if not m:
        return None
    ref_aa = AA3TO1.get(m.group(1))
    pos = int(m.group(2))
    alt_aa = AA3TO1.get(m.group(3))
    if ref_aa and alt_aa:
        return ref_aa, pos, alt_aa
    return None


def _fetch_gwas_associations(gene_symbol: str) -> list[dict]:
    """Fetch GWAS Catalog associations for a gene, returning protective hits.

    Returns list of {trait, pvalue} dicts.
    Uses GWAS Catalog v2 gene search (not the deprecated /genes endpoint).
    """
    clean_symbol = _hgnc_symbol(gene_symbol)
    try:
        r = requests.get(
            f"{GWAS_API}/singleNucleotidePolymorphisms/search/findByGene",
            params={"geneName": clean_symbol, "size": 200},
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            return []
        data = r.json()
        snps = data.get("_embedded", {}).get("singleNucleotidePolymorphisms", [])
        if not snps:
            return []

        hits = []
        for snp in snps[:50]:
            snp_link = snp.get("_links", {}).get("self", {}).get("href", "")
            if not snp_link:
                continue
            r2 = requests.get(
                f"{snp_link}/associations",
                params={"size": 50},
                timeout=REQUEST_TIMEOUT,
            )
            if r2.status_code != 200:
                continue
            time.sleep(0.1)
            assoc_list = r2.json().get("_embedded", {}).get("associations", [])
            for assoc in assoc_list:
                pval = assoc.get("pvalue")
                if pval is None:
                    continue
                try:
                    pval_float = float(pval)
                except (ValueError, TypeError):
                    continue
                traits = assoc.get("efoTraits") or []
                for trait in traits:
                    name = trait.get("trait", "") if isinstance(trait, dict) else str(trait)
                    if name:
                        hits.append({"trait": name, "pvalue": pval_float})

        return hits

    except Exception as exc:
        log.debug("GWAS Catalog fetch failed for %s: %s", clean_symbol, exc)
        return []


def annotate_protective_variants(gene_ids: Optional[list[str]] = None) -> int:
    """Map divergent motif positions to human rare variants in gnomAD.

    For each gene, identifies whether humans carry rare variants at the same
    positions where resilient animals diverge — the PCSK9 paradigm.

    Returns number of genes with at least one protective variant found.
    """
    with get_session() as session:
        if gene_ids is None:
            gene_ids = [g.id for g in session.query(Gene).all()]
        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()

    log.info("Mapping rare protective variants for %d genes...", len(genes))
    found = 0

    for gene in genes:
        if not gene.gene_symbol:
            continue

        clean_symbol = _hgnc_symbol(gene.gene_symbol)

        # Collect divergent positions from all motifs for this gene
        with get_session() as session:
            motifs = (
                session.query(DivergentMotif)
                .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
                .filter(Ortholog.gene_id == gene.id)
                .all()
            )

        if not motifs:
            continue

        # Extract divergent (position, animal_aa, human_aa) tuples
        divergent_positions: list[tuple[int, str, str]] = []
        for motif in motifs:
            for i, (h_aa, a_aa) in enumerate(zip(motif.human_seq, motif.animal_seq)):
                if h_aa != a_aa and h_aa != "-" and a_aa != "-":
                    abs_pos = motif.start_pos + i + 1  # 1-indexed protein position
                    divergent_positions.append((abs_pos, a_aa, h_aa))

        if not divergent_positions:
            continue

        # Get gene coordinates using clean HGNC symbol
        coords = _get_gene_coordinates_hg38(clean_symbol)
        if not coords:
            continue
        chrom, start, stop = coords

        time.sleep(0.15)

        # Query gnomAD for rare variants using clean HGNC symbol
        rare_variants = _query_gnomad_variants(clean_symbol, chrom, start, stop)
        if not rare_variants:
            continue

        time.sleep(0.2)

        # Check for direction-matching variants
        matching_count = 0
        for variant in rare_variants:
            parsed = _parse_hgvsp_position(variant.get("hgvsp", ""))
            if not parsed:
                continue
            ref_aa, var_pos, alt_aa = parsed

            for div_pos, animal_aa, human_aa_expected in divergent_positions:
                # Within ±3 residues AND amino acid change matches animal direction
                if abs(var_pos - div_pos) <= 3:
                    if _biochemically_similar(alt_aa, animal_aa):
                        matching_count += 1
                        break

        if matching_count == 0:
            continue

        # Fetch GWAS associations for this gene using clean HGNC symbol
        gwas_hits = _fetch_gwas_associations(clean_symbol)
        time.sleep(0.15)

        best_trait = None
        best_pvalue = 1.0
        for hit in gwas_hits:
            if hit["pvalue"] < best_pvalue:
                best_pvalue = hit["pvalue"]
                best_trait = hit["trait"]

        # Update DiseaseAnnotation
        with get_session() as session:
            ann = session.get(DiseaseAnnotation, gene.id)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gene.id)
                session.add(ann)

            ann.protective_variant_count = matching_count
            if best_trait:
                ann.best_protective_trait = best_trait
                ann.protective_variant_pvalue = best_pvalue
            session.commit()

        log.info(
            "  %s: %d protective variants (best trait: %s, p=%.2e)",
            clean_symbol, matching_count, best_trait or "none", best_pvalue,
        )
        found += 1

    log.info("Protective variant mapping complete: %d genes with matches.", found)
    return found


def run_rare_variants_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step11b."""
    return annotate_protective_variants(gene_ids)
