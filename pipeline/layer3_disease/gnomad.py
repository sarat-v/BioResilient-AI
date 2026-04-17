"""gnomAD v4 GraphQL API — gene constraint scores (pLI and LOEUF).

pLI  (probability loss-of-function intolerant): classic metric, threshold ≥ 0.9.
LOEUF (loss-of-function observed/expected upper bound fraction): continuous metric.

Threshold used throughout this pipeline: LOEUF < 0.35 (highly LoF-intolerant).
This corresponds to haploinsufficient / dosage-sensitive genes and is the threshold
defined in Karczewski et al. (2020) *Nature* 581:434-443 (gnomAD v2.1.1) for the
"high-confidence LoF intolerant" gene set. It is intentionally more conservative
than the v4 general threshold (0.6) to restrict the loss_of_function motif class
in Step 4d to genes where any LoF variant is robustly depleted from the population.

This threshold is used consistently in:
  - pipeline/layer3_disease/gnomad.py          (this file — constraint fetching)
  - pipeline/layer1_sequence/variant_direction.py (Step 4d motif classification)
  - pipeline/scoring.py disease_score()           (LoF constraint component)

Reference: Karczewski et al. (2020) *Nature* 581:434–443. doi:10.1038/s41586-020-2308-7

Drug-target interpretation:
  - LOEUF < 0.35 (highly intolerant): modulate, do not eliminate; allosteric or
    partial inhibition / PROTAC dosing strategy preferred.
  - LOEUF ≥ 0.35 (tolerant): full LoF therapeutic approach is viable.
  - Presence of healthy LoF carriers in gnomAD: strongest possible safety evidence
    (natural human knockout validated).
  High pLI / low LOEUF does NOT automatically exclude a target — many successful
  drug targets (COX-1, HMGCR) score as LoF intolerant. What matters is the
  direction and magnitude of the desired therapeutic effect.
"""

import logging
import time
from typing import Optional

import requests

from db.models import DiseaseAnnotation, Gene
from db.session import get_session

log = logging.getLogger(__name__)

GNOMAD_GRAPHQL = "https://gnomad.broadinstitute.org/api"

_SPECIES_SUFFIXES = ("_HUMAN", "_MOUSE", "_RAT", "_BOVIN", "_DANRE", "_YEAST", "_CAEEL", "_DROME")


def _hgnc_symbol(gene_symbol: str) -> str:
    """Strip UniProt species suffix to get an HGNC-style gene symbol (AKT1_HUMAN → AKT1)."""
    upper = gene_symbol.upper()
    for suffix in _SPECIES_SUFFIXES:
        if upper.endswith(suffix):
            return gene_symbol[: -len(suffix)]
    return gene_symbol

# gnomAD uses CloudFlare; rapid bursts trigger 429. 0.2 s = ~5 req/s per IP.
# 60 genes × 2 calls each = 120 requests in ~24 s.
_RATE_SLEEP = 0.2

# Constraint thresholds — pipeline-wide constants.
# LOEUF < 0.35: "high-confidence LoF intolerant" as defined in Karczewski et al.
# (2020) Nature 581:434-443. Used identically in variant_direction.py (Step 4d)
# and scoring.py (disease_score). Do NOT change without rerunning Steps 4d, 11–15.
PLI_INTOLERANT_THRESHOLD   = 0.9    # pLI threshold (unchanged from v2/v4)
LOEUF_INTOLERANT_THRESHOLD = 0.35   # LOEUF < 0.35 → highly LoF-intolerant gene


def _symbol_to_ensembl_id(gene_symbol: str) -> Optional[str]:
    """Resolve gene symbol to Ensembl gene ID via gnomAD search.

    gnomAD v4 API requires reference_genome and returns ensembl_id (not gene_id).
    """
    try:
        r = requests.post(
            GNOMAD_GRAPHQL,
            json={
                "query": """
                query GeneSearch($query: String!) {
                  gene_search(query: $query, reference_genome: GRCh38) {
                    ensembl_id
                    symbol
                  }
                }
                """,
                "variables": {"query": gene_symbol},
            },
            timeout=15,
        )
        if r.status_code != 200:
            return None
        data = r.json()
        hits = (data.get("data") or {}).get("gene_search") or []
        for h in hits:
            if h.get("symbol", "").upper() == gene_symbol.upper():
                return h.get("ensembl_id")
        if hits:
            return hits[0].get("ensembl_id")
    except Exception as exc:
        log.debug("gnomAD search for %s: %s", gene_symbol, exc)
    return None


def fetch_gnomad_constraint(ensembl_gene_id: str) -> Optional[dict]:
    """Fetch pLI and LOEUF constraint scores for a gene from gnomAD v4.

    Returns dict with keys 'pli' and 'loeuf', or None on failure.
    gnomAD v4 exposes constraint via genome_version="GRCh38".
    """
    # Try v4 first (gnomAD v4 uses dataset_id "gnomad_r4")
    for dataset in ("gnomad_r4", "gnomad_r2_1"):
        try:
            r = requests.post(
                GNOMAD_GRAPHQL,
                json={
                    "query": """
                    query GeneConstraint($geneId: String!, $referenceGenome: ReferenceGenomeId!) {
                      gene(gene_id: $geneId, reference_genome: $referenceGenome) {
                        gnomad_constraint {
                          pli
                          oe_lof_upper
                        }
                      }
                    }
                    """,
                    "variables": {
                        "geneId": ensembl_gene_id,
                        "referenceGenome": "GRCh38" if dataset == "gnomad_r4" else "GRCh37",
                    },
                },
                timeout=15,
            )
            if r.status_code != 200:
                continue
            data = r.json()
            gene_data = (data.get("data") or {}).get("gene")
            if not gene_data:
                continue
            constraint = gene_data.get("gnomad_constraint") or {}
            pli = constraint.get("pli")
            loeuf = constraint.get("oe_lof_upper")  # gnomAD v4: oe_lof_upper = LOEUF
            if pli is not None or loeuf is not None:
                return {
                    "pli": round(float(pli), 4) if pli is not None else None,
                    "loeuf": round(float(loeuf), 4) if loeuf is not None else None,
                }
        except Exception as exc:
            log.debug("gnomAD constraint (%s) for %s: %s", dataset, ensembl_gene_id, exc)

    # Fallback: legacy pli-only query (v2/v3 schema)
    try:
        r = requests.post(
            GNOMAD_GRAPHQL,
            json={
                "query": """
                query GeneConstraint($geneId: String!) {
                  gene(gene_id: $geneId) { pli }
                }
                """,
                "variables": {"geneId": ensembl_gene_id},
            },
            timeout=15,
        )
        if r.status_code == 200:
            data = r.json()
            gene_data = (data.get("data") or {}).get("gene")
            if gene_data:
                pli = gene_data.get("pli")
                return {"pli": round(float(pli), 4) if pli is not None else None, "loeuf": None}
    except Exception as exc:
        log.debug("gnomAD legacy pli for %s: %s", ensembl_gene_id, exc)

    return None


# Backward-compatible wrapper
def fetch_gnomad_pli(ensembl_gene_id: str) -> Optional[float]:
    result = fetch_gnomad_constraint(ensembl_gene_id)
    return result["pli"] if result else None


def annotate_genes_gnomad(gene_ids: list[str]) -> int:
    """Populate DiseaseAnnotation.gnomad_pli and .gnomad_loeuf for the given genes.

    Idempotent: skips genes where gnomad_pli is already populated.
    """
    with get_session() as session:
        gene_map = {g.id: g for g in session.query(Gene).filter(Gene.id.in_(gene_ids)).all()}
        already_done = {
            r.gene_id
            for r in session.query(DiseaseAnnotation.gene_id)
            .filter(
                DiseaseAnnotation.gene_id.in_(gene_ids),
                DiseaseAnnotation.gnomad_pli.isnot(None),
            )
            .all()
        }

    todo = [gid for gid in gene_ids if gid not in already_done and gid in gene_map]
    if not todo:
        log.info("gnomAD: all %d genes already annotated, skipping.", len(already_done))
        return 0

    # Fetch constraints outside DB session (rate-limited; gnomAD CloudFlare is sensitive)
    fetch_results: dict[str, dict] = {}
    for i, gid in enumerate(todo):
        gene = gene_map[gid]
        ensembl_id = _symbol_to_ensembl_id(_hgnc_symbol(gene.gene_symbol))
        if ensembl_id:
            time.sleep(_RATE_SLEEP)  # between symbol-lookup and constraint call
            constraint = fetch_gnomad_constraint(ensembl_id)
            if constraint:
                fetch_results[gid] = constraint
        if i < len(todo) - 1:
            time.sleep(_RATE_SLEEP)  # between genes

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
        for gid, constraint in fetch_results.items():
            gene = gene_map[gid]
            ann = ann_map.get(gid)
            if ann is None:
                ann = DiseaseAnnotation(gene_id=gid)
                session.add(ann)
            if constraint.get("pli") is not None:
                ann.gnomad_pli = constraint["pli"]
            if constraint.get("loeuf") is not None:
                ann.gnomad_loeuf = constraint["loeuf"]
            updated += 1
            log.debug(
                "gnomAD %s: pLI=%.3f LOEUF=%s (intolerant=%s)",
                gene.gene_symbol,
                constraint.get("pli") or 0,
                constraint.get("loeuf"),
                (constraint.get("loeuf") or 1.0) < LOEUF_INTOLERANT_THRESHOLD,
            )

    log.info("gnomAD: updated %d genes (pLI + LOEUF v4).", updated)
    return updated
