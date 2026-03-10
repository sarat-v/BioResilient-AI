"""Variant direction inference: GoF / LoF / neutral classification per divergent motif.

Combines three existing signals that are already computed and stored on DivergentMotif:
  - AlphaMissense consequence_score  (0 = benign, 1 = pathogenic)
  - ESM-1v LLR score                  (negative = destabilising, positive = neutral/gain)
  - LOEUF from gnomAD                 (fetched once per human gene; low = LoF intolerant)

Classification logic:
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ destabilising (ESM-1v LLR < -2)  + high AM score (> 0.6)               │
  │   + LoF-tolerant LOEUF > 0.5     → loss_of_function (protective LoF)   │
  │ destabilising                    + LoF-intolerant LOEUF ≤ 0.5          │
  │                                   → likely_pathogenic (filter out)      │
  │ neutral/stabilising ESM-1v        + high AM  → gain_of_function         │
  │ all other                          → neutral                             │
  └──────────────────────────────────────────────────────────────────────────┘

Result stored as DivergentMotif.motif_direction (string enum: 'gain_of_function',
'loss_of_function', 'likely_pathogenic', 'neutral').
"""

import logging
import time
from typing import Optional

import requests

from db.session import get_session

log = logging.getLogger(__name__)

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_GENE_QUERY = """
{
  gene(gene_symbol: "%s", reference_genome: GRCh38) {
    pLI
    LOEUF: gnomad_constraint { loe_uf }
  }
}
"""


# ---------------------------------------------------------------------------
# gnomAD LOEUF fetch (cached in-process for a single run)
# ---------------------------------------------------------------------------

_loeuf_cache: dict[str, Optional[float]] = {}


def _fetch_loeuf(gene_symbol: str) -> Optional[float]:
    """Fetch LOEUF (LoF o/e upper CI) for a human gene from gnomAD GraphQL API.

    Returns a float in [0, 1+]. Lower = more LoF intolerant.
    Caches results to avoid repeated requests.
    """
    if gene_symbol in _loeuf_cache:
        return _loeuf_cache[gene_symbol]

    try:
        r = requests.post(
            GNOMAD_API,
            json={"query": GNOMAD_GENE_QUERY % gene_symbol},
            timeout=15,
        )
        if r.status_code == 200:
            data = r.json()
            gene_data = data.get("data", {}).get("gene", {})
            # Try gnomad_constraint first, then pLI as fallback
            constraint = gene_data.get("LOEUF", {})
            if isinstance(constraint, dict):
                loeuf = constraint.get("loe_uf")
            else:
                loeuf = None
            # pLI fallback: invert (pLI close to 1 = intolerant → low LOEUF equivalent)
            if loeuf is None:
                pli = gene_data.get("pLI")
                if pli is not None:
                    loeuf = 1.0 - float(pli)   # approximate
            _loeuf_cache[gene_symbol] = float(loeuf) if loeuf is not None else None
            return _loeuf_cache[gene_symbol]
    except Exception as exc:
        log.debug("gnomAD LOEUF fetch error for %s: %s", gene_symbol, exc)

    _loeuf_cache[gene_symbol] = None
    return None


# ---------------------------------------------------------------------------
# Direction classifier
# ---------------------------------------------------------------------------

def classify_motif_direction(
    esm1v_score: Optional[float],
    consequence_score: Optional[float],
    loeuf: Optional[float],
) -> str:
    """Return 'gain_of_function', 'loss_of_function', 'likely_pathogenic', or 'neutral'.

    Args:
        esm1v_score: ESM-1v log-likelihood ratio (negative = destabilising)
        consequence_score: AlphaMissense score in [0,1] (high = pathogenic)
        loeuf: gnomAD LOEUF value (low = LoF intolerant)
    """
    am_high = (consequence_score is not None and consequence_score > 0.6)
    esm_destab = (esm1v_score is not None and esm1v_score < -2.0)
    lof_tolerant = (loeuf is None or loeuf > 0.5)   # unknown treated as tolerant

    if esm_destab and am_high:
        if lof_tolerant:
            return "loss_of_function"
        else:
            return "likely_pathogenic"
    elif am_high and not esm_destab:
        return "gain_of_function"
    else:
        return "neutral"


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def annotate_variant_directions(gene_ids: Optional[list[str]] = None) -> int:
    """Classify each DivergentMotif and store motif_direction.

    Args:
        gene_ids: Optional list of gene IDs to limit processing.

    Returns:
        Number of motifs annotated.
    """
    from db.models import DivergentMotif, Gene, Ortholog

    updated = 0
    with get_session() as session:
        q = (
            session.query(DivergentMotif)
            .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
            .join(Gene, Ortholog.gene_id == Gene.id)
        )
        if gene_ids:
            q = q.filter(Gene.id.in_(gene_ids))

        motifs = q.all()
        log.info("Classifying variant direction for %d motifs...", len(motifs))

        # Pre-fetch LOEUF values per gene symbol
        symbol_map: dict[str, Optional[str]] = {}
        gene_symbols_needed: set[str] = set()
        for motif in motifs:
            gs = motif.ortholog.gene.gene_symbol if motif.ortholog and motif.ortholog.gene else None
            if gs:
                gene_symbols_needed.add(gs)

        for gs in gene_symbols_needed:
            loeuf = _fetch_loeuf(gs)
            symbol_map[gs] = loeuf
            time.sleep(0.25)  # gnomAD rate limit

        for motif in motifs:
            gs = (motif.ortholog.gene.gene_symbol
                  if motif.ortholog and motif.ortholog.gene else None)
            loeuf = symbol_map.get(gs) if gs else None
            direction = classify_motif_direction(
                esm1v_score=getattr(motif, "esm1v_score", None),
                consequence_score=getattr(motif, "consequence_score", None),
                loeuf=loeuf,
            )
            motif.motif_direction = direction
            updated += 1

        session.commit()

    dist: dict[str, int] = {}
    for v in ["gain_of_function", "loss_of_function", "likely_pathogenic", "neutral"]:
        dist[v] = 0
    for motif in motifs:
        dist[getattr(motif, "motif_direction", "neutral")] = dist.get(
            getattr(motif, "motif_direction", "neutral"), 0) + 1
    log.info("Variant direction distribution: %s", dist)

    return updated


def run_variant_direction_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    return annotate_variant_directions(gene_ids)
