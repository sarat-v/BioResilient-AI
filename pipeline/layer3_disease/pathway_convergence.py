"""Pathway-level convergence scoring.

Maps convergent candidate genes to biological pathways (Reactome / KEGG via
pathway_ids stored on Gene) and computes a weighted enrichment score for each
pathway.

The signal: if multiple independent species independently evolve genes that
all belong to the same pathway, that pathway is likely under selection for
the trait — even if no single gene reaches the convergence threshold alone.

Algorithm:
  For each pathway P:
    observed = set of candidate genes mapped to P
    background = total genes in P (Reactome API)
    log_pvalue   = log10 hypergeometric right-tail p-value
    adjusted_pvalue = Benjamini-Hochberg FDR correction across all tested pathways
                      (Benjamini & Hochberg 1995, J Royal Stat Soc B 57:289-300)
    pathway_score = combined ranking: -log_pvalue * 0.5 + evo_weight * 0.5

Note on multiple-testing: the number of pathways tested is small (O(10-100)
for Tier1/Tier2 candidates), so the BH correction has moderate impact.
Pathways with adjusted_pvalue < 0.05 should be treated as statistically enriched.
Pathways with adjusted_pvalue < 0.20 are reportable as trends.

Stores results in PathwayConvergence table (gene_id FK removed; pathway-level).
Exposed via GET /research/pathway-convergence.
"""

import logging
import math
from collections import defaultdict
from typing import Optional

import requests

from db.session import get_session
from pipeline.stats import apply_bh_correction

log = logging.getLogger(__name__)

REACTOME_PATHWAY_API = "https://reactome.org/ContentService/data/pathways/low/diagram/entity/{gene_id}/allForms"
REACTOME_PATHWAY_GENES_API = "https://reactome.org/ContentService/data/pathway/{pathway_id}/containedEvents"


# ---------------------------------------------------------------------------
# Hypergeometric enrichment (log p-value for stability)
# ---------------------------------------------------------------------------

def _log_hypergeometric(k: int, K: int, n: int, N: int) -> float:
    """Log10 of the hypergeometric p-value (right tail, k or more successes).

    k = observed pathway genes in candidates
    K = total pathway gene count (background)
    n = total candidates
    N = total background genes
    """
    if K == 0 or N == 0 or n == 0:
        return 0.0

    # Use cumulative log-probability; fall back to approximation for large numbers
    try:
        from math import lgamma, log

        def log_comb(a: int, b: int) -> float:
            if b < 0 or b > a:
                return float("-inf")
            return lgamma(a + 1) - lgamma(b + 1) - lgamma(a - b + 1)

        log_total = log_comb(N, n)
        if log_total == float("-inf"):
            return 0.0

        log_pvalue = float("-inf")
        for i in range(k, min(K, n) + 1):
            lp = log_comb(K, i) + log_comb(N - K, n - i) - log_total
            if lp != float("-inf"):
                if log_pvalue == float("-inf"):
                    log_pvalue = lp
                else:
                    log_pvalue = lp + math.log1p(math.exp(log_pvalue - lp))

        return max(log_pvalue / math.log(10), -300.0) if log_pvalue != float("-inf") else 0.0
    except Exception:
        return 0.0


# ---------------------------------------------------------------------------
# Reactome gene → pathway mapping
# ---------------------------------------------------------------------------

def _fetch_reactome_pathways_for_gene(human_gene_id: str, timeout: int = 10) -> list[dict]:
    """Return list of Reactome pathways containing a human gene (by Ensembl / UniProt ID)."""
    url = REACTOME_PATHWAY_API.format(gene_id=human_gene_id)
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code == 200:
            return r.json()
    except Exception as exc:
        log.debug("Reactome API error for %s: %s", human_gene_id, exc)
    return []


def _fetch_reactome_pathway_gene_count(pathway_id: str, timeout: int = 10) -> int:
    """Return the number of human genes annotated to a Reactome pathway."""
    url = f"https://reactome.org/ContentService/data/pathway/{pathway_id}/participants/count"
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code == 200:
            data = r.json()
            return int(data.get("proteins", data.get("total", 50)))
    except Exception:
        pass
    return 50   # conservative fallback


def _fetch_reactome_pathway_name(pathway_id: str, timeout: int = 10) -> str:
    """Return the human-readable displayName for a Reactome pathway ID."""
    url = f"https://reactome.org/ContentService/data/query/{pathway_id}"
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code == 200:
            data = r.json()
            # ContentService returns either an object or a list
            if isinstance(data, list) and data:
                data = data[0]
            name = data.get("displayName") or data.get("name") or ""
            if name:
                return name
    except Exception:
        pass
    return pathway_id   # fallback to raw ID


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_pathway_convergence(
    min_genes: int = 2,
    total_background_genes: int = 20000,
) -> list[dict]:
    """Compute enrichment of convergent candidates in biological pathways.

    Scoped to Tier1/Tier2 candidates only (typically 40-60 genes).
    Uses pre-stored pathway_ids from the Gene table (populated by step 11)
    to avoid O(N) Reactome API calls per gene.

    Args:
        min_genes: Minimum candidate genes in a pathway to report it.
        total_background_genes: Approximate human protein-coding gene count.

    Returns:
        List of pathway dicts sorted by enrichment score descending.
    """
    from db.models import CandidateScore, Gene

    with get_session() as session:
        # Scope to Tier1/Tier2 only — avoids processing all 12k scored genes
        rows = (
            session.query(CandidateScore, Gene)
            .join(Gene, CandidateScore.gene_id == Gene.id)
            .filter(CandidateScore.tier.in_(["Tier1", "Tier2"]))
            .all()
        )
        if not rows:
            log.info("No Tier1/Tier2 candidates found; pathway convergence skipped.")
            return []

        n_candidates = len(rows)
        log.info("Computing pathway convergence over %d Tier1/Tier2 candidates...", n_candidates)

        # Build gene → pathways map using pre-stored pathway_ids only.
        # The live Reactome gene lookup is skipped here because:
        #  (a) pathway_ids are already populated by step 11 (UniProt/Reactome),
        #  (b) gene.human_gene_id uses OMA format (human|AKT1_HUMAN) which Reactome
        #      API doesn't accept directly.
        pathway_genes: dict[str, list[tuple[str, str, float]]] = defaultdict(list)
        pathway_names: dict[str, str] = {}

        for cs, gene in rows:
            conv_score = cs.convergence_score or 0.0

            pathway_ids_raw = getattr(gene, "pathway_ids", None)
            if not pathway_ids_raw:
                continue
            pid_list = pathway_ids_raw if isinstance(pathway_ids_raw, list) else [pathway_ids_raw]
            clean_symbol = gene.gene_symbol or gene.id
            # Strip _HUMAN suffix for display
            if clean_symbol.endswith("_HUMAN"):
                clean_symbol = clean_symbol[:-6]
            for pid in pid_list:
                if pid:
                    pathway_genes[pid].append((gene.id, clean_symbol, conv_score))

    # Score each pathway with in-memory cache for Reactome gene-count lookups
    _gene_count_cache: dict[str, int] = {}
    results = []
    for pid, gene_entries in pathway_genes.items():
        # Deduplicate by gene_id, keeping highest convergence score
        seen: dict[str, tuple[str, float]] = {}
        for gid, gsym, cscore in gene_entries:
            if gid not in seen or cscore > seen[gid][1]:
                seen[gid] = (gsym, cscore)

        unique_genes = list(seen.items())
        k = len(unique_genes)
        if k < min_genes:
            continue

        if pid not in _gene_count_cache:
            import time as _time
            _time.sleep(0.05)  # gentle rate limit (~20 req/s)
            _gene_count_cache[pid] = _fetch_reactome_pathway_gene_count(pid)
            # Populate human-readable name on first encounter
            if pid not in pathway_names:
                pathway_names[pid] = _fetch_reactome_pathway_name(pid)
        K = _gene_count_cache[pid]
        log_p = _log_hypergeometric(k, K, n_candidates, total_background_genes)
        evo_weight = sum(c for _, (_, c) in unique_genes)
        pathway_score = round(-log_p * 0.5 + evo_weight * 0.5, 4)

        results.append({
            "pathway_id": pid,
            "pathway_name": pathway_names.get(pid, pid),
            "gene_count": K,
            "candidate_count": k,
            "log_pvalue": round(log_p, 4),
            "adjusted_pvalue": None,   # filled in by BH step below
            "evolutionary_weight": round(evo_weight, 4),
            "pathway_score": pathway_score,
            "gene_symbols": sorted({gsym for _, (gsym, _) in unique_genes}),
        })

    # Benjamini-Hochberg FDR correction across all tested pathways.
    # Convert log10 p-values back to linear scale for BH; clamp to [0, 1].
    # BH reference: Benjamini & Hochberg (1995) J Royal Stat Soc B 57:289-300.
    if results:
        raw_pvalues = [
            min(1.0, max(0.0, 10 ** r["log_pvalue"])) for r in results
        ]
        adjusted = apply_bh_correction(raw_pvalues)
        for r, adj in zip(results, adjusted):
            r["adjusted_pvalue"] = round(adj, 6) if adj is not None else None

        n_sig = sum(1 for r in results if r["adjusted_pvalue"] is not None and r["adjusted_pvalue"] < 0.05)
        log.info("  Pathway BH-FDR: %d / %d pathways at adjusted p < 0.05.", n_sig, len(results))

    results.sort(key=lambda x: x["pathway_score"], reverse=True)
    log.info("Pathway convergence: %d pathways scored (≥%d candidates).", len(results), min_genes)
    return results


# ---------------------------------------------------------------------------
# DB persistence
# ---------------------------------------------------------------------------

def persist_pathway_convergence(results: list[dict]) -> int:
    """Store pathway convergence results in PathwayConvergence table."""
    from db.models import PathwayConvergence

    with get_session() as session:
        session.query(PathwayConvergence).delete()
        for r in results:
            row = PathwayConvergence(
                pathway_id=r["pathway_id"],
                pathway_name=r["pathway_name"],
                gene_count=r["gene_count"],
                candidate_count=r["candidate_count"],
                log_pvalue=r["log_pvalue"],
                adjusted_pvalue=r.get("adjusted_pvalue"),
                evolutionary_weight=r["evolutionary_weight"],
                pathway_score=r["pathway_score"],
                gene_symbols=r["gene_symbols"],
            )
            session.add(row)
        session.commit()

    log.info("Stored %d pathway convergence rows.", len(results))
    return len(results)


def run_pathway_convergence_pipeline() -> None:
    """Entry point called from orchestrator."""
    results = compute_pathway_convergence()
    if results:
        persist_pathway_convergence(results)
        log.info("Top pathways: %s",
                 [r["pathway_name"] for r in results[:5]])
