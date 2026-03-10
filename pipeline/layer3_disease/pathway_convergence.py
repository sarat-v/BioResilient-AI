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
    pathway_score = hypergeometric(observed, background, n_candidates, n_genes)
                  + evolutionary_weight (sum of convergence_scores for observed genes)

Stores results in PathwayConvergence table (gene_id FK removed; pathway-level).
Exposed via GET /research/pathway-convergence.
"""

import logging
import math
from collections import defaultdict
from typing import Optional

import requests

from db.session import get_session

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


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_pathway_convergence(
    min_genes: int = 2,
    total_background_genes: int = 20000,
) -> list[dict]:
    """Compute enrichment of convergent candidates in biological pathways.

    Args:
        min_genes: Minimum candidate genes in a pathway to report it.
        total_background_genes: Approximate human protein-coding gene count.

    Returns:
        List of pathway dicts sorted by enrichment score descending:
        [{"pathway_id", "pathway_name", "gene_count", "candidate_count",
          "log_pvalue", "evolutionary_weight", "pathway_score", "gene_symbols"}, ...]
    """
    from db.models import CandidateScore, Gene

    with get_session() as session:
        # Fetch all scored candidates
        rows = (
            session.query(CandidateScore, Gene)
            .join(Gene, CandidateScore.gene_id == Gene.id)
            .filter(CandidateScore.composite_score > 0)
            .all()
        )
        if not rows:
            log.info("No candidate scores found; pathway convergence skipped.")
            return []

        n_candidates = len(rows)
        log.info("Computing pathway convergence over %d candidates...", n_candidates)

        # Build gene → pathways map (using stored pathway_ids on Gene)
        # Also use Reactome API for genes with human_gene_id
        pathway_genes: dict[str, list[tuple[str, str, float]]] = defaultdict(list)
        pathway_names: dict[str, str] = {}

        for cs, gene in rows:
            conv_score = cs.convergence_score or 0.0

            # Use pre-stored pathway IDs (from step 11 Reactome annotation)
            pathway_ids_str = getattr(gene, "pathway_ids", None)
            if pathway_ids_str:
                pathway_id_list = pathway_ids_str if isinstance(pathway_ids_str, list) else [pathway_ids_str]
                for pid in pathway_id_list:
                    pathway_genes[pid].append((gene.id, gene.gene_symbol or gene.id, conv_score))

            # Supplement with live Reactome lookup
            if gene.human_gene_id:
                pathways = _fetch_reactome_pathways_for_gene(gene.human_gene_id)
                for p in pathways[:10]:  # cap to avoid overwhelming memory
                    pid = p.get("stId") or p.get("dbId", "")
                    pname = p.get("displayName", pid)
                    if pid:
                        pathway_genes[pid].append((gene.id, gene.gene_symbol or gene.id, conv_score))
                        pathway_names[pid] = pname

    # Score each pathway
    results = []
    for pid, gene_entries in pathway_genes.items():
        # Deduplicate
        seen = {}
        for gid, gsym, cscore in gene_entries:
            if gid not in seen or cscore > seen[gid][1]:
                seen[gid] = (gsym, cscore)

        unique_genes = list(seen.items())
        k = len(unique_genes)
        if k < min_genes:
            continue

        K = _fetch_reactome_pathway_gene_count(pid)
        log_p = _log_hypergeometric(k, K, n_candidates, total_background_genes)
        evo_weight = sum(c for _, (_, c) in unique_genes)
        pathway_score = round(-log_p * 0.5 + evo_weight * 0.5, 4)

        results.append({
            "pathway_id": pid,
            "pathway_name": pathway_names.get(pid, pid),
            "gene_count": K,
            "candidate_count": k,
            "log_pvalue": round(log_p, 4),
            "evolutionary_weight": round(evo_weight, 4),
            "pathway_score": pathway_score,
            "gene_symbols": sorted({gsym for _, (gsym, _) in unique_genes}),
        })

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
