#!/usr/bin/env python3
"""BioResilient AI — Benchmark Recall Evaluator.

Tests whether the pipeline rediscovers known evolutionary adaptations by
checking what percentile each benchmark gene scored in the full candidate
ranking.

Usage:
    python scripts/benchmark_recall.py
    python scripts/benchmark_recall.py --trait cancer_resistance
    python scripts/benchmark_recall.py --output benchmark_results.json

Success criterion (from ChatGPT review):
    >50% of benchmark genes appear in the top 5% of candidates.

Exit code:
    0 — benchmark passed (≥50% recall at top 5%)
    1 — benchmark failed
"""

import argparse
import json
import sys
from pathlib import Path

# Allow running from repo root or scripts/ directory
_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))


# ---------------------------------------------------------------------------
# Benchmark gene sets by trait
# ---------------------------------------------------------------------------

BENCHMARK_GENES: dict[str, list[dict]] = {
    "cancer_resistance": [
        {"symbol": "TP53",   "notes": "Tumour suppressor; duplicated in elephant (LIF6 pathway)"},
        {"symbol": "LIF6",   "notes": "Elephant-specific p53-like apoptosis inducer"},
        {"symbol": "ATM",    "notes": "DNA damage checkpoint kinase; mutated in many cancers"},
        {"symbol": "BRCA1",  "notes": "DNA repair; cancer predisposition gene"},
        {"symbol": "BRCA2",  "notes": "DNA repair; cancer predisposition gene"},
        {"symbol": "CDKN2A", "notes": "Cell cycle inhibitor (p16/p14); cancer suppressor"},
        {"symbol": "MDM2",   "notes": "p53 negative regulator; amplified in sarcomas"},
        {"symbol": "PTEN",   "notes": "Tumour suppressor; PI3K pathway"},
    ],
    "longevity": [
        {"symbol": "FOXO3",  "notes": "Forkhead TF; longevity association in multiple populations"},
        {"symbol": "IGF1R",  "notes": "Insulin/IGF1 receptor; dwarfism variants extend lifespan"},
        {"symbol": "SIRT6",  "notes": "NAD-dependent deacylase; DNA repair and metabolism"},
        {"symbol": "SIRT1",  "notes": "Caloric restriction mediator"},
        {"symbol": "MTOR",   "notes": "Nutrient sensor; rapamycin target extends lifespan"},
        {"symbol": "TERT",   "notes": "Telomerase reverse transcriptase"},
        {"symbol": "TERC",   "notes": "Telomerase RNA component"},
    ],
    "hypoxia_tolerance": [
        {"symbol": "EPAS1",  "notes": "HIF-2α; Tibetan high-altitude adaptation"},
        {"symbol": "HIF1A",  "notes": "HIF-1α; master hypoxia regulator"},
        {"symbol": "EGLN1",  "notes": "PHD2; EPAS1 regulator; Tibetan adaptation variant"},
        {"symbol": "PPARA",  "notes": "Fatty acid oxidation; hypoxia fuel switch"},
    ],
    "dna_repair": [
        {"symbol": "ERCC1",  "notes": "Nucleotide excision repair; bowhead whale evolution"},
        {"symbol": "XRCC5",  "notes": "NHEJ; Ku80 — DSB repair"},
        {"symbol": "XRCC6",  "notes": "NHEJ; Ku70 — DSB repair"},
        {"symbol": "RAD51",  "notes": "Homologous recombination — core recombinase"},
        {"symbol": "MLH1",   "notes": "Mismatch repair; microsatellite stability"},
    ],
    "viral_immunity": [
        {"symbol": "IFIH1",  "notes": "MDA5 — cytosolic RNA sensor; bat innate immune"},
        {"symbol": "STING1", "notes": "cGAS-STING innate immune pathway"},
        {"symbol": "IRF3",   "notes": "Interferon regulatory factor; antiviral"},
        {"symbol": "IRF7",   "notes": "Master regulator of type I IFN; bat-evolved"},
        {"symbol": "MAVS",   "notes": "Mitochondrial antiviral signaling"},
        {"symbol": "TLR7",   "notes": "Single-stranded RNA sensor; bat tolerance"},
        {"symbol": "ACE2",   "notes": "Coronavirus receptor; bat-adapted variants"},
    ],
}

# Genes expected in ANY trait run (core resilience biology)
UNIVERSAL_BENCHMARKS = [
    {"symbol": "TP53",   "trait": "cancer_resistance"},
    {"symbol": "FOXO3",  "trait": "longevity"},
    {"symbol": "EPAS1",  "trait": "hypoxia_tolerance"},
    {"symbol": "ERCC1",  "trait": "dna_repair"},
    {"symbol": "SIRT6",  "trait": "longevity"},
]


# ---------------------------------------------------------------------------
# Recall computation
# ---------------------------------------------------------------------------

def compute_recall(trait_id: str = "", top_pcts: list[float] = None) -> dict:
    """Query the DB and compute recall of benchmark genes.

    Args:
        trait_id: CandidateScore.trait_id to rank within; "" = default run.
        top_pcts: Percentile thresholds to evaluate (default: [1, 5, 10, 20]).

    Returns:
        Dict with per-gene results and summary recall metrics.
    """
    from db.models import CandidateScore, Gene
    from db.session import get_session

    if top_pcts is None:
        top_pcts = [1.0, 5.0, 10.0, 20.0]

    # Determine which benchmark set to use
    genes_to_check = BENCHMARK_GENES.get(trait_id, [])
    if not genes_to_check:
        # Default: use all benchmarks regardless of trait
        genes_to_check = [g for gs in BENCHMARK_GENES.values() for g in gs]

    with get_session() as session:
        all_scores = (
            session.query(CandidateScore, Gene)
            .join(Gene, CandidateScore.gene_id == Gene.id)
            .filter(CandidateScore.trait_id == trait_id)
            .order_by(CandidateScore.composite_score.desc())
            .all()
        )

        if not all_scores:
            print(f"  No CandidateScore rows found for trait_id={repr(trait_id)}.")
            print("  Run the pipeline first, then re-run this script.")
            return {}

        total = len(all_scores)
        symbol_to_rank: dict[str, int] = {}
        symbol_to_score: dict[str, float] = {}
        symbol_to_tier: dict[str, str] = {}

        for rank, (cs, gene) in enumerate(all_scores, 1):
            sym = (gene.gene_symbol or "").upper()
            if sym:
                symbol_to_rank[sym] = rank
                symbol_to_score[sym] = cs.composite_score or 0.0
                symbol_to_tier[sym] = cs.tier or "Tier3"

    # Evaluate each benchmark gene
    gene_results = []
    for entry in genes_to_check:
        sym = entry["symbol"].upper()
        if sym not in symbol_to_rank:
            gene_results.append({
                "symbol": sym,
                "notes": entry["notes"],
                "rank": None,
                "total": total,
                "percentile": None,
                "composite_score": None,
                "tier": None,
                "found": False,
            })
        else:
            rank = symbol_to_rank[sym]
            pct = rank / total * 100.0
            gene_results.append({
                "symbol": sym,
                "notes": entry["notes"],
                "rank": rank,
                "total": total,
                "percentile": round(pct, 2),
                "composite_score": round(symbol_to_score[sym], 4),
                "tier": symbol_to_tier[sym],
                "found": True,
            })

    # Summary recall at each threshold
    found_genes = [g for g in gene_results if g["found"]]
    recall_summary = {}
    for pct in top_pcts:
        n_in_top = sum(1 for g in found_genes if g["percentile"] is not None and g["percentile"] <= pct)
        recall_pct = n_in_top / len(genes_to_check) * 100 if genes_to_check else 0
        recall_summary[f"top_{pct:.0f}pct"] = {
            "n_genes_found": n_in_top,
            "n_benchmark_genes": len(genes_to_check),
            "recall_pct": round(recall_pct, 1),
            "passed": recall_pct >= 50.0,
        }

    return {
        "trait_id": trait_id,
        "total_candidates": total,
        "benchmark_genes": len(genes_to_check),
        "genes_in_candidates": len(found_genes),
        "recall": recall_summary,
        "gene_results": sorted(gene_results, key=lambda x: (x["percentile"] or 999, x["symbol"])),
    }


def print_results(results: dict) -> None:
    """Print a formatted recall report to stdout."""
    if not results:
        return

    print()
    print("=" * 65)
    print(f"  BioResilient Benchmark Recall — trait: {results['trait_id'] or 'all'}")
    print("=" * 65)
    print(f"  Total candidates scored : {results['total_candidates']}")
    print(f"  Benchmark genes checked : {results['benchmark_genes']}")
    print(f"  Benchmark genes present : {results['genes_in_candidates']}")
    print()

    # Per-threshold summary
    print("  Recall by rank threshold:")
    print(f"  {'Threshold':<12} {'Found':<8} {'Total':<8} {'Recall':<10} {'Status'}")
    print("  " + "-" * 50)
    for key, r in results["recall"].items():
        status = "✅ PASS" if r["passed"] else "❌ FAIL"
        print(f"  {key:<12} {r['n_genes_found']:<8} {r['n_benchmark_genes']:<8} "
              f"{r['recall_pct']:>5.1f}%    {status}")

    print()
    print("  Per-gene breakdown:")
    print(f"  {'Gene':<10} {'Rank':>6} {'Top%':>7} {'Score':>7} {'Tier':<8} Notes")
    print("  " + "-" * 70)
    for g in results["gene_results"]:
        if g["found"]:
            pct_str = f"{g['percentile']:>6.1f}%"
            rank_str = f"{g['rank']:>6,}"
            score_str = f"{g['composite_score']:>7.4f}"
        else:
            pct_str = "not found"
            rank_str = "    —"
            score_str = "      —"
        tier = g["tier"] or "—"
        notes = g["notes"][:40] if g["notes"] else ""
        print(f"  {g['symbol']:<10} {rank_str} {pct_str} {score_str} {tier:<8} {notes}")

    print()
    # Overall verdict
    top5_recall = results["recall"].get("top_5pct", {}).get("recall_pct", 0)
    if top5_recall >= 50:
        print(f"  ✅ BENCHMARK PASSED: {top5_recall:.1f}% of known genes in top 5% (target ≥50%)")
    else:
        print(f"  ❌ BENCHMARK FAILED: {top5_recall:.1f}% of known genes in top 5% (target ≥50%)")
        print()
        print("  Diagnostic: check which sub-score is lowest for missing genes:")
        print("    python -c \"from db.models import CandidateScore, Gene; ...\"")
        print("  Or use: GET /candidates/<gene_id>/scores via the API")
    print("=" * 65)
    print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate pipeline recall against known evolutionary adaptation genes."
    )
    parser.add_argument(
        "--trait",
        default="",
        choices=["", "cancer_resistance", "longevity", "hypoxia_tolerance",
                 "dna_repair", "viral_immunity"],
        help="Which trait_id to evaluate (default: '' = default pipeline run)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Optional JSON output file path",
    )
    parser.add_argument(
        "--top-pct",
        nargs="+",
        type=float,
        default=[1.0, 5.0, 10.0, 20.0],
        metavar="PCT",
        help="Rank percentile thresholds to report (default: 1 5 10 20)",
    )
    args = parser.parse_args()

    print(f"Evaluating benchmark recall for trait: {repr(args.trait) or 'default'}")

    results = compute_recall(trait_id=args.trait, top_pcts=args.top_pct)
    if not results:
        return 1

    print_results(results)

    if args.output:
        out_path = Path(args.output)
        out_path.write_text(json.dumps(results, indent=2))
        print(f"Results saved to: {out_path}")

    # Exit 0 if top-5% recall ≥ 50%, else exit 1
    top5 = results.get("recall", {}).get("top_5pct", {}).get("passed", False)
    return 0 if top5 else 1


if __name__ == "__main__":
    sys.exit(main())
