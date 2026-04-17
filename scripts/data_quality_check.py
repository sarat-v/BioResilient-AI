#!/usr/bin/env python3
"""Data quality checks for pipeline steps being SKIPPED in the rerun plan.

Runs against the live RDS database (must run from inside the VPC — use AWS Batch
or an EC2 bastion). Outputs a structured Markdown report with PASS/FAIL verdicts.

Usage (local with VPN / EC2 bastion):
    python scripts/data_quality_check.py

Usage (AWS Batch via Nextflow):
    nextflow run nextflow/main.nf -profile aws -entry data_quality_check
"""

import json
import os
import sys
import textwrap
from datetime import datetime
from typing import Any

import psycopg2
import psycopg2.extras

# ---------------------------------------------------------------------------
# DB connection
# ---------------------------------------------------------------------------

DATABASE_URL = os.environ.get("DATABASE_URL")
if not DATABASE_URL:
    rds_host = os.environ.get("RDS_HOST", "")
    rds_pass = os.environ.get("RDS_PASSWORD", "")
    if rds_host and rds_pass:
        DATABASE_URL = f"postgresql://bioresilient:{rds_pass}@{rds_host}:5432/bioresilient?sslmode=require"

if not DATABASE_URL:
    sys.exit("ERROR: DATABASE_URL not set. Export it or set RDS_HOST + RDS_PASSWORD.")


def get_conn():
    return psycopg2.connect(DATABASE_URL)


def run_query(sql: str) -> list[dict]:
    with get_conn() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(sql)
            return [dict(r) for r in cur.fetchall()]


def scalar(sql: str) -> Any:
    rows = run_query(sql)
    if not rows:
        return None
    first = rows[0]
    return list(first.values())[0]


# ---------------------------------------------------------------------------
# Report builder
# ---------------------------------------------------------------------------

PASS = "✅ PASS"
FAIL = "❌ FAIL"
WARN = "⚠️  WARN"
INFO = "ℹ️  INFO"

results: list[dict] = []
overall_fails: list[str] = []


def check(label: str, verdict: str, detail: str, step: str):
    results.append({"step": step, "label": label, "verdict": verdict, "detail": detail})
    if verdict == FAIL:
        overall_fails.append(f"[{step}] {label}: {detail}")
    print(f"  {verdict}  {label}: {detail}")


def section(title: str):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# STEP 3b — Orthologs
# ---------------------------------------------------------------------------

section("STEP 3b — Orthologs")

total_genes = scalar("SELECT COUNT(*) FROM gene")
genes_with_orthologs = scalar("SELECT COUNT(DISTINCT gene_id) FROM ortholog")
pct_coverage = round(100 * genes_with_orthologs / max(total_genes, 1), 1)
check(
    "Gene coverage",
    PASS if pct_coverage >= 80 else FAIL,
    f"{genes_with_orthologs}/{total_genes} genes have orthologs ({pct_coverage}%)",
    "step3b",
)

median_species = scalar("""
    SELECT PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY cnt)
    FROM (
        SELECT gene_id, COUNT(DISTINCT species_id) AS cnt
        FROM ortholog
        WHERE species_id NOT IN (SELECT id FROM species WHERE is_control = TRUE)
        GROUP BY gene_id
    ) t
""")
median_species = round(float(median_species), 1) if median_species else 0
check(
    "Median resilient species per gene",
    PASS if median_species >= 4 else (WARN if median_species >= 3 else FAIL),
    f"{median_species} species (need ≥4 for robust convergence)",
    "step3b",
)

low_identity_count = scalar("""
    SELECT COUNT(*) FROM ortholog WHERE sequence_identity_pct < 25
""")
total_orthologs = scalar("SELECT COUNT(*) FROM ortholog")
pct_low_id = round(100 * (low_identity_count or 0) / max(total_orthologs or 1, 1), 1)
median_identity = scalar("""
    SELECT PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY sequence_identity_pct)
    FROM ortholog WHERE sequence_identity_pct IS NOT NULL
""")
check(
    "Sequence identity distribution",
    PASS if pct_low_id < 15 else FAIL,
    f"Median identity={round(float(median_identity or 0), 1)}%, "
    f"Low-identity (<25%) rows: {low_identity_count} ({pct_low_id}%)",
    "step3b",
)

missing_identity = scalar("SELECT COUNT(*) FROM ortholog WHERE sequence_identity_pct IS NULL")
pct_missing = round(100 * (missing_identity or 0) / max(total_orthologs or 1, 1), 1)
check(
    "Missing identity values",
    PASS if pct_missing < 5 else WARN,
    f"{missing_identity} rows have NULL sequence_identity_pct ({pct_missing}%)",
    "step3b",
)

# ---------------------------------------------------------------------------
# STEP 3c — Nucleotide Regions
# ---------------------------------------------------------------------------

section("STEP 3c — Nucleotide Regions")

for region_type in ("cds", "promoter", "downstream"):
    cov = scalar(f"""
        SELECT COUNT(DISTINCT gene_id) FROM nucleotide_region
        WHERE region_type = '{region_type}' AND LENGTH(COALESCE(sequence,'')) > 0
    """)
    threshold = 70 if region_type == "cds" else 50
    pct = round(100 * (cov or 0) / max(total_genes, 1), 1)
    check(
        f"{region_type.upper()} sequence coverage",
        PASS if pct >= threshold else FAIL,
        f"{cov}/{total_genes} genes ({pct}%)",
        "step3c",
    )

cds_quality = run_query("""
    SELECT
        ROUND(AVG(percent_identity)::numeric, 1)  AS avg_pct_id,
        ROUND(AVG(gap_fraction)::numeric, 3)       AS avg_gap_frac,
        COUNT(*) FILTER (WHERE gap_fraction > 0.5) AS high_gap_count,
        COUNT(*) FILTER (WHERE percent_identity < 30) AS low_id_count
    FROM nucleotide_score WHERE region_type = 'cds'
""")
if cds_quality and cds_quality[0]["avg_pct_id"] is not None:
    q = cds_quality[0]
    avg_pct = float(q["avg_pct_id"] or 0)
    avg_gap = float(q["avg_gap_frac"] or 0)
    check(
        "CDS alignment quality",
        PASS if avg_pct >= 40 and avg_gap < 0.30 else FAIL,
        f"avg_pct_identity={avg_pct}%, avg_gap_frac={avg_gap:.3f}, "
        f"low_identity_rows={q['low_id_count']}, high_gap_rows={q['high_gap_count']}",
        "step3c",
    )
else:
    check("CDS alignment quality", FAIL, "No nucleotide_score rows found for CDS", "step3c")

# ---------------------------------------------------------------------------
# STEP 3d — PhyloConservation Scores
# ---------------------------------------------------------------------------

section("STEP 3d — PhyloConservation Scores")

cds_cov = scalar("SELECT COUNT(*) FILTER (WHERE cds_phylo_score IS NOT NULL) FROM phylo_conservation_score")
phylo_total = scalar("SELECT COUNT(*) FROM phylo_conservation_score")
pct_cds = round(100 * (cds_cov or 0) / max(phylo_total or total_genes or 1, 1), 1)
check(
    "CDS phyloP coverage",
    PASS if pct_cds >= 60 else FAIL,
    f"{cds_cov}/{phylo_total} genes with CDS phyloP score ({pct_cds}%)",
    "step3d",
)

phylo_dist = run_query("""
    SELECT
        ROUND(PERCENTILE_CONT(0.1) WITHIN GROUP (ORDER BY cds_phylo_score)::numeric, 3) AS p10,
        ROUND(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY cds_phylo_score)::numeric, 3) AS median,
        ROUND(PERCENTILE_CONT(0.9) WITHIN GROUP (ORDER BY cds_phylo_score)::numeric, 3) AS p90,
        COUNT(*) FILTER (WHERE cds_phylo_score < -3) AS highly_accelerated,
        COUNT(*) FILTER (WHERE cds_phylo_score > 10) AS implausibly_high
    FROM phylo_conservation_score WHERE cds_phylo_score IS NOT NULL
""")
if phylo_dist and phylo_dist[0]["median"] is not None:
    d = phylo_dist[0]
    median_cds = float(d["median"])
    check(
        "CDS phyloP distribution",
        PASS if median_cds > 0 and int(d["implausibly_high"] or 0) == 0 else FAIL,
        f"p10={d['p10']}, median={d['median']}, p90={d['p90']}, "
        f"highly_accelerated={d['highly_accelerated']}, implausibly_high={d['implausibly_high']}",
        "step3d",
    )
else:
    check("CDS phyloP distribution", FAIL, "No phylo_conservation_score data found", "step3d")

# ---------------------------------------------------------------------------
# STEP 4 — Divergent Motifs
# ---------------------------------------------------------------------------

section("STEP 4 — Divergent Motifs")

total_motifs = scalar("SELECT COUNT(*) FROM divergent_motif")
check("Total motifs exist", PASS if (total_motifs or 0) > 0 else FAIL,
      f"{total_motifs} motifs in DB", "step4")

genes_zero_motifs = scalar("""
    SELECT COUNT(*) FROM gene g
    WHERE NOT EXISTS (
        SELECT 1 FROM ortholog o JOIN divergent_motif dm ON dm.ortholog_id = o.id
        WHERE o.gene_id = g.id
    )
""")
pct_zero = round(100 * (genes_zero_motifs or 0) / max(total_genes, 1), 1)
check(
    "Genes with zero motifs",
    PASS if pct_zero < 15 else FAIL,
    f"{genes_zero_motifs} genes ({pct_zero}%) have no motifs — expect <15%",
    "step4",
)

div_dist = run_query("""
    SELECT
        ROUND(PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY divergence_score)::numeric, 3) AS p25,
        ROUND(PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY divergence_score)::numeric, 3) AS median,
        ROUND(PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY divergence_score)::numeric, 3) AS p75,
        COUNT(*) FILTER (WHERE divergence_score = 0)   AS exact_zero,
        COUNT(*) FILTER (WHERE divergence_score = 1.0) AS perfect_one,
        COUNT(*) FILTER (WHERE divergence_score > 0.9) AS very_high
    FROM divergent_motif WHERE divergence_score IS NOT NULL
""")
if div_dist and div_dist[0]["median"] is not None:
    d = div_dist[0]
    median_div = float(d["median"])
    exact_zero_pct = round(100 * int(d["exact_zero"] or 0) / max(total_motifs or 1, 1), 1)
    check(
        "Divergence score distribution",
        PASS if 0.10 <= median_div <= 0.65 and exact_zero_pct < 5 else FAIL,
        f"p25={d['p25']}, median={d['median']}, p75={d['p75']}, "
        f"exact_zero={d['exact_zero']} ({exact_zero_pct}%), perfect_one={d['perfect_one']}, very_high={d['very_high']}",
        "step4",
    )
else:
    check("Divergence score distribution", FAIL, "No divergence_score data", "step4")

# ---------------------------------------------------------------------------
# STEP 4b — Pfam Domains + AlphaMissense
# ---------------------------------------------------------------------------

section("STEP 4b — Pfam Domains + AlphaMissense")

am_scored = scalar("SELECT COUNT(*) FROM divergent_motif WHERE consequence_score IS NOT NULL")
am_pct = round(100 * (am_scored or 0) / max(total_motifs or 1, 1), 1)
check(
    "AlphaMissense coverage",
    PASS if am_pct >= 50 else FAIL,
    f"{am_scored}/{total_motifs} motifs ({am_pct}%) — expect ≥50%",
    "step4b",
)

flat_middle = scalar("""
    SELECT COUNT(*) FROM divergent_motif
    WHERE consequence_score BETWEEN 0.4 AND 0.6
""")
flat_pct = round(100 * (flat_middle or 0) / max(am_scored or 1, 1), 1)
check(
    "AlphaMissense not stuck at 0.5",
    PASS if flat_pct < 30 else FAIL,
    f"{flat_middle} motifs ({flat_pct}%) in flat 0.4–0.6 band — expect <30% (flat = scoring failed)",
    "step4b",
)

am_dist = run_query("""
    SELECT
        ROUND(PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY consequence_score)::numeric, 3) AS p25,
        ROUND(PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY consequence_score)::numeric, 3) AS median,
        ROUND(PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY consequence_score)::numeric, 3) AS p75,
        ROUND(AVG(consequence_score)::numeric, 3) AS mean
    FROM divergent_motif WHERE consequence_score IS NOT NULL
""")
if am_dist and am_dist[0]["median"] is not None:
    d = am_dist[0]
    check("AlphaMissense distribution", INFO,
          f"p25={d['p25']}, median={d['median']}, p75={d['p75']}, mean={d['mean']}", "step4b")

domain_pct = scalar("""
    SELECT ROUND(100.0 * COUNT(*) FILTER (WHERE domain_name IS NOT NULL) / COUNT(*), 1)
    FROM divergent_motif
""")
in_func_domain = scalar("SELECT COUNT(*) FROM divergent_motif WHERE in_functional_domain = TRUE")
check(
    "Domain annotation coverage",
    PASS if float(domain_pct or 0) > 0 else FAIL,
    f"{domain_pct}% motifs with domain_name; {in_func_domain} in functional domains",
    "step4b",
)

# ---------------------------------------------------------------------------
# STEP 4d — Variant Direction
# ---------------------------------------------------------------------------

section("STEP 4d — Variant Direction")

direction_rows = run_query("""
    SELECT
        motif_direction,
        COUNT(*) AS cnt,
        ROUND(100.0 * COUNT(*) / SUM(COUNT(*)) OVER (), 1) AS pct
    FROM divergent_motif
    WHERE motif_direction IS NOT NULL
    GROUP BY motif_direction
    ORDER BY cnt DESC
""")
direction_map = {r["motif_direction"]: {"cnt": r["cnt"], "pct": float(r["pct"])} for r in direction_rows}

classified = scalar("SELECT COUNT(*) FROM divergent_motif WHERE motif_direction IS NOT NULL")
classified_pct = round(100 * (classified or 0) / max(total_motifs or 1, 1), 1)
check(
    "Direction classification coverage",
    PASS if classified_pct >= 80 else FAIL,
    f"{classified}/{total_motifs} motifs ({classified_pct}%) classified — expect ≥80%",
    "step4d",
)

neutral_pct = direction_map.get("neutral", {}).get("pct", 0)
check(
    "Neutral % not stuck at 90 %",
    PASS if neutral_pct <= 88 else FAIL,
    f"neutral={neutral_pct}% — expect ≤88% (90%+ indicates old LOEUF bug still active)",
    "step4d",
)

func_shift = direction_map.get("functional_shift", {}).get("cnt", 0)
check(
    "functional_shift present (LOEUF fix applied)",
    PASS if int(func_shift or 0) > 0 else FAIL,
    f"functional_shift count={func_shift} — must be >0 after LOEUF fix",
    "step4d",
)

# Internal consistency: GoF AM score
gof_am = scalar("""
    SELECT AVG(consequence_score)
    FROM divergent_motif
    WHERE motif_direction = 'gain_of_function' AND consequence_score IS NOT NULL
""")
if gof_am is not None:
    check(
        "GoF motifs have high AlphaMissense",
        PASS if float(gof_am) > 0.55 else FAIL,
        f"GoF avg AM score = {round(float(gof_am), 3)} — expect >0.55",
        "step4d",
    )

lof_esm = scalar("""
    SELECT AVG(esm1v_score)
    FROM divergent_motif
    WHERE motif_direction = 'loss_of_function' AND esm1v_score IS NOT NULL
""")
if lof_esm is not None:
    check(
        "LoF motifs have destabilising ESM-1v",
        PASS if float(lof_esm) < -1.0 else FAIL,
        f"LoF avg ESM-1v = {round(float(lof_esm), 3)} — expect <-1.0",
        "step4d",
    )

for d, v in direction_map.items():
    print(f"    {d:25s}  {v['cnt']:>7,}  ({v['pct']}%)")

# ---------------------------------------------------------------------------
# STEP 5 — Phylogenetic Tree (file check via S3)
# ---------------------------------------------------------------------------

section("STEP 5 — Phylogenetic Tree (S3 artifact check)")

import subprocess

s3_bucket = os.environ.get("S3_BUCKET", "bioresilient-data")
treefile_key = "cache/species.treefile"
try:
    result = subprocess.run(
        ["aws", "s3", "ls", f"s3://{s3_bucket}/{treefile_key}"],
        capture_output=True, text=True, timeout=15,
    )
    tree_exists = result.returncode == 0
    check(
        "species.treefile exists in S3",
        PASS if tree_exists else FAIL,
        f"s3://{s3_bucket}/{treefile_key} — {'found' if tree_exists else 'NOT FOUND'}",
        "step5",
    )

    if tree_exists:
        # Download and count leaves
        dl = subprocess.run(
            ["aws", "s3", "cp", f"s3://{s3_bucket}/{treefile_key}", "-"],
            capture_output=True, text=True, timeout=15,
        )
        if dl.returncode == 0:
            from Bio import Phylo
            import io
            tree = Phylo.read(io.StringIO(dl.stdout), "newick")
            leaves = tree.get_terminals()
            db_species_count = scalar("SELECT COUNT(*) FROM species WHERE is_control = FALSE")
            check(
                "Tree leaf count matches DB species",
                PASS if len(leaves) == db_species_count else FAIL,
                f"Tree leaves={len(leaves)}, DB non-control species={db_species_count}",
                "step5",
            )
except Exception as e:
    check("species.treefile S3 check", WARN, f"Could not check: {e}", "step5")


# ---------------------------------------------------------------------------
# STEP 7b — True Convergent AA Substitutions
# ---------------------------------------------------------------------------

section("STEP 7b — True Convergent AA Substitutions")

conv_nonzero = scalar("SELECT COUNT(*) FROM divergent_motif WHERE convergent_aa_count > 0")
conv_max = scalar("SELECT MAX(convergent_aa_count) FROM divergent_motif")
check(
    "Convergent AA data written",
    PASS if (conv_nonzero or 0) > 0 else FAIL,
    f"{conv_nonzero} motifs with convergent_aa_count > 0, max={conv_max}",
    "step7b",
)

# Cross-check: motifs with high conv AA should be in genes with lower convergence_pval
cross_7b = run_query("""
    SELECT
        CASE WHEN dm.convergent_aa_count > 2 THEN 'high_conv' ELSE 'low_conv' END AS grp,
        ROUND(AVG(es.convergence_pval)::numeric, 4) AS avg_pval,
        COUNT(*) AS n
    FROM divergent_motif dm
    JOIN ortholog o ON o.id = dm.ortholog_id
    JOIN evolution_score es ON es.gene_id = o.gene_id
    WHERE es.convergence_pval IS NOT NULL
    GROUP BY 1
""")
cross_map = {r["grp"]: float(r["avg_pval"] or 1.0) for r in cross_7b}
high_pval = cross_map.get("high_conv", 1.0)
low_pval = cross_map.get("low_conv", 1.0)
check(
    "High conv-AA motifs in lower p-val genes",
    PASS if high_pval <= low_pval else FAIL,
    f"avg convergence_pval: high_conv_aa={high_pval:.4f}, low_conv_aa={low_pval:.4f} "
    f"(high_conv should be ≤ low_conv)",
    "step7b",
)

# ---------------------------------------------------------------------------
# STEP 8 — Functional Evidence (Expression)
# ---------------------------------------------------------------------------

section("STEP 8 — Functional Evidence")

tier12_genes = scalar("""
    SELECT COUNT(DISTINCT gene_id) FROM candidate_score WHERE tier IN ('Tier1','Tier2')
""")
expr_genes = scalar("SELECT COUNT(DISTINCT gene_id) FROM expression_result")
expr_cov_pct = round(100 * (expr_genes or 0) / max(tier12_genes or 1, 1), 1)
check(
    "Expression coverage vs Tier1/2 genes",
    PASS if expr_cov_pct >= 80 else FAIL,
    f"{expr_genes} genes with expression data, {tier12_genes} Tier1/2 genes ({expr_cov_pct}%)",
    "step8",
)

expr_dist = run_query("""
    SELECT
        ROUND(MIN(log2fc)::numeric, 3)                                            AS min,
        ROUND(PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY log2fc)::numeric, 3)  AS p25,
        ROUND(PERCENTILE_CONT(0.50) WITHIN GROUP (ORDER BY log2fc)::numeric, 3)  AS median,
        ROUND(PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY log2fc)::numeric, 3)  AS p75,
        ROUND(MAX(log2fc)::numeric, 3)                                            AS max,
        COUNT(*) FILTER (WHERE log2fc = 0)   AS zero_count,
        COUNT(*) FILTER (WHERE log2fc IS NULL) AS null_count,
        COUNT(*)                               AS total
    FROM expression_result
""")
if expr_dist and expr_dist[0]["median"] is not None:
    d = expr_dist[0]
    zero_pct = round(100 * int(d["zero_count"] or 0) / max(int(d["total"] or 1), 1), 1)
    check(
        "Expression scores not mostly zero",
        PASS if zero_pct < 30 else FAIL,
        f"min={d['min']}, p25={d['p25']}, median={d['median']}, p75={d['p75']}, max={d['max']}, "
        f"zero={d['zero_count']} ({zero_pct}%), null={d['null_count']}",
        "step8",
    )
else:
    check("Expression scores", FAIL, "No expression_result rows found", "step8")

tissue_count = scalar("SELECT COUNT(DISTINCT comparison) FROM expression_result")
check(
    "Tissue/condition diversity",
    PASS if (tissue_count or 0) >= 3 else FAIL,
    f"{tissue_count} distinct tissues/conditions in expression_result — expect ≥3",
    "step8",
)

cs_zero_pct = scalar("""
    SELECT ROUND(100.0 * COUNT(*) FILTER (WHERE expression_score = 0) / COUNT(*), 1)
    FROM candidate_score
""")
check(
    "CandidateScore expression_score not mostly zero",
    PASS if float(cs_zero_pct or 100) < 50 else FAIL,
    f"{cs_zero_pct}% of CandidateScore rows have expression_score=0 — expect <50%",
    "step8",
)

# ---------------------------------------------------------------------------
# CROSS-STEP INTEGRITY
# ---------------------------------------------------------------------------

section("Cross-Step Integrity")

orphan_candidates = scalar("""
    SELECT COUNT(*) FROM candidate_score cs
    LEFT JOIN evolution_score es ON es.gene_id = cs.gene_id
    WHERE es.gene_id IS NULL
""")
check(
    "No orphan candidate_score rows",
    PASS if (orphan_candidates or 0) == 0 else FAIL,
    f"{orphan_candidates} candidate_score rows with no matching evolution_score",
    "integrity",
)

tier1_zero_both = scalar("""
    SELECT COUNT(*) FROM candidate_score
    WHERE tier = 'Tier1'
      AND convergence_score = 0 AND selection_score = 0
""")
check(
    "No Tier1 genes with zero convergence+selection",
    PASS if (tier1_zero_both or 0) == 0 else FAIL,
    f"{tier1_zero_both} Tier1 genes have BOTH convergence_score=0 AND selection_score=0",
    "integrity",
)

conv_stddev = scalar("""
    SELECT ROUND(STDDEV(convergence_score)::numeric, 4)
    FROM candidate_score
""")
check(
    "Convergence score has real spread (stddev > 0.05)",
    PASS if float(conv_stddev or 0) > 0.05 else FAIL,
    f"stddev={conv_stddev} — if ≈0 all genes got same score (scoring collapse)",
    "integrity",
)

suspicious_tier1 = scalar("""
    SELECT COUNT(*) FROM candidate_score
    WHERE tier = 'Tier1' AND composite_score < 0.10
""")
check(
    "No Tier1 genes with composite_score < 0.10",
    PASS if (suspicious_tier1 or 0) == 0 else FAIL,
    f"{suspicious_tier1} Tier1 genes have composite_score < 0.10",
    "integrity",
)

# Top Tier1 genes (sanity display)
top_genes = run_query("""
    SELECT g.gene_symbol, cs.composite_score, cs.convergence_score, cs.selection_score, cs.expression_score
    FROM candidate_score cs
    JOIN gene g ON g.id = cs.gene_id
    WHERE cs.tier = 'Tier1'
    ORDER BY cs.composite_score DESC
    LIMIT 10
""")
print("\n  Top Tier1 genes:")
print(f"  {'gene':12s} {'composite':>10} {'convergence':>12} {'selection':>10} {'expression':>11}")
for r in top_genes:
    print(f"  {r['gene_symbol']:12s} {float(r['composite_score'] or 0):>10.4f} "
          f"{float(r['convergence_score'] or 0):>12.4f} "
          f"{float(r['selection_score'] or 0):>10.4f} "
          f"{float(r['expression_score'] or 0):>11.4f}")

# ---------------------------------------------------------------------------
# Final report
# ---------------------------------------------------------------------------

section("SUMMARY")

fail_count = sum(1 for r in results if r["verdict"] == FAIL)
warn_count = sum(1 for r in results if r["verdict"] == WARN)
pass_count = sum(1 for r in results if r["verdict"] == PASS)

print(f"\n  Total checks : {len(results)}")
print(f"  {PASS}     : {pass_count}")
print(f"  {WARN}     : {warn_count}")
print(f"  {FAIL}     : {fail_count}")

if overall_fails:
    print(f"\n  ⛔ Steps that must be MOVED TO RERUN:\n")
    steps_to_rerun = set()
    for f in overall_fails:
        print(f"    {f}")
        steps_to_rerun.add(f.split("]")[0].replace("[", ""))
    print(f"\n  Affected steps: {sorted(steps_to_rerun)}")
else:
    print(f"\n  ✅ All quality checks passed. All skipped steps are safe to keep as SKIP.")

# Write JSON for downstream consumption
report = {
    "generated_at": datetime.utcnow().isoformat() + "Z",
    "total_checks": len(results),
    "pass": pass_count,
    "warn": warn_count,
    "fail": fail_count,
    "steps_to_rerun": sorted(steps_to_rerun) if overall_fails else [],
    "details": results,
}
with open("data_quality_report.json", "w") as f:
    json.dump(report, f, indent=2)
print("\n  Full report written to data_quality_report.json")

# Upload to S3 if running inside AWS
s3_bucket = os.environ.get("S3_BUCKET")
if s3_bucket:
    try:
        import boto3
        s3 = boto3.client("s3", region_name=os.environ.get("AWS_DEFAULT_REGION", "ap-south-1"))
        s3.upload_file("data_quality_report.json", s3_bucket, "reports/data_quality_report.json")
        print(f"  Uploaded to s3://{s3_bucket}/reports/data_quality_report.json")
    except Exception as e:
        print(f"  S3 upload skipped: {e}")

if fail_count > 0:
    sys.exit(1)
