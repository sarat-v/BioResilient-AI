#!/usr/bin/env python3
"""BioResilient AI — Full pipeline health check.

Queries every DB table, validates row counts and data quality for each step,
and outputs a comprehensive markdown report + JSON summary.

Designed to run inside AWS Batch (via Nextflow) or locally.

Usage:
    DATABASE_URL=... python scripts/pipeline_health_check.py
    DATABASE_URL=... python scripts/pipeline_health_check.py --output /tmp/health_check.md
"""
import json
import logging
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-5s %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("health_check")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def _query(conn, sql):
    cur = conn.cursor()
    cur.execute(sql)
    rows = cur.fetchall()
    cols = [d[0] for d in cur.description] if cur.description else []
    cur.close()
    return rows, cols


def run_health_check(output_path: str = None):
    import psycopg2

    db_url = os.environ.get("DATABASE_URL")
    if not db_url:
        log.error("DATABASE_URL not set")
        sys.exit(1)

    conn = psycopg2.connect(db_url)
    results = {}
    issues = []
    report_lines = []

    def section(title):
        report_lines.append(f"\n## {title}\n")

    def line(text):
        report_lines.append(text)

    def check(name, value, expected_min=None, expected_max=None, warn_msg=None):
        status = "OK"
        if expected_min is not None and value < expected_min:
            status = "FAIL"
            issues.append(f"{name}: {value} < expected minimum {expected_min}")
        if expected_max is not None and value > expected_max:
            status = "WARN"
            issues.append(f"{name}: {value} > expected maximum {expected_max}")
        icon = {"OK": "PASS", "WARN": "WARN", "FAIL": "FAIL"}[status]
        line(f"- [{icon}] **{name}**: {value:,}" if isinstance(value, int) else f"- [{icon}] **{name}**: {value}")
        results[name] = {"value": value, "status": status}
        return status

    report_lines.append(f"# BioResilient Pipeline Health Check")
    report_lines.append(f"\n**Generated**: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}\n")

    # ── STEP 1: Environment ────────────────────────────────────────────────
    section("Step 1: Environment")
    rows, _ = _query(conn, "SELECT 1")
    check("DB reachable", "yes" if rows else "no")

    # ── STEP 2: Species + Proteomes ────────────────────────────────────────
    section("Step 2: Species & Proteomes")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM species")
    check("Species count", rows[0][0], expected_min=16)
    rows, _ = _query(conn, "SELECT id, scientific_name, lineage_group, is_control FROM species ORDER BY id")
    line("\n| Species | Scientific Name | Lineage | Control |")
    line("|---------|----------------|---------|---------|")
    for r in rows:
        line(f"| {r[0]} | {r[1]} | {r[2] or '-'} | {'yes' if r[3] else 'no'} |")

    # ── STEP 3: OrthoFinder ────────────────────────────────────────────────
    section("Step 3: OrthoFinder & Orthologs")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM gene")
    check("Gene count (human)", rows[0][0], expected_min=10000)
    rows, _ = _query(conn, "SELECT COUNT(*) FROM ortholog")
    check("Ortholog count", rows[0][0], expected_min=100000)
    rows, _ = _query(conn, "SELECT COUNT(DISTINCT orthofinder_og) FROM ortholog WHERE orthofinder_og IS NOT NULL")
    check("Distinct OrthoGroups", rows[0][0], expected_min=5000)
    rows, _ = _query(conn, """
        SELECT species_id, COUNT(*) AS cnt
        FROM ortholog GROUP BY species_id ORDER BY cnt DESC
    """)
    line("\n| Species | Ortholog Count |")
    line("|---------|---------------|")
    for r in rows:
        line(f"| {r[0]} | {r[1]:,} |")

    # ── STEP 3c: Nucleotide Conservation ───────────────────────────────────
    section("Step 3c: Nucleotide Conservation")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM nucleotide_region")
    check("NucleotideRegion rows", rows[0][0])
    rows, _ = _query(conn, "SELECT COUNT(*) FROM nucleotide_score")
    check("NucleotideScore rows", rows[0][0])

    # ── STEP 4: Divergent Motifs ───────────────────────────────────────────
    section("Step 4: Divergent Motifs")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM divergent_motif")
    total_motifs = rows[0][0]
    check("Total motif rows", total_motifs, expected_min=100000)

    rows, _ = _query(conn, """
        SELECT COUNT(*) FROM (
            SELECT ortholog_id, start_pos, COUNT(*) AS cnt
            FROM divergent_motif
            GROUP BY ortholog_id, start_pos
            HAVING COUNT(*) > 1
        ) dups
    """)
    dup_count = rows[0][0]
    check("Duplicate motif positions", dup_count, expected_max=0)

    rows, _ = _query(conn, "SELECT COUNT(DISTINCT ortholog_id) FROM divergent_motif")
    check("Orthologs with motifs", rows[0][0], expected_min=50000)

    # Step 4b: domains
    section("Step 4b: Domain Annotations")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM divergent_motif WHERE in_functional_domain = true")
    check("Motifs in functional domains", rows[0][0])
    rows, _ = _query(conn, "SELECT COUNT(*) FROM divergent_motif WHERE domain_name IS NOT NULL")
    check("Motifs with domain name", rows[0][0])

    # Step 4c: ESM-1v
    section("Step 4c: ESM-1v Scores")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM divergent_motif WHERE esm1v_score IS NOT NULL")
    check("Motifs with ESM-1v score", rows[0][0])

    # Step 4d: Variant direction
    section("Step 4d: Variant Direction")
    rows, _ = _query(conn, "SELECT motif_direction, COUNT(*) FROM divergent_motif GROUP BY motif_direction ORDER BY COUNT(*) DESC")
    line("\n| Direction | Count |")
    line("|-----------|-------|")
    has_direction = 0
    for r in rows:
        line(f"| {r[0] or 'NULL'} | {r[1]:,} |")
        if r[0] is not None:
            has_direction += r[1]
    check("Motifs with direction assigned", has_direction)

    # ── STEP 3d: Phylo Conservation ───────────────────────────────────────
    section("Step 3d: Phylo Conservation")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM phylo_conservation_score")
    check("PhyloConservationScore rows", rows[0][0])

    # ── STEP 5: Species Tree ──────────────────────────────────────────────
    section("Step 5: Species Tree")
    line("- Species tree is a file artifact (species.treefile), not a DB table.")
    line("- Check S3: `aws s3 ls s3://bioresilient-data/cache/species.treefile`")

    # ── STEP 6: MEME/Selection ────────────────────────────────────────────
    section("Step 6: MEME / aBSREL Selection")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM evolution_score")
    check("EvolutionScore rows", rows[0][0], expected_min=1)
    rows, _ = _query(conn, "SELECT selection_model, COUNT(*) FROM evolution_score WHERE selection_model IS NOT NULL GROUP BY selection_model")
    if rows:
        line("\n| Model | Count |")
        line("|-------|-------|")
        for r in rows:
            line(f"| {r[0]} | {r[1]:,} |")
    else:
        line("\n**WARNING: No selection models found — step 6 likely did not produce results.**")
        issues.append("evolution_score has no selection_model data")

    rows, _ = _query(conn, "SELECT COUNT(*) FROM evolution_score WHERE dnds_ratio IS NOT NULL")
    check("Genes with dN/dS ratio", rows[0][0])

    # Step 6b: FEL/BUSTED
    section("Step 6b: FEL + BUSTED")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM evolution_score WHERE fel_sites IS NOT NULL")
    check("Genes with FEL sites", rows[0][0])
    rows, _ = _query(conn, "SELECT COUNT(*) FROM evolution_score WHERE busted_pvalue IS NOT NULL")
    check("Genes with BUSTED p-value", rows[0][0])

    # Step 6c: RELAX
    section("Step 6c: RELAX")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM evolution_score WHERE relax_k IS NOT NULL")
    check("Genes with RELAX k", rows[0][0])

    # ── STEP 7: Convergence ───────────────────────────────────────────────
    section("Step 7: Convergence Scoring")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM evolution_score WHERE convergence_count > 0")
    check("Genes with convergence_count > 0", rows[0][0])
    rows, _ = _query(conn, "SELECT COUNT(*) FROM evolution_score WHERE phylop_score IS NOT NULL")
    check("Genes with phyloP score", rows[0][0])

    # Step 7b: Convergent AA
    section("Step 7b: Convergent Amino Acids")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM divergent_motif WHERE convergent_aa_count > 0")
    check("Motifs with convergent AA > 0", rows[0][0])

    # ── STEP 8: Expression ────────────────────────────────────────────────
    section("Step 8: Expression Analysis")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM expression_result")
    check("ExpressionResult rows", rows[0][0])

    rows, _ = _query(conn, """
        SELECT COUNT(*) FROM (
            SELECT gene_id, geo_accession, comparison, COUNT(*) AS cnt
            FROM expression_result
            GROUP BY gene_id, geo_accession, comparison
            HAVING COUNT(*) > 1
        ) dups
    """)
    check("Duplicate expression results", rows[0][0], expected_max=0)

    rows, _ = _query(conn, "SELECT COUNT(DISTINCT geo_accession) FROM expression_result")
    check("Distinct GEO accessions", rows[0][0])

    # Step 8b: Bgee
    section("Step 8b: Bgee Expression")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM expression_result WHERE geo_accession LIKE 'bgee:%'")
    check("Bgee expression rows", rows[0][0])

    # ── STEP 9: Composite Scores ──────────────────────────────────────────
    section("Step 9: Composite Scoring")
    rows, _ = _query(conn, "SELECT COUNT(*) FROM candidate_score")
    check("CandidateScore rows", rows[0][0])
    rows, _ = _query(conn, "SELECT tier, COUNT(*) FROM candidate_score WHERE tier IS NOT NULL GROUP BY tier ORDER BY tier")
    if rows:
        line("\n| Tier | Count |")
        line("|------|-------|")
        for r in rows:
            line(f"| {r[0]} | {r[1]:,} |")
    else:
        line("\n**WARNING: No tier assignments found.**")
        issues.append("candidate_score has no tier assignments")

    rows, _ = _query(conn, """
        SELECT g.gene_symbol, cs.composite_score, cs.tier,
               cs.convergence_score, cs.selection_score, cs.expression_score
        FROM candidate_score cs
        JOIN gene g ON g.id = cs.gene_id
        WHERE cs.tier IN ('Tier1', 'Tier2')
        ORDER BY cs.composite_score DESC
        LIMIT 20
    """)
    if rows:
        line("\n### Top Tier1/Tier2 Candidates")
        line("| Gene | Composite | Tier | Convergence | Selection | Expression |")
        line("|------|-----------|------|-------------|-----------|------------|")
        for r in rows:
            line(f"| {r[0]} | {r[1]:.3f} | {r[2]} | {r[3]:.3f} | {r[4]:.3f} | {r[5]:.3f} |")

    # ── SUMMARY ───────────────────────────────────────────────────────────
    section("Summary")
    if issues:
        line(f"\n**{len(issues)} issue(s) found:**\n")
        for i, issue in enumerate(issues, 1):
            line(f"{i}. {issue}")
    else:
        line("\n**All checks passed.**")

    conn.close()

    report = "\n".join(report_lines)

    out = Path(output_path) if output_path else Path("/tmp/bioresilient/health_check.md")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report)
    log.info("Health check report written to %s", out)

    json_out = out.with_suffix(".json")
    json_out.write_text(json.dumps({"results": results, "issues": issues}, indent=2, default=str))
    log.info("Health check JSON written to %s", json_out)

    print(report)
    return len(issues)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", help="Output markdown file path")
    args = parser.parse_args()
    sys.exit(run_health_check(args.output))
