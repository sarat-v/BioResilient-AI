/*
 * Phase 1 — Functional Evidence & Scoring Module
 * Steps: 8, 8b (no-op), 9
 *
 * Step 8 replaces GEO/DESeq2 + Bgee with three well-curated human databases:
 *   - Open Targets Platform  (gene-disease association scores, any phenotype)
 *   - GTEx v10               (tissue expression in phenotype-relevant tissues)
 *   - DepMap CRISPR          (selective essentiality, cancer/dna_repair only)
 * Config: config/functional_evidence_config.json
 */

process functional_evidence {
    label 'base'
    cpus 2
    memory '8 GB'
    time '3h'

    input:
    val convergent_aa_done

    output:
    val true, emit: expression_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step8 \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step8 || true
    """
}

process bgee_expression {
    label 'base'
    cpus 1
    memory '2 GB'
    time '5m'

    input:
    val expression_done

    output:
    val true, emit: bgee_done

    script:
    """
    # Step 8b is now a no-op — functional evidence is unified in step 8.
    echo "Step 8b: pass-through (functional evidence unified in step 8)"
    python -m scripts.nf_wrappers.run_step \
        --step step8b \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' || true
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step8b || true
    """
}

// Pre-step B: remove orphan candidate_score rows before step9 rescores.
// Rows that have no matching evolution_score are stale genes from earlier
// runs (e.g. genes removed from the target set). If left in, they inflate
// tier counts and cause unexpected joins during composite scoring.
// Idempotent: safe to run even when there are no orphans.
process purge_orphan_scores {
    label 'base'
    cpus 1
    memory '2 GB'
    time '10m'

    input:
    val bgee_done

    output:
    val true, emit: scores_clean

    script:
    """
    export DATABASE_URL='${params.db_url}'
    python3 << 'PYEOF'
import psycopg2, os

conn = psycopg2.connect(os.environ['DATABASE_URL'])
cur  = conn.cursor()

cur.execute(
    'SELECT COUNT(*) FROM candidate_score cs'
    ' WHERE NOT EXISTS (SELECT 1 FROM evolution_score es WHERE es.gene_id = cs.gene_id)'
)
orphan_count = cur.fetchone()[0]

if orphan_count > 0:
    cur.execute(
        'DELETE FROM candidate_score'
        ' WHERE gene_id NOT IN (SELECT gene_id FROM evolution_score)'
    )
    conn.commit()
    print(f'Pre-step B: deleted {orphan_count} orphan candidate_score rows')
else:
    print('Pre-step B: no orphan candidate_score rows found (already clean)')

conn.close()
PYEOF
    """
}

process composite_score_phase1 {
    label 'base'
    cpus 4
    memory '8 GB'
    time '1h'

    input:
    val scores_clean
    val phylo_done

    output:
    val true, emit: scored

    script:
    """
    echo "rank_product_v3_qval_composite"
    python -m scripts.nf_wrappers.run_step \
        --step step9 \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step9 || true
    """
}

// Runs after every Phase 1 step completes.
// Generates per-step MD/JSON reports AND a comprehensive health check.
process generate_phase1_reports {
    label 'base'
    cpus 1
    memory '4 GB'
    time '30m'

    input:
    val scored

    output:
    val true, emit: reports_done
    path 'health_check.md', emit: health_report, optional: true
    path 'health_check.json', emit: health_json, optional: true

    script:
    def allSteps = "step1 step2 step3 step3b step3c step4 step4b step4c step4d step3d step5 step6 step7 step7b step8 step8b step9"
    """
    export DATABASE_URL='${params.db_url}'
    export BIORESILIENT_STORAGE_ROOT='${params.storage_root}'

    echo "=== Per-step reports ==="
    for step in ${allSteps}; do
        echo "--- \$step ---"
        python -m pipeline.step_reporter --step "\$step" \
            || echo "WARN: report for \$step skipped (step may not have run)"
    done

    echo ""
    echo "=== Pipeline Health Check ==="
    python /app/scripts/pipeline_health_check.py --output health_check.md
    echo "=== Reports complete ==="
    """
}

workflow PHASE1_EXPRESSION {
    take:
    convergent_aa_done
    phylo_done

    main:
    def EXPR_STEPS = ['step8','step8b','step9']
    def fromStep = params.from_step ?: 'step1'
    def untilStep = params['until'] ?: 'step9'
    def fromIdx = EXPR_STEPS.indexOf(fromStep)
    def untilIdx = EXPR_STEPS.indexOf(untilStep)
    if (fromIdx < 0) fromIdx = 0
    if (untilIdx < 0) untilIdx = EXPR_STEPS.size() - 1

    if (fromIdx <= EXPR_STEPS.indexOf('step8') && untilIdx >= EXPR_STEPS.indexOf('step8')) {
        functional_evidence(convergent_aa_done)
        expression_done_ch = functional_evidence.out.expression_done
    } else {
        expression_done_ch = Channel.value(true)
    }

    if (fromIdx <= EXPR_STEPS.indexOf('step8b') && untilIdx >= EXPR_STEPS.indexOf('step8b')) {
        bgee_expression(expression_done_ch)
        bgee_done_ch = bgee_expression.out.bgee_done
    } else {
        bgee_done_ch = Channel.value(true)
    }

    if (fromIdx <= EXPR_STEPS.indexOf('step9') && untilIdx >= EXPR_STEPS.indexOf('step9')) {
        purge_orphan_scores(bgee_done_ch)
        composite_score_phase1(purge_orphan_scores.out.scores_clean, phylo_done)
        scored_ch = composite_score_phase1.out.scored
    } else {
        scored_ch = Channel.value(true)
    }

    // Final consolidated summary runs after the last completed step
    generate_phase1_reports(scored_ch)

    emit:
    scored = scored_ch
}
