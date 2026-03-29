/*
 * Phase 1 — Expression & Scoring Module
 * Steps: 8, 8b, 9
 */

process expression_analysis {
    label 'clinical'
    cpus 4
    memory '8 GB'
    time '2h'

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
        --storage-root '${params.storage_root}' \
        --ncbi-api-key '${params.ncbi_api_key ?: ""}'
    """
}

process bgee_expression {
    label 'base'
    cpus 2
    memory '4 GB'
    time '1h'

    input:
    val expression_done

    output:
    val true, emit: bgee_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step8b \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process composite_score_phase1 {
    label 'base'
    cpus 2
    memory '4 GB'
    time '30m'

    input:
    val bgee_done
    val phylo_done

    output:
    val true, emit: scored

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step9 \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process generate_phase1_reports {
    label 'base'
    cpus 1
    memory '4 GB'
    time '30m'

    input:
    val scored

    output:
    val true, emit: reports_done

    script:
    def steps = "step5 step6 step6b step6c step7 step7b step8 step9"
    """
    export DATABASE_URL='${params.db_url}'
    export BIORESILIENT_STORAGE_ROOT='${params.storage_root}'
    for step in ${steps}; do
        echo "=== Generating report: \$step ==="
        python -m pipeline.step_reporter --step "\$step" \
            || echo "WARN: report for \$step failed (non-fatal)"
    done
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
        expression_analysis(convergent_aa_done)
        expression_done_ch = expression_analysis.out.expression_done
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
        composite_score_phase1(bgee_done_ch, phylo_done)
        scored_ch = composite_score_phase1.out.scored
    } else {
        scored_ch = Channel.value(true)
    }

    // Always generate reports after the last requested step completes
    generate_phase1_reports(scored_ch)

    emit:
    scored = scored_ch
}
