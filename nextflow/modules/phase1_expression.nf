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

workflow PHASE1_EXPRESSION {
    take:
    convergent_aa_done
    phylo_done

    main:
    expression_analysis(convergent_aa_done)
    bgee_expression(expression_analysis.out.expression_done)
    composite_score_phase1(bgee_expression.out.bgee_done, phylo_done)

    emit:
    scored = composite_score_phase1.out.scored
}
