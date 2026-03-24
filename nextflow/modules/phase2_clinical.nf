/*
 * Phase 2 — Clinical Translation Module
 * Steps: 10b, 11, 11b, 11c, 11d, 12, 12b, 13, 14, 14b, 15
 *
 * All Phase 2 steps operate on Tier1/Tier2 genes from Phase 1 scoring.
 * Many annotation steps (11, 11b, 11c) are API-bound and can run in parallel.
 */

process alphagenome_regulatory {
    label 'base'
    cpus 2
    memory '4 GB'
    time '2h'

    input:
    val scored

    output:
    val true, emit: alphagenome_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step10b \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --ncbi-api-key '${params.ncbi_api_key ?: ""}'
    """
}

process disease_annotation {
    label 'base'
    cpus 4
    memory '8 GB'
    time '2h'

    input:
    val scored

    output:
    val true, emit: disease_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step11 \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process rare_variants {
    label 'base'
    cpus 2
    memory '4 GB'
    time '1h'

    input:
    val disease_done

    output:
    val true, emit: rare_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step11b \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process literature_search {
    label 'base'
    cpus 2
    memory '4 GB'
    time '1h'

    input:
    val disease_done

    output:
    val true, emit: lit_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step11c \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --ncbi-api-key '${params.ncbi_api_key ?: ""}'
    """
}

process pathway_convergence {
    label 'base'
    cpus 2
    memory '4 GB'
    time '30m'

    input:
    val rare_done
    val lit_done

    output:
    val true, emit: pathways_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step11d \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process druggability {
    label 'clinical'
    cpus 4
    memory '8 GB'
    time '2h'

    input:
    val pathways_done
    val alphagenome_done

    output:
    val true, emit: drug_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step12 \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process p2rank_pockets {
    label 'clinical'
    cpus 4
    memory '8 GB'
    time '1h'

    input:
    val drug_done

    output:
    val true, emit: p2rank_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step12b \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process gene_therapy {
    label 'base'
    cpus 2
    memory '4 GB'
    time '1h'

    input:
    val p2rank_done

    output:
    val true, emit: therapy_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step13 \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process safety_screen {
    label 'base'
    cpus 2
    memory '4 GB'
    time '1h'

    input:
    val p2rank_done

    output:
    val true, emit: safety_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step14 \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process depmap_gtex {
    label 'base'
    cpus 2
    memory '4 GB'
    time '1h'

    input:
    val safety_done

    output:
    val true, emit: depmap_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step14b \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process final_rescore {
    label 'base'
    cpus 2
    memory '4 GB'
    time '30m'

    input:
    val therapy_done
    val depmap_done

    output:
    val true, emit: pipeline_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step15 \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

workflow PHASE2_CLINICAL {
    take:
    scored

    main:
    def CLINICAL_STEPS = [
        'step10b','step11','step11b','step11c','step11d',
        'step12','step12b','step13','step14','step14b','step15'
    ]
    def fromStep = params.from_step ?: 'step1'
    def untilStep = params['until'] ?: 'step15'
    def fromIdx = CLINICAL_STEPS.indexOf(fromStep)
    def untilIdx = CLINICAL_STEPS.indexOf(untilStep)
    if (fromIdx < 0) fromIdx = 0
    if (untilIdx < 0) untilIdx = CLINICAL_STEPS.size() - 1

    // step10b and step11 can run in parallel after scoring
    if (fromIdx <= CLINICAL_STEPS.indexOf('step10b') && untilIdx >= CLINICAL_STEPS.indexOf('step10b')) {
        alphagenome_regulatory(scored)
        alphagenome_done_ch = alphagenome_regulatory.out.alphagenome_done
    } else {
        alphagenome_done_ch = Channel.value(true)
    }

    if (fromIdx <= CLINICAL_STEPS.indexOf('step11') && untilIdx >= CLINICAL_STEPS.indexOf('step11')) {
        disease_annotation(scored)
        disease_done_ch = disease_annotation.out.disease_done
    } else {
        disease_done_ch = Channel.value(true)
    }

    // step11b and step11c can run in parallel after disease annotation
    if (fromIdx <= CLINICAL_STEPS.indexOf('step11b') && untilIdx >= CLINICAL_STEPS.indexOf('step11b')) {
        rare_variants(disease_done_ch)
        rare_done_ch = rare_variants.out.rare_done
    } else {
        rare_done_ch = Channel.value(true)
    }

    if (fromIdx <= CLINICAL_STEPS.indexOf('step11c') && untilIdx >= CLINICAL_STEPS.indexOf('step11c')) {
        literature_search(disease_done_ch)
        lit_done_ch = literature_search.out.lit_done
    } else {
        lit_done_ch = Channel.value(true)
    }

    // step11d waits for 11b + 11c
    if (fromIdx <= CLINICAL_STEPS.indexOf('step11d') && untilIdx >= CLINICAL_STEPS.indexOf('step11d')) {
        pathway_convergence(rare_done_ch, lit_done_ch)
        pathways_done_ch = pathway_convergence.out.pathways_done
    } else {
        pathways_done_ch = Channel.value(true)
    }

    // step12 waits for pathways + alphagenome
    if (fromIdx <= CLINICAL_STEPS.indexOf('step12') && untilIdx >= CLINICAL_STEPS.indexOf('step12')) {
        druggability(pathways_done_ch, alphagenome_done_ch)
        drug_done_ch = druggability.out.drug_done
    } else {
        drug_done_ch = Channel.value(true)
    }

    if (fromIdx <= CLINICAL_STEPS.indexOf('step12b') && untilIdx >= CLINICAL_STEPS.indexOf('step12b')) {
        p2rank_pockets(drug_done_ch)
        p2rank_done_ch = p2rank_pockets.out.p2rank_done
    } else {
        p2rank_done_ch = Channel.value(true)
    }

    // step13 and step14 can run in parallel
    if (fromIdx <= CLINICAL_STEPS.indexOf('step13') && untilIdx >= CLINICAL_STEPS.indexOf('step13')) {
        gene_therapy(p2rank_done_ch)
        therapy_done_ch = gene_therapy.out.therapy_done
    } else {
        therapy_done_ch = Channel.value(true)
    }

    if (fromIdx <= CLINICAL_STEPS.indexOf('step14') && untilIdx >= CLINICAL_STEPS.indexOf('step14')) {
        safety_screen(p2rank_done_ch)
        safety_done_ch = safety_screen.out.safety_done
    } else {
        safety_done_ch = Channel.value(true)
    }

    if (fromIdx <= CLINICAL_STEPS.indexOf('step14b') && untilIdx >= CLINICAL_STEPS.indexOf('step14b')) {
        depmap_gtex(safety_done_ch)
        depmap_done_ch = depmap_gtex.out.depmap_done
    } else {
        depmap_done_ch = Channel.value(true)
    }

    // step15 waits for therapy + depmap
    if (fromIdx <= CLINICAL_STEPS.indexOf('step15') && untilIdx >= CLINICAL_STEPS.indexOf('step15')) {
        final_rescore(therapy_done_ch, depmap_done_ch)
        pipeline_done_ch = final_rescore.out.pipeline_done
    } else {
        pipeline_done_ch = Channel.value(true)
    }

    emit:
    pipeline_done = pipeline_done_ch
}
