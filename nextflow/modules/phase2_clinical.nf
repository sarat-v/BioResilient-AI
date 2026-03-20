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
    // step10b and step11 can run in parallel after scoring
    alphagenome_regulatory(scored)
    disease_annotation(scored)

    // step11b and step11c can run in parallel after disease annotation
    rare_variants(disease_annotation.out.disease_done)
    literature_search(disease_annotation.out.disease_done)

    // step11d waits for 11b + 11c
    pathway_convergence(rare_variants.out.rare_done, literature_search.out.lit_done)

    // step12 waits for pathways + alphagenome
    druggability(pathway_convergence.out.pathways_done, alphagenome_regulatory.out.alphagenome_done)
    p2rank_pockets(druggability.out.drug_done)

    // step13 and step14 can run in parallel
    gene_therapy(p2rank_pockets.out.p2rank_done)
    safety_screen(p2rank_pockets.out.p2rank_done)
    depmap_gtex(safety_screen.out.safety_done)

    // step15 waits for therapy + depmap
    final_rescore(gene_therapy.out.therapy_done, depmap_gtex.out.depmap_done)

    emit:
    pipeline_done = final_rescore.out.pipeline_done
}
