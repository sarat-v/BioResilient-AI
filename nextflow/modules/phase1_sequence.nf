/*
 * Phase 1 — Sequence Analysis Module
 * Steps: 1, 2, 3, 3b, 3c, 4, 4b, 4c, 4d
 */

process validate_environment {
    label 'base'
    cpus 1
    memory '2 GB'

    output:
    val true, emit: validated

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step1 \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process download_proteomes {
    label 'base'
    cpus 2
    memory '4 GB'

    input:
    val ready

    output:
    path 'proteomes/', emit: proteomes_dir

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step2 \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --ncbi-api-key '${params.ncbi_api_key ?: ""}'
    ln -sf \$(python -c "from pipeline.config import get_local_storage_root; print(get_local_storage_root())")/proteomes proteomes
    """
}

process run_orthofinder {
    label 'orthofinder'
    cpus params.orthofinder_cpus
    memory params.orthofinder_memory
    time '4h'

    input:
    path proteomes

    output:
    path 'orthofinder_results/', emit: results_dir

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step3 \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    STORAGE=\$(python -c "from pipeline.config import get_local_storage_root; print(get_local_storage_root())")
    ln -sf \$STORAGE/orthofinder_out orthofinder_results
    """
}

process load_orthologs {
    label 'base'
    cpus 2
    memory '8 GB'
    time '30m'

    input:
    path results_dir

    output:
    val true, emit: orthologs_loaded

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step3b \
        --input-dir '${results_dir}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process nucleotide_conservation {
    label 'align'
    cpus params.align_cpus
    memory params.align_memory
    time '4h'

    input:
    val orthologs_loaded

    output:
    val true, emit: nuc_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step3c \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --ncbi-api-key '${params.ncbi_api_key ?: ""}'
    """
}

process align_and_divergence {
    label 'align'
    cpus params.align_cpus
    memory params.align_memory
    time '6h'

    input:
    val orthologs_loaded

    output:
    path 'aligned_orthogroups.pkl', emit: aligned_pkl
    val true, emit: motifs_loaded

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step4 \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --output-dir '.'
    """
}

process domain_and_consequence {
    label 'base'
    cpus 4
    memory '8 GB'
    time '2h'

    input:
    val motifs_loaded

    output:
    val true, emit: domains_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step4b \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process esm1v_scoring {
    label 'esm'
    cpus params.esm_cpus
    memory params.esm_memory
    time '2h'

    input:
    val motifs_loaded

    output:
    val true, emit: esm_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step4c \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process variant_direction {
    label 'base'
    cpus 2
    memory '4 GB'
    time '30m'
    val domains_done
    val esm_done

    output:
    val true, emit: directions_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step4d \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

workflow PHASE1_SEQUENCE {
    take:
    start_signal

    main:
    validate_environment()
    download_proteomes(validate_environment.out.validated)
    run_orthofinder(download_proteomes.out.proteomes_dir)
    load_orthologs(run_orthofinder.out.results_dir)

    // step3c and step4 can run in parallel after orthologs are loaded
    nucleotide_conservation(load_orthologs.out.orthologs_loaded)
    align_and_divergence(load_orthologs.out.orthologs_loaded)

    // step4b and step4c can run in parallel after motifs exist
    domain_and_consequence(align_and_divergence.out.motifs_loaded)
    esm1v_scoring(align_and_divergence.out.motifs_loaded)

    // step4d waits for both 4b and 4c
    variant_direction(domain_and_consequence.out.domains_done, esm1v_scoring.out.esm_done)

    emit:
    aligned_pkl     = align_and_divergence.out.aligned_pkl
    nuc_done        = nucleotide_conservation.out.nuc_done
    directions_done = variant_direction.out.directions_done
}
