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

    input:
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
    // Step ordering within this phase — drives skip logic when resuming
    def SEQ_STEPS = ['step1','step2','step3','step3b','step3c',
                     'step4','step4b','step4c','step4d']
    def fromStep = params.from_step ?: 'step1'
    def untilStep = params['until'] ?: 'step4d'
    // fromIdx == -1 means from_step is a later phase; treat as 0 (run nothing here — main.nf handles it)
    def fromIdx = SEQ_STEPS.indexOf(fromStep)
    def untilIdx = SEQ_STEPS.indexOf(untilStep)
    if (fromIdx < 0) fromIdx = 0
    if (untilIdx < 0) untilIdx = SEQ_STEPS.size() - 1

    // ── Steps 1–3b: download + OrthoFinder ───────────────────────────────
    // Skip if resuming from step3c or later (orthologs already in DB)
    if (fromIdx <= SEQ_STEPS.indexOf('step3b') && untilIdx >= SEQ_STEPS.indexOf('step1')) {
        validate_environment()
        download_proteomes(validate_environment.out.validated)
        run_orthofinder(download_proteomes.out.proteomes_dir)
        load_orthologs(run_orthofinder.out.results_dir)
        orthologs_loaded_ch = load_orthologs.out.orthologs_loaded
    } else {
        orthologs_loaded_ch = Channel.value(true)
    }

    // ── Step 3c: nucleotide conservation ─────────────────────────────────
    // Skip if resuming from step4 or later
    if (fromIdx <= SEQ_STEPS.indexOf('step3c') && untilIdx >= SEQ_STEPS.indexOf('step3c')) {
        nucleotide_conservation(orthologs_loaded_ch)
        nuc_done_ch = nucleotide_conservation.out.nuc_done
    } else {
        nuc_done_ch = Channel.value(true)
    }

    // ── Step 4: alignment + divergence scoring ────────────────────────────
    // Skip if resuming from step4b or later (load pre-computed pkl from S3)
    if (fromIdx <= SEQ_STEPS.indexOf('step4') && untilIdx >= SEQ_STEPS.indexOf('step4')) {
        align_and_divergence(orthologs_loaded_ch)
        aligned_pkl_ch = align_and_divergence.out.aligned_pkl
        motifs_done_ch = align_and_divergence.out.motifs_loaded
    } else {
        aligned_pkl_ch = Channel.fromPath("s3://${params.s3_bucket}/cache/aligned_orthogroups.pkl")
        motifs_done_ch = Channel.value(true)
    }

    // ── Steps 4b + 4c: domain annotation + ESM scoring (parallel) ────────
    // Each is independently skippable
    if (fromIdx <= SEQ_STEPS.indexOf('step4b') && untilIdx >= SEQ_STEPS.indexOf('step4b')) {
        domain_and_consequence(motifs_done_ch)
        domains_done_ch = domain_and_consequence.out.domains_done
    } else {
        domains_done_ch = Channel.value(true)
    }

    if (fromIdx <= SEQ_STEPS.indexOf('step4c') && untilIdx >= SEQ_STEPS.indexOf('step4c')) {
        esm1v_scoring(motifs_done_ch)
        esm_done_ch = esm1v_scoring.out.esm_done
    } else {
        esm_done_ch = Channel.value(true)
    }

    // ── Step 4d: variant direction — waits for both 4b and 4c ────────────
    if (untilIdx >= SEQ_STEPS.indexOf('step4d')) {
        variant_direction(domains_done_ch, esm_done_ch)
        directions_done_ch = variant_direction.out.directions_done
    } else {
        directions_done_ch = Channel.value(true)
    }

    emit:
    aligned_pkl     = aligned_pkl_ch
    nuc_done        = nuc_done_ch
    directions_done = directions_done_ch
}
