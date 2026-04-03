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
        --skip-if-done \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step1 || true
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
        --skip-if-done \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --ncbi-api-key '${params.ncbi_api_key ?: ""}'
    # Ensure proteomes dir exists (may be empty symlink if step was skipped)
    STORAGE=\$(python -c "from pipeline.config import get_local_storage_root; print(get_local_storage_root())" 2>/dev/null || echo ".")
    if [ -d "\$STORAGE/proteomes" ]; then
        ln -sf "\$STORAGE/proteomes" proteomes
    else
        mkdir -p proteomes
    fi
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step2 || true
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
        --skip-if-done \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    STORAGE=\$(python -c "from pipeline.config import get_local_storage_root; print(get_local_storage_root())" 2>/dev/null || echo ".")
    if [ -d "\$STORAGE/orthofinder_out" ]; then
        ln -sf "\$STORAGE/orthofinder_out" orthofinder_results
    else
        mkdir -p orthofinder_results
    fi
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step3 || true
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
        --skip-if-done \
        --input-dir '${results_dir}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step3b || true
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
        --skip-if-done \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --ncbi-api-key '${params.ncbi_api_key ?: ""}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step3c || true
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
        --skip-if-done \
        --phenotype '${params.phenotype}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --output-dir '.'
    # If skip-if-done exits 0 (step already complete), the pkl may already be in S3.
    # Ensure the output file exists so Nextflow's output declaration is satisfied.
    if [ ! -f aligned_orthogroups.pkl ]; then
        STORAGE=\$(python -c "from pipeline.config import get_local_storage_root; print(get_local_storage_root())" 2>/dev/null || echo ".")
        if [ -f "\$STORAGE/cache/aligned_orthogroups.pkl" ]; then
            cp "\$STORAGE/cache/aligned_orthogroups.pkl" aligned_orthogroups.pkl
        else
            echo "ERROR: aligned_orthogroups.pkl not found after step4 skip" >&2
            exit 1
        fi
    fi
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step4 || true
    """
}

process domain_and_consequence {
    label 'base'
    cpus 8
    memory '16 GB'
    time '8h'   // 1.9M motifs: UniProt batch + InterPro hybrid fallback needs headroom

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
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step4b || true
    """
}

process esm1v_scatter {
    label 'base'
    cpus 1
    memory '4 GB'
    time '15m'

    input:
    val motifs_loaded

    output:
    path 'chunk_*.txt', emit: chunks, optional: true
    val true, emit: scatter_done

    script:
    """
    export DATABASE_URL='${params.db_url}'
    export BIORESILIENT_STORAGE_ROOT='${params.storage_root}'
    python3 << 'PYEOF'
import math, sys
from db.session import get_session
from db.models import Gene, DivergentMotif, Ortholog
from sqlalchemy import func

with get_session() as s:
    # All genes with at least one motif
    all_genes = set(str(r[0]) for r in (
        s.query(Ortholog.gene_id)
        .join(DivergentMotif, DivergentMotif.ortholog_id == Ortholog.id)
        .group_by(Ortholog.gene_id)
        .having(func.count(DivergentMotif.id) > 0)
        .all()
    ))

    # Genes where EVERY motif already has esm1v_score — fully done
    fully_scored = set()
    for gid in all_genes:
        total = (
            s.query(func.count(DivergentMotif.id))
            .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
            .filter(Ortholog.gene_id == gid)
            .scalar()
        ) or 0
        scored = (
            s.query(func.count(DivergentMotif.id))
            .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
            .filter(Ortholog.gene_id == gid, DivergentMotif.esm1v_score.isnot(None))
            .scalar()
        ) or 0
        if total > 0 and scored == total:
            fully_scored.add(gid)

gene_ids = sorted(all_genes - fully_scored)
sys.stderr.write(
    f'ESM-1v scatter: {len(all_genes)} total genes, '
    f'{len(fully_scored)} fully scored, {len(gene_ids)} remaining\\n'
)

if not gene_ids:
    sys.stderr.write('All genes already ESM-scored — no chunks needed.\\n')
    sys.exit(0)

n_chunks = min(${params.esm_chunks}, max(1, len(gene_ids)))
chunk_size = math.ceil(len(gene_ids) / n_chunks)

for i in range(n_chunks):
    chunk = gene_ids[i * chunk_size : (i + 1) * chunk_size]
    if chunk:
        with open(f'chunk_{i:03d}.txt', 'w') as f:
            f.write('\\n'.join(chunk))
sys.stderr.write(f'Wrote {len(gene_ids)} genes into {n_chunks} chunk files\\n')
PYEOF
    # Create a sentinel chunk if none were written (all scored)
    if ! ls chunk_*.txt 2>/dev/null; then
        echo "# all_scored" > chunk_sentinel.txt
    fi
    """
}

process esm1v_score_chunk {
    label 'esm'
    cpus params.esm_cpus
    memory params.esm_memory
    time '2h'   // 2h: model load ~5 min + ~120 genes × <1 min each = well within limit
    tag "${chunk_file.baseName}"
    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'finish' }
    maxRetries 3

    input:
    path chunk_file

    output:
    val true

    script:
    """
    # Skip the no-op sentinel file created when all genes are already scored
    if grep -q '^# all_scored' '${chunk_file}'; then
        echo "Sentinel chunk — all genes already scored, skipping."
        exit 0
    fi
    python -m scripts.nf_wrappers.run_step \
        --step step4c_chunk \
        --gene-id-file '${chunk_file}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process esm1v_report {
    label 'base'
    cpus 1
    memory '2 GB'
    time '5m'

    input:
    val all_done

    output:
    val true, emit: esm_done

    script:
    """
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step4c || true
    """
}

process variant_direction {
    label 'base'
    cpus 4
    memory '16 GB'  // 1.9M motifs + gnomAD TSV in memory → 16 GB safe
    time '6h'       // DB reads + classifying 1.9M rows + bulk UPDATE

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
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step4d || true
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
        esm1v_scatter(motifs_done_ch)
        esm1v_score_chunk(esm1v_scatter.out.chunks.flatten(), )
        esm1v_report(esm1v_score_chunk.out.collect())
        esm_done_ch = esm1v_report.out.esm_done
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
