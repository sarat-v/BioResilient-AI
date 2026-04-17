/*
 * Phase 1 — Evolution Module
 * Steps: 5, 3d, 6 (per-OG scatter), 7, 7b
 *
 * Architecture:
 *   - PAML branch-site model A (step6) is the sole positive-selection test.
 *     It scatters over orthogroups, writing dnds_pvalue (LRT, df=1) and
 *     dnds_ratio (foreground ω) to evolution_score.
 *   - FEL, BUSTED, and RELAX (formerly step6b / step6c) have been removed.
 *     PAML directly answers the biological question (positive selection in
 *     cancer-resistant foreground branches) and is 10-50× faster.
 *   - Convergence scoring (step7) and convergent-AA annotation (step7b)
 *     proceed directly after PAML collect.
 */

process build_species_tree {
    label 'phylo'
    cpus params.iqtree_cpus
    memory params.iqtree_memory
    time '1h'

    input:
    path aligned_pkl

    output:
    path 'species.treefile', emit: treefile

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step5 \
        --skip-if-done \
        --input-pkl '${aligned_pkl}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step5 || true
    """
}

process phylo_conservation {
    label 'align'
    cpus 16
    memory '32 GB'
    time '2h'

    input:
    path treefile
    val  nuc_done

    output:
    val true, emit: phylo_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step3d \
        --skip-if-done \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step3d || true
    """
}

// extract_og_ids queries the DB for the authoritative OG list and writes
// og_batch_NNN.txt files for PAML scatter. storeDir pins the outputs to a
// FIXED, versioned S3 path so that -resume always receives identical input
// file paths → deterministic run_paml task hashes → proper cache hits.
// To force a fresh OG list (e.g. after DB state changes), bump
// params.paml_og_cache_key (e.g. --paml_og_cache_key v2) or delete the
// storeDir prefix from S3.
process extract_og_ids {
    label 'base'
    storeDir "${params.storage_root}/step_cache/${params.phenotype}/paml_og_batches_b${params.paml_batch_size}_${params.paml_og_cache_key}"
    cpus 1
    memory '16 GB'
    time '20m'

    input:
    path aligned_pkl

    output:
    path 'og_batch_*.txt', emit: og_batches

    script:
    """
    python3 << 'PYEOF'
import pickle, sys, math

db_url    = '${params.db_url}'
batch_size = int('${params.paml_batch_size}')

with open('${aligned_pkl}', 'rb') as f:
    data = pickle.load(f)
aligned = data.get('aligned', data)
motifs  = data.get('motifs_by_og', {})

sys.stderr.write(f'pkl: {len(aligned)} aligned OGs, {len(motifs)} in motifs_by_og\\n')

import psycopg2
conn = psycopg2.connect(db_url)
cur  = conn.cursor()

# OGs that have divergent motifs
cur.execute('''
    SELECT DISTINCT o.orthofinder_og
    FROM divergent_motif dm
    JOIN ortholog o ON o.id = dm.ortholog_id
    WHERE o.orthofinder_og IS NOT NULL
''')
db_og_ids = sorted({row[0] for row in cur.fetchall()})

# OGs with a completed PAML or HyPhy run — skip these to avoid redundant recomputation.
# Includes paml_no_signal (PAML ran, found no selection) in addition to paml_branch_site
# (PAML ran, found selection) and busted_ph (legacy HyPhy). Proxy results are intentionally
# excluded so re-runs can upgrade them to a real PAML result.
cur.execute('''
    SELECT DISTINCT o.orthofinder_og
    FROM evolution_score es
    JOIN ortholog o ON o.gene_id = es.gene_id
    WHERE es.selection_model IN (\'paml_branch_site\', \'busted_ph\', \'paml_no_signal\')
      AND o.orthofinder_og IS NOT NULL
''')
already_scored_ogs = {row[0] for row in cur.fetchall()}
conn.close()

sys.stderr.write(f'DB has {len(db_og_ids)} OGs with divergent motifs\\n')
sys.stderr.write(f'DB has {len(already_scored_ogs)} OGs already scored (skip these)\\n')

og_list = [og for og in db_og_ids if og in aligned and og not in already_scored_ogs]
sys.stderr.write(f'OGs remaining to process: {len(og_list)}\\n')

if not og_list:
    sys.stderr.write('All OGs already scored — writing empty sentinel batch.\\n')
    with open('og_batch_0000.txt', 'w') as f:
        f.write('')
    sys.exit(0)

# Pre-filter: skip OGs with < 4 unique species (PAML minimum for reliable LRT)
filtered = []
skipped = 0
for og in og_list:
    seqs = aligned[og]
    species = set()
    for label_seq_pair in seqs:
        label = label_seq_pair[0] if isinstance(label_seq_pair, (list, tuple)) else label_seq_pair
        species.add(label.split('|')[0])
    if len(species) >= 4:
        filtered.append(og)
    else:
        skipped += 1
sys.stderr.write(f'Pre-filter: {len(filtered)} OGs with >=4 species, {skipped} skipped (<4 sp)\\n')
og_list = filtered

if not og_list:
    sys.stderr.write('No OGs pass pre-filter — writing empty sentinel batch.\\n')
    with open('og_batch_0000.txt', 'w') as f:
        f.write('')
    sys.exit(0)

n_batches = math.ceil(len(og_list) / batch_size)
for i in range(n_batches):
    chunk = og_list[i * batch_size : (i + 1) * batch_size]
    with open(f'og_batch_{i:04d}.txt', 'w') as f:
        f.write('\\n'.join(chunk) + '\\n')

sys.stderr.write(f'Wrote {len(og_list)} OGs into {n_batches} batch files (batch_size={batch_size})\\n')
PYEOF
    echo "Batch files: \$(ls og_batch_*.txt | wc -l)"
    """
}

process prefetch_cds {
    label 'base'
    cpus 2
    memory '16 GB'
    time '2h'

    input:
    path aligned_pkl

    output:
    path 'cds_cache.pkl', emit: cds_cache

    script:
    """
    echo "cds_prefetch_v10_shared_cache"
    python -m scripts.nf_wrappers.run_step \
        --step step6_prefetch_cds \
        --input-pkl '${aligned_pkl}' \
        --output-dir . \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --ncbi-api-key '${params.ncbi_api_key ?: ""}'
    ls -lh cds_cache.pkl
    """
}

process run_paml {
    label 'paml'
    cpus params.paml_cpus
    memory { (params.paml_memory as nextflow.util.MemoryUnit) * task.attempt }
    time '3h'
    tag "${og_batch.baseName}"
    maxForks params.paml_max_forks
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    input:
    tuple path(og_batch), path(cds_cache), path(aligned_pkl), path(treefile)

    output:
    path "${og_batch.baseName}", emit: paml_result

    script:
    """
    echo "paml_v2_timeout600s_batch20"
    export NF_TASK_CPUS=${task.cpus}  # actual allocated CPUs so Python runs 1 OG per CPU
    mkdir -p '${og_batch.baseName}'
    python -m scripts.nf_wrappers.run_step \
        --step step6_batch \
        --og-id-file '${og_batch}' \
        --input-pkl '${aligned_pkl}' \
        --treefile '${treefile}' \
        --cds-cache '${cds_cache}' \
        --output-dir '${og_batch.baseName}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --ncbi-api-key '${params.ncbi_api_key ?: ""}'
    """
}

process collect_all_paml_results {
    label 'base'
    cpus 2
    memory '8 GB'

    input:
    path 'results/*'

    output:
    val true, emit: all_paml_done

    script:
    """
    echo "collect_v9_all_genes_per_og"
    python -m scripts.nf_wrappers.run_step \
        --step step6_all_collect \
        --input-dir results/ \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step6 || true
    """
}

process convergence_scoring {
    label 'base'
    cpus 4
    memory '8 GB'
    time '3h'

    input:
    val paml_done

    output:
    val true, emit: convergence_done

    script:
    """
    echo "Running step7 convergence (after PAML collect)"
    python -m scripts.nf_wrappers.run_step \
        --step step7 \
        --skip-if-done \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step7 || true
    """
}

process convergent_aa {
    label 'base'
    cpus 4
    memory '16 GB'
    time '3h'

    input:
    val convergence_done

    output:
    val true, emit: convergent_aa_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step7b \
        --skip-if-done \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step7b || true
    """
}

workflow PHASE1_EVOLUTION {
    take:
    aligned_pkl
    nuc_done
    directions_done

    main:
    def EVO_STEPS = ['step5','step3d','step6','step7','step7b']
    def fromStep   = params.from_step ?: 'step1'
    def untilStep  = params['until'] ?: 'step7b'
    def fromIdx    = EVO_STEPS.indexOf(fromStep)
    def untilIdx   = EVO_STEPS.indexOf(untilStep)
    def step5Idx   = EVO_STEPS.indexOf('step5')
    if (fromIdx < 0) fromIdx = 0
    if (untilIdx < 0) untilIdx = EVO_STEPS.size() - 1

    aligned_pkl_val = aligned_pkl

    // ── Step 5: build species tree ────────────────────────────────────────
    if (fromIdx < 0 || fromIdx <= step5Idx) {
        build_species_tree(aligned_pkl_val)
        treefile_ch  = build_species_tree.out.treefile
    } else {
        treefile_ch  = Channel.fromPath("s3://${params.s3_bucket}/cache/species.treefile").first()
    }

    if (fromIdx <= EVO_STEPS.indexOf('step3d') && untilIdx >= EVO_STEPS.indexOf('step3d')) {
        phylo_conservation(treefile_ch, nuc_done)
        phylo_done_ch = phylo_conservation.out.phylo_done
    } else {
        phylo_done_ch = Channel.value(true)
    }

    if (fromIdx <= EVO_STEPS.indexOf('step6') && untilIdx >= EVO_STEPS.indexOf('step6')) {
        // ── PAML branch-site model A scatter ─────────────────────────────
        extract_og_ids(aligned_pkl_val)
        og_batches_ch = extract_og_ids.out.og_batches.flatten()

        prefetch_cds(aligned_pkl_val)
        cds_cache_ch = prefetch_cds.out.cds_cache

        run_paml_input_ch = og_batches_ch
            .combine(cds_cache_ch)
            .combine(aligned_pkl_val)
            .combine(treefile_ch)

        run_paml(run_paml_input_ch)
        collect_all_paml_results(run_paml.out.paml_result.collect())
        all_paml_done_ch = collect_all_paml_results.out.all_paml_done
    } else {
        all_paml_done_ch = Channel.value(true)
    }

    if (fromIdx <= EVO_STEPS.indexOf('step7') && untilIdx >= EVO_STEPS.indexOf('step7')) {
        convergence_scoring(all_paml_done_ch)
        convergence_done_ch = convergence_scoring.out.convergence_done
    } else {
        convergence_done_ch = Channel.value(true)
    }

    if (fromIdx <= EVO_STEPS.indexOf('step7b') && untilIdx >= EVO_STEPS.indexOf('step7b')) {
        convergent_aa(convergence_done_ch)
        convergent_aa_done_ch = convergent_aa.out.convergent_aa_done
    } else {
        convergent_aa_done_ch = Channel.value(true)
    }

    emit:
    treefile          = treefile_ch
    phylo_done        = phylo_done_ch
    convergent_aa_done = convergent_aa_done_ch
}
