/*
 * Phase 1 — Evolution Module
 * Steps: 5, 3d, 6 (per-OG scatter), 6b, 6c, 7, 7b
 *
 * Key design: steps 6/6b/6c scatter across orthogroups.
 * Each OG is an independent job — max spot fault tolerance.
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
        --input-pkl '${aligned_pkl}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    STORAGE=\$(python -c "from pipeline.config import get_local_storage_root; print(get_local_storage_root())")
    cp \$STORAGE/phylo/species.treefile species.treefile
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step5 || true
    """
}

process phylo_conservation {
    label 'align'
    cpus 8
    memory '16 GB'
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
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step3d || true
    """
}

// extract_og_ids always queries the DB for authoritative OG list (the pkl's
// motifs_by_og can be stale). The pkl is still needed for aligned sequences.
// Writes og_batch_NNN.txt files with params.hyphy_batch_size OGs each for
// efficient scatter — one pkl load per task instead of one per OG.
process extract_og_ids {
    label 'base'
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
batch_size = int('${params.hyphy_batch_size}')

with open('${aligned_pkl}', 'rb') as f:
    data = pickle.load(f)
aligned = data.get('aligned', data)
motifs  = data.get('motifs_by_og', {})

sys.stderr.write(f'pkl: {len(aligned)} aligned OGs, {len(motifs)} in motifs_by_og\\n')

import psycopg2
conn = psycopg2.connect(db_url)
cur  = conn.cursor()
cur.execute('''
    SELECT DISTINCT o.orthofinder_og
    FROM divergent_motif dm
    JOIN ortholog o ON o.id = dm.ortholog_id
    WHERE o.orthofinder_og IS NOT NULL
''')
db_og_ids = sorted({row[0] for row in cur.fetchall()})
conn.close()

sys.stderr.write(f'DB has {len(db_og_ids)} OGs with divergent motifs\\n')
og_list = [og for og in db_og_ids if og in aligned]
sys.stderr.write(f'OGs in both DB and aligned pkl: {len(og_list)}\\n')

if not og_list:
    sys.stderr.write('WARNING: No OGs found. Falling back to pkl motifs_by_og\\n')
    og_list = sorted(og for og in motifs if og in aligned)
    sys.stderr.write(f'Fallback OGs from pkl: {len(og_list)}\\n')

# Pre-filter: skip OGs with < 4 unique species (HyPhy minimum)
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
sys.stderr.write(f'Pre-filter: {len(filtered)} OGs with >=4 species, {skipped} skipped\\n')
og_list = filtered

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

process run_meme {
    label 'hyphy'
    cpus params.hyphy_cpus
    memory { (params.hyphy_memory as nextflow.util.MemoryUnit) * task.attempt }
    time '3h'
    tag "${og_batch.baseName}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    input:
    path og_batch
    path aligned_pkl
    path treefile
    path cds_cache

    output:
    path "${og_batch.baseName}", emit: meme_result

    script:
    """
    echo "hyphy_v18_spot_bustedph_batch10"
    export HYPHY_CPUS=1          # 1 CPU per tool: 4 tools × 1 = 4 total (BUSTED-PH+FEL+BUSTED+RELAX)
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

process collect_all_hyphy_results {
    label 'base'
    cpus 2
    memory '8 GB'

    input:
    path 'results/*'

    output:
    val true, emit: all_hyphy_done

    script:
    """
    echo "collect_v7_merged_all_hyphy"
    python -m scripts.nf_wrappers.run_step \
        --step step6_all_collect \
        --input-dir results/ \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step6 || true
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step6b || true
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step6c || true
    """
}

process run_fel_busted {
    label 'hyphy'
    cpus params.hyphy_cpus
    memory { (params.hyphy_memory as nextflow.util.MemoryUnit) * task.attempt }
    time '12h'
    tag "${meme_out.baseName}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    input:
    path meme_out
    path treefile

    output:
    path "${meme_out.baseName}_fb", emit: fb_result

    script:
    """
    mkdir -p '${meme_out.baseName}_fb'
    python -m scripts.nf_wrappers.run_step \
        --step step6b_batch \
        --input-dir '${meme_out}' \
        --treefile '${treefile}' \
        --output-dir '${meme_out.baseName}_fb' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process run_relax {
    label 'hyphy'
    cpus params.hyphy_cpus
    memory { (params.hyphy_memory as nextflow.util.MemoryUnit) * task.attempt }
    time '12h'
    tag "${meme_out.baseName}"
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    input:
    path meme_out
    path treefile

    output:
    path "${meme_out.baseName}_relax", emit: relax_result

    script:
    """
    mkdir -p '${meme_out.baseName}_relax'
    python -m scripts.nf_wrappers.run_step \
        --step step6c_batch \
        --input-dir '${meme_out}' \
        --treefile '${treefile}' \
        --output-dir '${meme_out.baseName}_relax' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process collect_fel_busted_results {
    label 'base'
    cpus 2
    memory '8 GB'

    input:
    path 'results/*'

    output:
    val true, emit: fb_done

    script:
    """
    echo "collect_v6_all_fixes"
    python -m scripts.nf_wrappers.run_step \
        --step step6b_collect \
        --input-dir results/ \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step6b || true
    """
}

process collect_relax_results {
    label 'base'
    cpus 2
    memory '8 GB'

    input:
    path 'results/*'

    output:
    val true, emit: relax_done

    script:
    """
    echo "collect_v6_all_fixes"
    python -m scripts.nf_wrappers.run_step \
        --step step6c_collect \
        --input-dir results/ \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    DATABASE_URL='${params.db_url}' BIORESILIENT_STORAGE_ROOT='${params.storage_root}' \
        python -m pipeline.step_reporter --step step6c || true
    """
}

process convergence_scoring {
    label 'base'
    cpus 4
    memory '8 GB'
    time '2h'

    input:
    val meme_done

    output:
    val true, emit: convergence_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step7 \
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
    // Phase execution order — step3d runs after step5 in this workflow
    def EVO_STEPS = ['step5','step3d','step6','step6b','step6c','step7','step7b']
    def fromStep   = params.from_step ?: 'step1'
    def untilStep  = params['until'] ?: 'step7b'
    def fromIdx    = EVO_STEPS.indexOf(fromStep)
    def untilIdx   = EVO_STEPS.indexOf(untilStep)
    def step5Idx   = EVO_STEPS.indexOf('step5')
    if (fromIdx < 0) fromIdx = 0
    if (untilIdx < 0) untilIdx = EVO_STEPS.size() - 1

    // Convert aligned_pkl to a value channel so it can be reused across scatter
    aligned_pkl_val = aligned_pkl.first()

    // ── Step 5: build species tree ────────────────────────────────────────
    // Skip when resuming from step6 or later — tree already in S3.
    if (fromIdx < 0 || fromIdx <= step5Idx) {
        build_species_tree(aligned_pkl_val)
        treefile_ch  = build_species_tree.out.treefile.first()
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
        // Pre-fetch ALL CDS once (single task) — eliminates NCBI bottleneck in scatter
        prefetch_cds(aligned_pkl_val)
        cds_cache_ch = prefetch_cds.out.cds_cache.first()

        extract_og_ids(aligned_pkl_val)
        og_batches_ch = extract_og_ids.out.og_batches.flatten()

        // Step 6: MEME+FEL+BUSTED per OG batch (RELAX deferred)
        run_meme(og_batches_ch, aligned_pkl_val, treefile_ch, cds_cache_ch)
        collect_all_hyphy_results(run_meme.out.meme_result.collect())
        all_hyphy_done_ch = collect_all_hyphy_results.out.all_hyphy_done
    } else {
        all_hyphy_done_ch = Channel.value(true)
    }

    if (fromIdx <= EVO_STEPS.indexOf('step7') && untilIdx >= EVO_STEPS.indexOf('step7')) {
        convergence_scoring(all_hyphy_done_ch)
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
