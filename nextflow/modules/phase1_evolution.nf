/*
 * Phase 1 — Evolution Module
 * Steps: 5, 3d, 6 (per-OG scatter), 6b, 6c, 7, 7b
 *
 * Architecture (2024-redesign):
 *   - PAML (step6) scatters over OGs that need a fresh run.
 *   - FEL/BUSTED (step6b) and RELAX (step6c) are FULLY INDEPENDENT of PAML outputs.
 *     They rebuild codon alignments from aligned_pkl + cds_cache directly, targeting
 *     all paml_branch_site + paml_no_signal OGs that don't yet have fel_sites.
 *   - This eliminates the fragile codon_aln.fna intermediate-file dependency and
 *     makes HyPhy runs robust to Nextflow cache hits on the PAML step.
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

// Backfill busted_pvalue = dnds_pvalue for paml_branch_site genes.
// The new parse_paml_results() emits busted_pvalue, but cached PAML result.json
// files from previous runs don't have this field. This lightweight SQL UPDATE
// ensures consistency without requiring a full PAML rerun.
process backfill_busted_pvalue {
    label 'base'
    cpus 1
    memory '2 GB'
    time '10m'

    input:
    val paml_done

    output:
    val true, emit: backfill_done

    script:
    """
    python3 << 'PYEOF'
import psycopg2
conn = psycopg2.connect('${params.db_url}')
cur = conn.cursor()
cur.execute("""
    UPDATE evolution_score
    SET busted_pvalue = dnds_pvalue
    WHERE selection_model = 'paml_branch_site'
      AND busted_pvalue IS NULL
      AND dnds_pvalue IS NOT NULL
""")
updated = cur.rowcount
conn.commit()
print(f'Backfilled busted_pvalue for {updated} paml_branch_site genes')
cur.execute("""
    SELECT COUNT(*) FROM evolution_score
    WHERE selection_model = 'paml_branch_site' AND busted_pvalue IS NULL
""")
remaining = cur.fetchone()[0]
print(f'Remaining paml_branch_site with NULL busted_pvalue: {remaining}')
conn.close()
PYEOF
    """
}

// extract_og_ids_for_hyphy queries ALL OGs with paml_branch_site or paml_no_signal
// that do not yet have fel_sites filled. Uses storeDir with a versioned key so runs
// are deterministic on -resume. Bump params.hyphy_og_cache_key to force a fresh query
// (e.g. after adding new PAML results to the DB).
process extract_og_ids_for_hyphy {
    label 'base'
    storeDir "${params.storage_root}/step_cache/${params.phenotype}/hyphy_og_batches_b${params.paml_batch_size}_${params.hyphy_og_cache_key}"
    cpus 1
    memory '16 GB'
    time '20m'

    input:
    path aligned_pkl
    val backfill_done

    output:
    path 'hyphy_batch_*.txt', emit: hyphy_batches

    script:
    """
    python3 << 'PYEOF'
import pickle, sys, math

db_url     = '${params.db_url}'
batch_size = int('${params.paml_batch_size}')

with open('${aligned_pkl}', 'rb') as f:
    data = pickle.load(f)
aligned = data.get('aligned', data)

import psycopg2
conn = psycopg2.connect(db_url)
cur = conn.cursor()

# All OGs with a real PAML result that don't have FEL computed yet.
# Include paml_no_signal: their scoring uses the HyPhy pathway, so FEL adds value.
cur.execute('''
    SELECT DISTINCT o.orthofinder_og
    FROM evolution_score es
    JOIN ortholog o ON o.gene_id = es.gene_id
    WHERE es.selection_model IN (\'paml_branch_site\', \'paml_no_signal\')
      AND es.fel_sites IS NULL
      AND o.orthofinder_og IS NOT NULL
''')
og_ids = sorted({row[0] for row in cur.fetchall() if row[0] in aligned})
conn.close()

sys.stderr.write(f'HyPhy targets: {len(og_ids)} OGs with PAML results but no fel_sites\\n')

if not og_ids:
    sys.stderr.write('All paml OGs already have fel_sites — writing empty sentinel.\\n')
    with open('hyphy_batch_0000.txt', 'w') as f:
        f.write('')
    sys.exit(0)

n_batches = math.ceil(len(og_ids) / batch_size)
for i in range(n_batches):
    chunk = og_ids[i * batch_size : (i + 1) * batch_size]
    with open(f'hyphy_batch_{i:04d}.txt', 'w') as f:
        f.write('\\n'.join(chunk) + '\\n')

sys.stderr.write(f'Wrote {len(og_ids)} OGs into {n_batches} HyPhy batch files\\n')
PYEOF
    echo "HyPhy batch files: \$(ls hyphy_batch_*.txt | wc -l)"
    """
}

// run_fel_busted: self-sufficient FEL + BUSTED per batch.
// Receives the same raw inputs as run_paml and rebuilds codon alignments independently.
// No dependency on codon_aln.fna or any other PAML intermediate.
process run_fel_busted {
    label 'paml'
    cpus params.paml_cpus
    memory { (params.paml_memory as nextflow.util.MemoryUnit) * task.attempt }
    time '12h'
    tag "${og_batch.baseName}"
    maxForks params.paml_max_forks
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    input:
    tuple path(og_batch), path(cds_cache), path(aligned_pkl), path(treefile)

    output:
    path "${og_batch.baseName}_fb", emit: fb_result

    script:
    """
    echo "fel_busted_v1_self_sufficient"
    export NF_TASK_CPUS=${task.cpus}
    mkdir -p '${og_batch.baseName}_fb'
    python -m scripts.nf_wrappers.run_step \
        --step step6b_batch \
        --og-id-file '${og_batch}' \
        --input-pkl '${aligned_pkl}' \
        --cds-cache '${cds_cache}' \
        --treefile '${treefile}' \
        --output-dir '${og_batch.baseName}_fb' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

// run_relax: self-sufficient RELAX per batch.
process run_relax {
    label 'paml'
    cpus params.paml_cpus
    memory { (params.paml_memory as nextflow.util.MemoryUnit) * task.attempt }
    time '12h'
    tag "${og_batch.baseName}"
    maxForks params.paml_max_forks
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3

    input:
    tuple path(og_batch), path(cds_cache), path(aligned_pkl), path(treefile)

    output:
    path "${og_batch.baseName}_relax", emit: relax_result

    script:
    """
    echo "relax_v1_self_sufficient"
    export NF_TASK_CPUS=${task.cpus}
    mkdir -p '${og_batch.baseName}_relax'
    python -m scripts.nf_wrappers.run_step \
        --step step6c_batch \
        --og-id-file '${og_batch}' \
        --input-pkl '${aligned_pkl}' \
        --cds-cache '${cds_cache}' \
        --treefile '${treefile}' \
        --output-dir '${og_batch.baseName}_relax' \
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
    echo "collect_v7_self_sufficient"
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
    echo "collect_v7_self_sufficient"
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
    time '3h'

    input:
    val paml_done
    val fb_done
    val relax_done

    output:
    val true, emit: convergence_done

    script:
    """
    echo "Running step7 convergence (after PAML + FEL/BUSTED + RELAX collects)"
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
    def EVO_STEPS = ['step5','step3d','step6','step6b','step6c','step7','step7b']
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
        // ── PAML scatter ─────────────────────────────────────────────────
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

        // ── HyPhy (FEL/BUSTED + RELAX) — self-sufficient scatter ─────────
        // Runs AFTER PAML collect so the DB has the full paml_branch_site set.
        // Targets ALL paml_branch_site + paml_no_signal OGs without fel_sites,
        // not just the ones processed in this PAML scatter — ensuring backfill
        // of existing rows from prior runs.
        if (untilIdx >= EVO_STEPS.indexOf('step6b')) {
            // Backfill busted_pvalue for cached PAML rows that predate the field
            backfill_busted_pvalue(all_paml_done_ch)

            // Extract HyPhy targets from DB (all PAML-scored OGs without FEL)
            extract_og_ids_for_hyphy(aligned_pkl_val, backfill_busted_pvalue.out.backfill_done)
            hyphy_batches_ch = extract_og_ids_for_hyphy.out.hyphy_batches.flatten()

            hyphy_input_ch = hyphy_batches_ch
                .combine(cds_cache_ch)
                .combine(aligned_pkl_val)
                .combine(treefile_ch)

            run_fel_busted(hyphy_input_ch)
            run_relax(hyphy_input_ch)
            collect_fel_busted_results(run_fel_busted.out.fb_result.collect())
            collect_relax_results(run_relax.out.relax_result.collect())
            fb_done_ch    = collect_fel_busted_results.out.fb_done
            relax_done_ch = collect_relax_results.out.relax_done
        } else {
            fb_done_ch    = Channel.value(true)
            relax_done_ch = Channel.value(true)
        }
    } else {
        all_paml_done_ch = Channel.value(true)
        fb_done_ch       = Channel.value(true)
        relax_done_ch    = Channel.value(true)
    }

    if (fromIdx <= EVO_STEPS.indexOf('step7') && untilIdx >= EVO_STEPS.indexOf('step7')) {
        convergence_scoring(all_paml_done_ch, fb_done_ch, relax_done_ch)
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
