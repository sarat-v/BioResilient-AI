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
    time '4h'

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
    """
}

process extract_og_ids {
    label 'base'
    cpus 1
    memory '8 GB'

    input:
    path aligned_pkl

    output:
    path 'og_ids.txt', emit: og_list

    script:
    """
    python3 -c "
import pickle, sys
with open('${aligned_pkl}', 'rb') as f:
    data = pickle.load(f)
aligned = data.get('aligned', data)
motifs = data.get('motifs_by_og', {})
# Only scatter OGs that have motifs (candidates for HyPhy)
for og_id in sorted(motifs.keys()):
    if og_id in aligned:
        print(og_id)
" > og_ids.txt
    echo "OGs to scatter: \$(wc -l < og_ids.txt)"
    """
}

process run_meme {
    label 'hyphy'
    cpus params.hyphy_cpus
    memory params.hyphy_memory
    time '30m'
    tag "${og_id}"

    input:
    val og_id
    path aligned_pkl
    path treefile

    output:
    tuple val(og_id), path("${og_id}/"), emit: meme_result

    script:
    """
    mkdir -p '${og_id}'
    python -m scripts.nf_wrappers.run_step \
        --step step6_single_og \
        --og-id '${og_id}' \
        --input-pkl '${aligned_pkl}' \
        --treefile '${treefile}' \
        --output-dir '${og_id}' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}' \
        --ncbi-api-key '${params.ncbi_api_key ?: ""}'
    """
}

process collect_meme_results {
    label 'base'
    cpus 2
    memory '8 GB'

    input:
    path 'results/*'

    output:
    val true, emit: meme_done

    script:
    """
    python -m scripts.nf_wrappers.run_step \
        --step step6_collect \
        --input-dir results/ \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process run_fel_busted {
    label 'hyphy'
    cpus params.hyphy_cpus
    memory params.hyphy_memory
    time '30m'
    tag "${og_id}"

    input:
    tuple val(og_id), path(meme_dir)
    path treefile

    output:
    tuple val(og_id), path("${og_id}_fb/"), emit: fb_result

    script:
    """
    mkdir -p '${og_id}_fb'
    CODON_ALN='${meme_dir}/codon_aln.fna'
    python -m scripts.nf_wrappers.run_step \
        --step step6b_single_og \
        --og-id '${og_id}' \
        --codon-aln "\$CODON_ALN" \
        --treefile '${treefile}' \
        --output-dir '${og_id}_fb' \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process run_relax {
    label 'hyphy'
    cpus params.hyphy_cpus
    memory params.hyphy_memory
    time '30m'
    tag "${og_id}"

    input:
    tuple val(og_id), path(meme_dir)
    path treefile

    output:
    tuple val(og_id), path("${og_id}_relax/"), emit: relax_result

    script:
    """
    mkdir -p '${og_id}_relax'
    CODON_ALN='${meme_dir}/codon_aln.fna'
    python -m scripts.nf_wrappers.run_step \
        --step step6c_single_og \
        --og-id '${og_id}' \
        --codon-aln "\$CODON_ALN" \
        --treefile '${treefile}' \
        --output-dir '${og_id}_relax' \
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
    python -m scripts.nf_wrappers.run_step \
        --step step6b_collect \
        --input-dir results/ \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
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
    python -m scripts.nf_wrappers.run_step \
        --step step6c_collect \
        --input-dir results/ \
        --db-url '${params.db_url}' \
        --storage-root '${params.storage_root}'
    """
}

process convergence_scoring {
    label 'base'
    cpus 4
    memory '8 GB'
    time '1h'

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
    """
}

process convergent_aa {
    label 'base'
    cpus 4
    memory '8 GB'

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
    """
}

workflow PHASE1_EVOLUTION {
    take:
    aligned_pkl
    nuc_done
    directions_done

    main:
    // Step 5: build species tree (needs aligned_pkl)
    build_species_tree(aligned_pkl)

    // Step 3d: phylo conservation (needs tree + nucleotide data)
    phylo_conservation(build_species_tree.out.treefile, nuc_done)

    // Extract OG IDs for scatter
    extract_og_ids(aligned_pkl)
    og_ids_ch = extract_og_ids.out.og_list
        .splitText()
        .map { it.trim() }
        .filter { it.length() > 0 }

    // Step 6: MEME — per-OG scatter
    run_meme(og_ids_ch, aligned_pkl, build_species_tree.out.treefile)

    // Collect MEME results, write to DB
    collect_meme_results(run_meme.out.meme_result.map { it[1] }.collect())

    // Step 6b: FEL+BUSTED — reuses codon alignments from MEME
    run_fel_busted(run_meme.out.meme_result, build_species_tree.out.treefile)
    collect_fel_busted_results(run_fel_busted.out.fb_result.map { it[1] }.collect())

    // Step 6c: RELAX — also reuses codon alignments
    run_relax(run_meme.out.meme_result, build_species_tree.out.treefile)
    collect_relax_results(run_relax.out.relax_result.map { it[1] }.collect())

    // Step 7: convergence (after all selection results in DB)
    convergence_scoring(
        collect_meme_results.out.meme_done
            .combine(collect_fel_busted_results.out.fb_done)
            .combine(collect_relax_results.out.relax_done)
            .map { true }
    )

    // Step 7b: convergent AA
    convergent_aa(convergence_scoring.out.convergence_done)

    emit:
    treefile          = build_species_tree.out.treefile
    phylo_done        = phylo_conservation.out.phylo_done
    convergent_aa_done = convergent_aa.out.convergent_aa_done
}
