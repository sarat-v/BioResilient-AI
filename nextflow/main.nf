#!/usr/bin/env nextflow
/*
 * BioResilient AI v2 — Main Nextflow Workflow
 *
 * A cloud-agnostic pipeline for discovering convergent evolution genes
 * that confer cancer resistance. Wraps the existing Python pipeline with
 * Nextflow orchestration for cloud scalability and per-OG parallelism.
 *
 * Usage:
 *   nextflow run nextflow/main.nf -profile local
 *   nextflow run nextflow/main.nf -profile aws
 *   nextflow run nextflow/main.nf -profile aws --from_step step6
 *   nextflow run nextflow/main.nf -profile gcp --phenotype cancer_resistance
 *
 * Profiles: local, aws, gcp, seqera
 *
 * Key features:
 *   - Per-orthogroup scatter for HyPhy steps (6, 6b, 6c) — spot-safe
 *   - Automatic retries on spot/preemptible interruptions
 *   - Parallel execution of independent steps within each phase
 *   - Cloud-agnostic: same workflow runs on AWS Batch, GCP Batch, or local
 */

nextflow.enable.dsl = 2

include { PHASE1_SEQUENCE   } from './modules/phase1_sequence'
include { PHASE1_EVOLUTION  } from './modules/phase1_evolution'
include { PHASE1_EXPRESSION } from './modules/phase1_expression'
include { PHASE2_CLINICAL   } from './modules/phase2_clinical'

log.info """
╔══════════════════════════════════════════════════════════╗
║           BioResilient AI v2 — Nextflow Pipeline        ║
╠══════════════════════════════════════════════════════════╣
║  Phenotype    : ${params.phenotype}
║  From step    : ${params.from_step}
║  DB URL       : ${params.db_url?.take(40)}...
║  Storage      : ${params.storage_root}
║  Profile      : ${workflow.profile}
║  Nextflow ver : ${nextflow.version}
╚══════════════════════════════════════════════════════════╝
"""

workflow {
    start = Channel.value(true)

    PHASE1_SEQUENCE(start)

    PHASE1_EVOLUTION(
        PHASE1_SEQUENCE.out.aligned_pkl,
        PHASE1_SEQUENCE.out.nuc_done,
        PHASE1_SEQUENCE.out.directions_done
    )

    PHASE1_EXPRESSION(
        PHASE1_EVOLUTION.out.convergent_aa_done,
        PHASE1_EVOLUTION.out.phylo_done
    )

    PHASE2_CLINICAL(
        PHASE1_EXPRESSION.out.scored
    )
}

workflow.onComplete {
    log.info """
    ══════════════════════════════════════════════════════════
    Pipeline complete!
    Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration  : ${workflow.duration}
    Completed : ${workflow.complete}
    Work dir  : ${workflow.workDir}
    ══════════════════════════════════════════════════════════
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline failed: ${workflow.errorMessage}"
}
