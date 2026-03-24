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

// Step order — used to decide which phases to skip when resuming
def STEP_ORDER = [
    'step1','step2','step3','step3b','step3c','step3d',
    'step4','step4b','step4c','step4d',
    'step5','step6','step6b','step6c','step7','step7b',
    'step8','step8b','step9',
]

workflow {
    start = Channel.value(true)

    def fromStep = params.from_step ?: 'step1'
    def fromIdx  = STEP_ORDER.indexOf(fromStep)
    if (fromIdx < 0) {
        error "Unknown from_step value '${fromStep}'. Valid values: ${STEP_ORDER.join(', ')}"
    }

    // ── Phase 1 Sequence (steps 1–4d) ────────────────────────────────────
    // Skip entirely when resuming from step5 or later.
    // aligned_orthogroups.pkl and step3c nucleotide data are loaded from S3.
    if (fromIdx <= STEP_ORDER.indexOf('step4d')) {
        PHASE1_SEQUENCE(start)
        aligned_pkl_ch     = PHASE1_SEQUENCE.out.aligned_pkl
        nuc_done_ch        = PHASE1_SEQUENCE.out.nuc_done
        directions_done_ch = PHASE1_SEQUENCE.out.directions_done
    } else {
        aligned_pkl_ch     = Channel.fromPath("s3://${params.s3_bucket}/cache/aligned_orthogroups.pkl")
        nuc_done_ch        = Channel.value(true)
        directions_done_ch = Channel.value(true)
    }

    // ── Phase 1 Evolution (steps 5–7b) ───────────────────────────────────
    // Skip entirely when resuming from step8 or later.
    // Note: when resuming from step6+, steps 5 and 3d re-run (~20 min each).
    // They are idempotent and fast compared to the HyPhy scatter jobs.
    if (fromIdx <= STEP_ORDER.indexOf('step7b')) {
        PHASE1_EVOLUTION(aligned_pkl_ch, nuc_done_ch, directions_done_ch)
        convergent_aa_done_ch = PHASE1_EVOLUTION.out.convergent_aa_done
        phylo_done_ch         = PHASE1_EVOLUTION.out.phylo_done
    } else {
        convergent_aa_done_ch = Channel.value(true)
        phylo_done_ch         = Channel.value(true)
    }

    // ── Phase 1 Expression + scoring (steps 8–9) ─────────────────────────
    PHASE1_EXPRESSION(convergent_aa_done_ch, phylo_done_ch)

    PHASE2_CLINICAL(PHASE1_EXPRESSION.out.scored)
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
