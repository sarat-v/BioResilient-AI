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
 *   nextflow run nextflow/main.nf -profile aws --until step4d
 *   nextflow run nextflow/main.nf -profile gcp --phenotype cancer_resistance
 *
 * Profiles: local, aws, gcp, seqera
 *
 * Key features:
 *   - Per-orthogroup scatter for PAML step 6 — spot-safe
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
║  Until step   : ${params.until ?: 'step18'}
║  DB URL       : ${params.db_url?.take(40)}...
║  Storage      : ${params.storage_root}
║  Profile      : ${workflow.profile}
║  Nextflow ver : ${nextflow.version}
╚══════════════════════════════════════════════════════════╝
"""

// Actual workflow execution order — used for from_step / until gating.
// Pipeline boundary: step18 (target dossier) is the final computational step.
// Steps 19-25 (virtual screening, lead optimisation) are a separate product
// requiring wet-lab inputs and are NOT part of this pipeline.
def STEP_ORDER = [
    'step1','step2','step3','step3b','step3c',
    'step4','step4b','step4c','step4d',
    'step5','step3d','step6','step7','step7b',
    'step8','step8b','step9','step9b',
    'step10b','step11','step11b','step11c','step11d',
    'step12','step12b','step13','step14','step14b','step15',
    'step16','step17','step18',
]

// ── Phase 3 processes (steps 16-18) ────────────────────────────────────────
// These are lightweight API / DB-write steps — 2 CPUs, 4 GB RAM each.

process translational_status {
    label 'base'
    cpus 2
    memory '4 GB'
    time '1h'

    input:
    val ready

    output:
    val true, emit: trans_done

    script:
    """
    export DATABASE_URL='${params.db_url}'
    export BIORESILIENT_STORAGE_ROOT='${params.storage_root}'
    python -m scripts.nf_wrappers.run_step \
        --step step16 \
        --phenotype '${params.phenotype}'
    """
}

process preclinical_readiness {
    label 'base'
    cpus 2
    memory '4 GB'
    time '1h'

    input:
    val ready

    output:
    val true, emit: preclinical_done

    script:
    """
    export DATABASE_URL='${params.db_url}'
    export BIORESILIENT_STORAGE_ROOT='${params.storage_root}'
    python -m scripts.nf_wrappers.run_step \
        --step step17 \
        --phenotype '${params.phenotype}'
    """
}

process target_dossier {
    label 'base'
    cpus 2
    memory '4 GB'
    time '1h'

    input:
    val ready

    output:
    val true

    script:
    """
    export DATABASE_URL='${params.db_url}'
    export BIORESILIENT_STORAGE_ROOT='${params.storage_root}'
    python -m scripts.nf_wrappers.run_step \
        --step step18 \
        --phenotype '${params.phenotype}'
    """
}

workflow {
    start = Channel.value(true)

    def fromStep  = params.from_step ?: 'step1'
    def fromIdx   = STEP_ORDER.indexOf(fromStep)
    def untilStep = params.until ?: STEP_ORDER[-1]
    def untilIdx  = STEP_ORDER.indexOf(untilStep)
    if (fromIdx < 0) {
        error "Unknown from_step value '${fromStep}'. Valid values: ${STEP_ORDER.join(', ')}"
    }
    if (untilIdx < 0) {
        error "Unknown until value '${untilStep}'. Valid values: ${STEP_ORDER.join(', ')}"
    }
    if (fromIdx > untilIdx) {
        error "from_step '${fromStep}' must come before or equal to until '${untilStep}'"
    }

    def overlaps = { String firstStep, String lastStep ->
        fromIdx <= STEP_ORDER.indexOf(lastStep) && untilIdx >= STEP_ORDER.indexOf(firstStep)
    }

    // ── Phase 1 Sequence (steps 1–4d) ────────────────────────────────────
    // Skip entirely when resuming from step5 or later.
    // aligned_orthogroups.pkl and step3c nucleotide data are loaded from S3.
    if (overlaps('step1', 'step4d')) {
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
    if (overlaps('step5', 'step7b')) {
        PHASE1_EVOLUTION(aligned_pkl_ch, nuc_done_ch, directions_done_ch)
        convergent_aa_done_ch = PHASE1_EVOLUTION.out.convergent_aa_done
        phylo_done_ch         = PHASE1_EVOLUTION.out.phylo_done
    } else {
        convergent_aa_done_ch = Channel.value(true)
        phylo_done_ch         = Channel.value(true)
    }

    // ── Phase 1 Expression + scoring (steps 8–9) ─────────────────────────
    if (overlaps('step8', 'step9')) {
        PHASE1_EXPRESSION(convergent_aa_done_ch, phylo_done_ch)
        scored_ch = PHASE1_EXPRESSION.out.scored
    } else {
        scored_ch = Channel.value(true)
    }

    // ── Phase 2 Clinical translation (steps 9b–15) ───────────────────────
    if (overlaps('step9b', 'step15')) {
        PHASE2_CLINICAL(scored_ch)
        phase2_done_ch = PHASE2_CLINICAL.out.pipeline_done
    } else {
        phase2_done_ch = Channel.value(true)
    }

    // ── Phase 3 Translational / Dossier (steps 16–18) ────────────────────
    // Pipeline boundary: step18 (target dossier) is the final computational
    // output. Steps 19-25 (virtual screening, ADMET, lead optimisation)
    // require wet-lab validation inputs and are a separate downstream product.
    if (overlaps('step16', 'step16')) {
        translational_status(phase2_done_ch)
        translational_done_ch = translational_status.out.trans_done
    } else {
        translational_done_ch = Channel.value(true)
    }

    if (overlaps('step17', 'step17')) {
        preclinical_readiness(translational_done_ch)
        preclinical_done_ch = preclinical_readiness.out.preclinical_done
    } else {
        preclinical_done_ch = Channel.value(true)
    }

    if (overlaps('step18', 'step18')) {
        target_dossier(preclinical_done_ch)
    }
}

// ── Standalone health check ────────────────────────────────────────────────
// Run with: nextflow run main.nf -entry health_check -profile aws
workflow health_check {
    process run_health_check {
        label 'base'
        cpus 1
        memory '4 GB'
        time '30m'

        output:
        path 'health_check.md'
        path 'health_check.json'

        script:
        """
        export DATABASE_URL='${params.db_url}'
        export BIORESILIENT_STORAGE_ROOT='${params.storage_root}'

        echo "=== Per-step reports ==="
        STEPS="step1 step2 step3 step3b step3c step4 step4b step4c step4d step3d step5 step6 step7 step7b step8 step8b step9"
        for step in \$STEPS; do
            echo "--- \$step ---"
            python -m pipeline.step_reporter --step "\$step" || true
        done

        echo ""
        echo "=== Pipeline Health Check ==="
        python /app/scripts/pipeline_health_check.py --output health_check.md || true
        """
    }
}

// ── Data quality check for SKIPPED steps ──────────────────────────────────
// Runs inside AWS Batch (VPC → RDS access). Reports PASS/FAIL per skipped step.
// Run with: nextflow run nextflow/main.nf -entry data_quality_check -profile aws

process run_data_quality_check {
    label 'base'
    cpus 1
    memory '4 GB'
    time '30m'

    script:
    """
    export DATABASE_URL='${params.db_url}'
    export AWS_DEFAULT_REGION='ap-south-1'
    export S3_BUCKET='${params.s3_bucket}'

    # Download the quality check script from S3 using Python/boto3 (no aws CLI needed)
    python3 - << 'PYEOF'
import boto3, os, subprocess, sys

s3 = boto3.client('s3', region_name='ap-south-1')
bucket = os.environ.get('S3_BUCKET', 'bioresilient-data')

# Download the script
s3.download_file(bucket, 'scripts/data_quality_check.py', '/tmp/data_quality_check.py')

# Run it
result = subprocess.run([sys.executable, '/tmp/data_quality_check.py'], text=True)

# Upload results to S3
import os as _os
if _os.path.exists('data_quality_report.json'):
    s3.upload_file('data_quality_report.json', bucket, 'reports/data_quality_report.json')
    print(f"Report uploaded to s3://{bucket}/reports/data_quality_report.json")

sys.exit(result.returncode)
PYEOF
    """
}

workflow data_quality_check {
    run_data_quality_check()
}

// ── Stage 1 (step4c–4d) success verification ──────────────────────────────
// Runs inside AWS Batch (VPC → RDS). Local machines cannot reach private RDS.
// Run with: nextflow run nextflow/main.nf -entry stage1_verify -profile aws

process stage1_success_verify {
    label 'base'
    cpus 1
    memory '2 GB'
    time '15m'

    script:
    """
    export DATABASE_URL='${params.db_url}'
    python3 << 'PYEOF'
import os
import sys
import psycopg2

conn = psycopg2.connect(os.environ['DATABASE_URL'])
cur = conn.cursor()

print('=== Stage 1 success checks (after step4c + step4d) ===')
# Note: exact 0.0 can be a genuine mean LLR (neutral substitution), not only
# the old sentinel. Pre-step A clears false zeros before 4c; 4c may write 0.0 again.
cur.execute('SELECT COUNT(*) FROM divergent_motif WHERE esm1v_score = 0.0')
zeros = cur.fetchone()[0]
cur.execute('SELECT COUNT(*) FROM divergent_motif WHERE esm1v_score IS NOT NULL')
scored = cur.fetchone()[0]
pct_zero = 100.0 * zeros / scored if scored else 0.0
print(f"1) esm1v_score = 0.0: {zeros:,} of {scored:,} scored ({pct_zero:.1f}%) — OK if neutral LLRs; investigate if >40%")

cur.execute('SELECT COUNT(*) FROM divergent_motif WHERE esm1v_score IS NULL')
nulls = cur.fetchone()[0]
cur.execute('SELECT COUNT(*) FROM divergent_motif')
total = cur.fetchone()[0]
print(f"2) esm1v_score IS NULL: {nulls} / {total} motifs")

cur.execute(
    'SELECT motif_direction, COUNT(*) AS n,'
    ' ROUND(100.0 * COUNT(*) / SUM(COUNT(*)) OVER (), 2) AS pct'
    ' FROM divergent_motif WHERE motif_direction IS NOT NULL'
    ' GROUP BY motif_direction ORDER BY n DESC'
)
rows = cur.fetchall()
print('')
print('3) motif_direction distribution:')
for r in rows:
    print(f"   {str(r[0]):22s}  n={r[1]:>10,}  pct={r[2]}%")

neutral_pct = None
for r in rows:
    if r[0] == 'neutral':
        neutral_pct = float(r[2])
        break
if neutral_pct is not None:
    ok = neutral_pct < 90.0
    print('')
    print(f"4) neutral < 90% (plan gate): {neutral_pct}%  {'PASS' if ok else 'FAIL (variant direction still dominated by neutral)'}")

cur.execute(
    "SELECT COUNT(*) FROM divergent_motif WHERE motif_direction = 'loss_of_function'"
)
lof = cur.fetchone()[0]
print(f"5) loss_of_function count: {lof:,}  {'PASS' if lof > 0 else 'WARN (zero LoF — check LOEUF / step4d logic)'}")

conn.close()
print('')
print('=== End Stage 1 checks ===')
# Exit 0: informational verify; inspect printed FAIL/WARN above for biology.
sys.exit(0)
PYEOF
    """
}

workflow stage1_verify {
    stage1_success_verify()
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
