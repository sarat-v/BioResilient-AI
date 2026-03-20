#!/usr/bin/env nextflow
/*
 * BioResilient AI — repository root entry point.
 *
 * Seqera Platform requires a main.nf at the repository root to validate
 * the repository before accepting a custom "Main script" path.
 *
 * The actual pipeline lives at nextflow/main.nf.
 * In the Seqera Launchpad pipeline settings, set:
 *   Main script: nextflow/main.nf
 */
nextflow.enable.dsl = 2
