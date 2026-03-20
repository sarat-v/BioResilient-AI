#!/usr/bin/env bash
#
# BioResilient v2 — Seqera Platform Setup
#
# Configures Seqera Platform (Nextflow Tower) to orchestrate the pipeline.
# Free tier: 250 runs, 3 concurrent, 100 GB Fusion storage.
#
# Prerequisites:
#   - Seqera account at https://cloud.seqera.io
#   - TOWER_ACCESS_TOKEN set (from Seqera → Settings → Access tokens)
#   - AWS Batch already set up (via setup_aws_batch.sh)
#
# Usage:
#   TOWER_ACCESS_TOKEN=your_token bash nextflow/infra/setup_seqera.sh
set -euo pipefail

TOWER_API="https://api.cloud.seqera.io"
TOWER_TOKEN="${TOWER_ACCESS_TOKEN:-}"
REGION="${AWS_REGION:-ap-south-1}"

GREEN="\033[0;32m"; CYAN="\033[0;36m"; BOLD="\033[1m"; RED="\033[0;31m"; RESET="\033[0m"
banner() { echo -e "\n${BOLD}${CYAN}━━━  $1  ━━━${RESET}\n"; }
ok()     { echo -e "${GREEN}  ✓ $1${RESET}"; }
die()    { echo -e "${RED}  ✗ $1${RESET}"; exit 1; }

banner "BioResilient v2 — Seqera Platform Setup"

[[ -n "$TOWER_TOKEN" ]] || die "TOWER_ACCESS_TOKEN not set. Get it from: https://cloud.seqera.io → Settings → Access tokens"

command -v tw &>/dev/null || {
    echo -e "  Installing Seqera CLI (tw) ..."
    curl -fsSL https://get.seqera.io | bash 2>/dev/null
    ok "Installed Seqera CLI"
}

export TOWER_ACCESS_TOKEN="$TOWER_TOKEN"

echo -e "  ${BOLD}Seqera Platform Configuration${RESET}"
echo ""
echo -e "  The Seqera Platform provides:"
echo -e "    - Web UI to monitor pipeline runs"
echo -e "    - Automatic retry of failed tasks"
echo -e "    - Run history and resource usage reports"
echo -e "    - Team collaboration features"
echo ""
echo -e "  ${BOLD}Free tier includes:${RESET}"
echo -e "    - 250 pipeline runs"
echo -e "    - 3 concurrent runs"
echo -e "    - 100 GB Fusion storage/month"
echo ""

banner "Manual Setup Steps"
echo -e "  ${BOLD}1. Create workspace:${RESET}"
echo -e "     Go to https://cloud.seqera.io → Create workspace 'bioresilient'"
echo ""
echo -e "  ${BOLD}2. Add AWS Batch compute environment:${RESET}"
echo -e "     Workspace → Compute Environments → Add"
echo -e "     - Name: aws-batch-spot"
echo -e "     - Platform: AWS Batch"
echo -e "     - Region: ${REGION}"
echo -e "     - Batch queues: bioresilient-spot, bioresilient-large, bioresilient-gpu"
echo -e "     - Work directory: s3://bioresilient-data/nextflow-work/"
echo -e "     - Enable Fusion: Yes"
echo -e "     - Enable Wave: Yes"
echo ""
echo -e "  ${BOLD}3. Add credentials:${RESET}"
echo -e "     Workspace → Credentials → Add"
echo -e "     - AWS Access Key + Secret (for S3 and Batch access)"
echo -e "     - Container Registry (ECR credentials)"
echo ""
echo -e "  ${BOLD}4. Launch pipeline:${RESET}"
echo -e "     Launchpad → Add Pipeline"
echo -e "     - Repository: https://github.com/sarat-v/BioResilient-AI"
echo -e "     - Main script: nextflow/main.nf"
echo -e "     - Config: nextflow/nextflow.config"
echo -e "     - Profile: aws"
echo ""
echo -e "  ${BOLD}5. Or launch from CLI:${RESET}"
echo -e "     export TOWER_ACCESS_TOKEN=${TOWER_TOKEN:0:8}..."
echo -e "     nextflow run nextflow/main.nf -profile aws,seqera"
echo ""

banner "Environment Variables for Nextflow"
echo -e "  Add these to your shell (.bashrc or .zshrc):"
echo ""
echo -e "  export TOWER_ACCESS_TOKEN='${TOWER_TOKEN:0:8}...'"
echo -e "  export NXF_VER=24.04.0"
echo ""
echo -e "  Then run with Seqera monitoring:"
echo -e "  nextflow run nextflow/main.nf -profile aws,seqera"
echo ""
