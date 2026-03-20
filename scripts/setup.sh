#!/usr/bin/env bash
# BioResilient AI — Full AWS setup in one command (ap-south-1 / Mumbai)
#
# What this does:
#   1. Creates S3 bucket
#   2. Creates RDS PostgreSQL 15 instance
#   3. Restores S3 data from your local backup
#   4. Creates ECR repos + Batch compute environments + job queues
#   5. Builds and pushes all 7 Docker images
#   6. Writes a .env file with all the env vars you need
#
# Usage:
#   RDS_PASSWORD="<your-password>" \
#   BACKUP_DIR="/Volumes/My Passport/BioResilient-backup" \
#   bash scripts/setup.sh
#
# Optional overrides:
#   REGION        (default: ap-south-1)
#   S3_BUCKET     (default: bioresilient-data)
#   DB_INSTANCE   (default: bioresilient)
#   SKIP_DOCKER=1 to skip image build/push (if already pushed)
#   SKIP_S3=1     to skip S3 restore (if already done)
set -euo pipefail

# ── Config ────────────────────────────────────────────────────────────────────
REGION="${REGION:-ap-south-1}"
S3_BUCKET="${S3_BUCKET:-bioresilient-data}"
DB_INSTANCE="${DB_INSTANCE:-bioresilient}"
DB_USER="bioresilient"
DB_NAME="bioresilient"
RDS_PASSWORD="${RDS_PASSWORD:?  ERROR: set RDS_PASSWORD env var before running}"
BACKUP_DIR="${BACKUP_DIR:-}"
SKIP_DOCKER="${SKIP_DOCKER:-0}"
SKIP_S3="${SKIP_S3:-0}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
ECR="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"

GREEN="\033[0;32m"; CYAN="\033[0;36m"; YELLOW="\033[1;33m"; BOLD="\033[1m"; RESET="\033[0m"
banner() { echo -e "\n${BOLD}${CYAN}━━━  $1  ━━━${RESET}\n"; }
ok()     { echo -e "${GREEN}  ✓ $1${RESET}"; }
warn()   { echo -e "${YELLOW}  ! $1${RESET}"; }

banner "BioResilient AI — AWS Setup (${REGION})"
echo -e "  Account  : ${BOLD}${ACCOUNT_ID}${RESET}"
echo -e "  Region   : ${BOLD}${REGION}${RESET}"
echo -e "  S3 bucket: ${BOLD}${S3_BUCKET}${RESET}"
echo -e "  ECR      : ${BOLD}${ECR}${RESET}"
echo ""

# ── 1. S3 bucket ──────────────────────────────────────────────────────────────
banner "Step 1/5: S3 bucket"
if aws s3api head-bucket --bucket "$S3_BUCKET" --region "$REGION" 2>/dev/null; then
    ok "Bucket s3://${S3_BUCKET} already exists"
else
    aws s3 mb "s3://${S3_BUCKET}" --region "$REGION"
    ok "Created s3://${S3_BUCKET}"
fi

if [[ "$SKIP_S3" == "1" ]]; then
    warn "Skipping S3 restore (SKIP_S3=1)"
elif [[ -n "$BACKUP_DIR" && -d "$BACKUP_DIR/s3" ]]; then
    echo "  Syncing ${BACKUP_DIR}/s3/ → s3://${S3_BUCKET}/ ..."
    aws s3 sync "${BACKUP_DIR}/s3/" "s3://${S3_BUCKET}/" --region "$REGION"
    ok "S3 data restored"
else
    warn "BACKUP_DIR not set or missing — skipping S3 restore. Run manually:"
    warn "  aws s3 sync '/Volumes/My Passport/BioResilient-backup/s3/' s3://${S3_BUCKET}/ --region ${REGION}"
fi

# ── 2. RDS ────────────────────────────────────────────────────────────────────
banner "Step 2/5: RDS PostgreSQL 15"
if aws rds describe-db-instances --db-instance-identifier "$DB_INSTANCE" \
       --region "$REGION" --query "DBInstances[0].DBInstanceStatus" \
       --output text 2>/dev/null | grep -qE "available|creating|backing-up"; then
    ok "RDS instance ${DB_INSTANCE} already exists"
else
    echo "  Creating RDS instance (takes ~5 min, running in background)..."
    DEFAULT_VPC=$(aws ec2 describe-vpcs \
        --filters Name=isDefault,Values=true \
        --region "$REGION" --query "Vpcs[0].VpcId" --output text)
    DEFAULT_SG=$(aws ec2 describe-security-groups \
        --filters "Name=vpc-id,Values=${DEFAULT_VPC}" "Name=group-name,Values=default" \
        --region "$REGION" --query "SecurityGroups[0].GroupId" --output text)

    aws rds create-db-instance \
        --db-instance-identifier "$DB_INSTANCE" \
        --db-instance-class db.t3.medium \
        --engine postgres \
        --engine-version "15" \
        --master-username "$DB_USER" \
        --master-user-password "$RDS_PASSWORD" \
        --db-name "$DB_NAME" \
        --allocated-storage 50 \
        --storage-type gp3 \
        --no-publicly-accessible \
        --vpc-security-group-ids "$DEFAULT_SG" \
        --region "$REGION" \
        --output text > /dev/null
    ok "RDS creation started (takes ~5 min in background)"
    echo "  Check status: aws rds describe-db-instances --db-instance-identifier ${DB_INSTANCE} --region ${REGION} --query 'DBInstances[0].DBInstanceStatus'"
fi

RDS_HOST=$(aws rds describe-db-instances \
    --db-instance-identifier "$DB_INSTANCE" \
    --region "$REGION" \
    --query "DBInstances[0].Endpoint.Address" \
    --output text 2>/dev/null || echo "")
if [[ -n "$RDS_HOST" && "$RDS_HOST" != "None" ]]; then
    ok "RDS endpoint: ${RDS_HOST}"
else
    warn "RDS endpoint not ready yet — will be in the .env file once available"
    RDS_HOST="<pending — check aws rds describe-db-instances>"
fi

# ── 3. AWS Batch infrastructure ───────────────────────────────────────────────
banner "Step 3/5: AWS Batch infrastructure"
bash "${REPO_ROOT}/nextflow/infra/setup_aws_batch.sh"

# ── 4. ECR login + Docker build + push ───────────────────────────────────────
banner "Step 4/5: Docker images → ECR"
aws ecr get-login-password --region "$REGION" | \
    docker login --username AWS --password-stdin "$ECR"
ok "ECR login successful"

if [[ "$SKIP_DOCKER" == "1" ]]; then
    warn "Skipping Docker build/push (SKIP_DOCKER=1)"
else
    echo "  Building and pushing 7 images (this takes ~20 min)..."
    REGISTRY="${ECR}/bioresilient" TAG="latest" \
        bash "${REPO_ROOT}/docker/build_all.sh" --push
fi

# ── 5. Write .env file ────────────────────────────────────────────────────────
banner "Step 5/5: Writing .env"
ENV_FILE="${REPO_ROOT}/.env"
cat > "$ENV_FILE" <<EOF
# BioResilient AI — environment variables
# Source this before starting the API: source .env
# Auto-generated by scripts/setup.sh on $(date)

# AWS
export AWS_REGION="${REGION}"
export S3_BUCKET="${S3_BUCKET}"

# RDS
export RDS_HOST="${RDS_HOST}"
export RDS_PASSWORD="${RDS_PASSWORD}"

# Nextflow / Seqera
export NXF_WORK="s3://${S3_BUCKET}/nextflow-work"
# Fill these in after setting up Seqera Platform (cloud.seqera.io):
export SEQERA_API_TOKEN=""
export SEQERA_WORKSPACE_ID=""
export SEQERA_PIPELINE_ID=""
export SEQERA_ORG_NAME=""
export SEQERA_WORKSPACE_NAME=""
EOF
ok "Wrote ${ENV_FILE}"

# ── Done ──────────────────────────────────────────────────────────────────────
banner "Setup Complete"
echo -e "  ${BOLD}Two manual steps remaining:${RESET}"
echo ""
echo -e "  ${BOLD}A) Restore DB from backup${RESET} (skip if fresh run)"
echo "     Wait for RDS to be 'available', then on a temp EC2 in the same VPC:"
echo "     pg_restore -h ${RDS_HOST} -U ${DB_USER} -d ${DB_NAME} -Fc backup.dump"
echo ""
echo -e "  ${BOLD}B) Set up Seqera Platform${RESET}"
echo "     1. cloud.seqera.io → Credentials → Add AWS"
echo "     2. Compute Environments → AWS Batch → region ${REGION}"
echo "        queues: bioresilient-spot, bioresilient-large, bioresilient-gpu"
echo "     3. Launchpad → Add pipeline → github.com/sarat-v/BioResilient-AI"
echo "        main script: nextflow/main.nf   profile: aws"
echo "     4. Fill SEQERA_* values in .env (printed above)"
echo ""
echo -e "  ${BOLD}Then start the API:${RESET}"
echo "     source .env && uvicorn api.main:app --host 0.0.0.0 --port 8000"
echo ""
echo -e "  ${BOLD}Or run pipeline directly:${RESET}"
echo "     source .env && nextflow run nextflow/main.nf -profile aws --from_step step6 -resume"
echo ""
