#!/usr/bin/env bash
#
# BioResilient v2 — AWS Batch Setup for Nextflow
#
# Creates:
#   1. ECR repositories for all Docker images
#   2. IAM roles for Batch execution
#   3. Three Batch compute environments: spot (CPU), large (high-CPU), gpu
#   4. Job queues linking to compute environments
#
# Prerequisites: aws CLI configured, Docker logged into ECR
#
# Usage:
#   bash nextflow/infra/setup_aws_batch.sh
#
# Skip flags:
#   SKIP_ECR=1    — skip ECR repo creation (if they already exist)
#   SKIP_IAM=1    — skip IAM role creation (if they already exist)
set -uo pipefail   # Note: NOT -e — we handle errors explicitly

REGION="${AWS_REGION:-ap-south-1}"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text 2>/dev/null || echo "UNKNOWN")
ECR_REGISTRY="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"
PREFIX="bioresilient"

GREEN="\033[0;32m"; CYAN="\033[0;36m"; YELLOW="\033[1;33m"; BOLD="\033[1m"; RESET="\033[0m"
banner() { echo -e "\n${BOLD}${CYAN}━━━  $1  ━━━${RESET}\n"; }
ok()     { echo -e "${GREEN}  ✓ $1${RESET}"; }
warn()   { echo -e "${YELLOW}  ! $1${RESET}"; }

banner "BioResilient v2 — AWS Batch Setup"
echo -e "  Region   : ${BOLD}${REGION}${RESET}"
echo -e "  Account  : ${BOLD}${ACCOUNT_ID}${RESET}"
echo -e "  ECR      : ${BOLD}${ECR_REGISTRY}${RESET}"
echo ""

# ── 1. Create ECR repositories ──────────────────────────────────────────────
banner "Step 1: Creating ECR repositories"

if [[ "${SKIP_ECR:-0}" == "1" ]]; then
    ok "Skipping ECR (SKIP_ECR=1)"
else
    IMAGES=(base orthofinder align phylo hyphy esm clinical)
    for img in "${IMAGES[@]}"; do
        REPO="${PREFIX}/${img}"
        if aws ecr describe-repositories --repository-names "$REPO" --region "$REGION" &>/dev/null; then
            ok "ECR repo ${REPO} already exists"
        else
            if aws ecr create-repository --repository-name "$REPO" --region "$REGION" --output text &>/dev/null; then
                ok "Created ECR repo ${REPO}"
            else
                warn "Could not create ${REPO} — check ECR permissions"
            fi
        fi
    done
fi

# ── 2. IAM roles ────────────────────────────────────────────────────────────
banner "Step 2: IAM roles"

if [[ "${SKIP_IAM:-0}" == "1" ]]; then
    ok "Skipping IAM (SKIP_IAM=1)"
else
    BATCH_ROLE="${PREFIX}-batch-execution"
    if aws iam get-role --role-name "$BATCH_ROLE" &>/dev/null; then
        ok "IAM role ${BATCH_ROLE} already exists"
    else
        if aws iam create-role \
            --role-name "$BATCH_ROLE" \
            --assume-role-policy-document '{
                "Version": "2012-10-17",
                "Statement": [{
                    "Effect": "Allow",
                    "Principal": {"Service": ["ecs-tasks.amazonaws.com", "batch.amazonaws.com"]},
                    "Action": "sts:AssumeRole"
                }]
            }' --output text &>/dev/null; then

            aws iam attach-role-policy --role-name "$BATCH_ROLE" \
                --policy-arn arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy
            aws iam attach-role-policy --role-name "$BATCH_ROLE" \
                --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess
            aws iam attach-role-policy --role-name "$BATCH_ROLE" \
                --policy-arn arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryReadOnly
            ok "Created IAM role ${BATCH_ROLE}"
        else
            warn "Could not create ${BATCH_ROLE} — check IAM permissions"
        fi
    fi

    SERVICE_ROLE="${PREFIX}-batch-service"
    if aws iam get-role --role-name "$SERVICE_ROLE" &>/dev/null; then
        ok "IAM role ${SERVICE_ROLE} already exists"
    else
        if aws iam create-role \
            --role-name "$SERVICE_ROLE" \
            --assume-role-policy-document '{
                "Version": "2012-10-17",
                "Statement": [{
                    "Effect": "Allow",
                    "Principal": {"Service": "batch.amazonaws.com"},
                    "Action": "sts:AssumeRole"
                }]
            }' --output text &>/dev/null; then

            aws iam attach-role-policy --role-name "$SERVICE_ROLE" \
                --policy-arn arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole
            ok "Created IAM service role ${SERVICE_ROLE}"
        else
            warn "Could not create ${SERVICE_ROLE} — check IAM permissions"
        fi
    fi
fi

# ── 3. Compute Environments ─────────────────────────────────────────────────
banner "Step 3: Batch Compute Environments"

DEFAULT_VPC=$(aws ec2 describe-vpcs --filters Name=isDefault,Values=true \
    --region "$REGION" --query "Vpcs[0].VpcId" --output text 2>/dev/null || echo "")

if [[ -z "$DEFAULT_VPC" || "$DEFAULT_VPC" == "None" ]]; then
    echo "  No default VPC found in ${REGION}."
    echo "  Set SUBNET_IDS and SG_IDS environment variables and re-run."
    exit 1
fi

SUBNET_IDS=$(aws ec2 describe-subnets \
    --filters Name=vpc-id,Values="$DEFAULT_VPC" \
    --region "$REGION" \
    --query "Subnets[*].SubnetId" \
    --output text 2>/dev/null | tr '\t' ',')

SG_IDS=$(aws ec2 describe-security-groups \
    --filters Name=vpc-id,Values="$DEFAULT_VPC" Name=group-name,Values=default \
    --region "$REGION" \
    --query "SecurityGroups[0].GroupId" \
    --output text 2>/dev/null)

echo -e "  VPC     : $DEFAULT_VPC"
echo -e "  Subnets : $SUBNET_IDS"
echo -e "  SG      : $SG_IDS"
echo ""

_create_or_skip_ce() {
    local NAME="$1"
    local JSON="$2"
    local STATUS
    STATUS=$(aws batch describe-compute-environments \
        --compute-environments "$NAME" \
        --region "$REGION" \
        --query "computeEnvironments[0].status" \
        --output text 2>/dev/null || echo "")
    if [[ "$STATUS" == "VALID" || "$STATUS" == "CREATING" ]]; then
        ok "Compute env ${NAME} already exists (${STATUS})"
    else
        if aws batch create-compute-environment \
            --compute-environment-name "$NAME" \
            --type MANAGED \
            --state ENABLED \
            --service-role "arn:aws:iam::${ACCOUNT_ID}:role/${PREFIX}-batch-service" \
            --compute-resources "$JSON" \
            --region "$REGION" --output text &>/dev/null; then
            ok "Created compute environment ${NAME}"
        else
            warn "Could not create ${NAME} — it may already exist or there is a permissions issue"
        fi
    fi
}

SUBNET_JSON=$(echo "$SUBNET_IDS" | tr ',' '\n' | sed 's/.*/"&"/' | paste -sd ',' -)

_create_or_skip_ce "${PREFIX}-spot" "{
    \"type\": \"SPOT\",
    \"allocationStrategy\": \"SPOT_PRICE_CAPACITY_OPTIMIZED\",
    \"minvCpus\": 0, \"maxvCpus\": 512, \"desiredvCpus\": 0,
    \"instanceTypes\": [\"c6i.large\",\"c6i.xlarge\",\"c6i.2xlarge\",\"c5.large\",\"c5.xlarge\",\"c5.2xlarge\",\"m6i.large\",\"m6i.xlarge\",\"m6i.2xlarge\"],
    \"subnets\": [${SUBNET_JSON}],
    \"securityGroupIds\": [\"${SG_IDS}\"],
    \"spotIamFleetRole\": \"arn:aws:iam::${ACCOUNT_ID}:role/aws-ec2-spot-fleet-tagging-role\"
}"

_create_or_skip_ce "${PREFIX}-large" "{
    \"type\": \"SPOT\",
    \"allocationStrategy\": \"SPOT_PRICE_CAPACITY_OPTIMIZED\",
    \"minvCpus\": 0, \"maxvCpus\": 128, \"desiredvCpus\": 0,
    \"instanceTypes\": [\"c6i.4xlarge\",\"c6i.8xlarge\",\"m6i.4xlarge\",\"m6i.8xlarge\",\"r6i.2xlarge\",\"r6i.4xlarge\"],
    \"subnets\": [${SUBNET_JSON}],
    \"securityGroupIds\": [\"${SG_IDS}\"]
}"

_create_or_skip_ce "${PREFIX}-gpu" "{
    \"type\": \"SPOT\",
    \"allocationStrategy\": \"SPOT_PRICE_CAPACITY_OPTIMIZED\",
    \"minvCpus\": 0, \"maxvCpus\": 32, \"desiredvCpus\": 0,
    \"instanceTypes\": [\"g4dn.2xlarge\",\"g4dn.4xlarge\"],
    \"subnets\": [${SUBNET_JSON}],
    \"securityGroupIds\": [\"${SG_IDS}\"]
}"

# ── 4. Job Queues ────────────────────────────────────────────────────────────
banner "Step 4: Job Queues"

for queue_name in spot large gpu; do
    CE="${PREFIX}-${queue_name}"
    QUEUE="${PREFIX}-${queue_name}"
    STATUS=$(aws batch describe-job-queues \
        --job-queues "$QUEUE" \
        --region "$REGION" \
        --query "jobQueues[0].status" \
        --output text 2>/dev/null || echo "")
    if [[ "$STATUS" == "VALID" || "$STATUS" == "CREATING" ]]; then
        ok "Job queue ${QUEUE} already exists (${STATUS})"
    else
        if aws batch create-job-queue \
            --job-queue-name "$QUEUE" \
            --state ENABLED \
            --priority 1 \
            --compute-environment-order "order=1,computeEnvironment=${CE}" \
            --region "$REGION" --output text &>/dev/null; then
            ok "Created job queue ${QUEUE}"
        else
            warn "Could not create job queue ${QUEUE}"
        fi
    fi
done

# ── Summary ──────────────────────────────────────────────────────────────────
banner "Setup Complete"
echo -e "  ${BOLD}ECR Registry:${RESET}  ${ECR_REGISTRY}"
echo -e "  ${BOLD}Compute Envs:${RESET}  ${PREFIX}-spot, ${PREFIX}-large, ${PREFIX}-gpu"
echo -e "  ${BOLD}Job Queues:${RESET}    ${PREFIX}-spot, ${PREFIX}-large, ${PREFIX}-gpu"
echo -e "  ${BOLD}Exec Role ARN:${RESET} arn:aws:iam::${ACCOUNT_ID}:role/${PREFIX}-batch-execution"
echo ""
echo -e "  ${BOLD}Next steps:${RESET}"
echo -e "    1. Build & push Docker images:"
echo -e "       REGISTRY=${ECR_REGISTRY} TAG=latest bash docker/build_all.sh --push"
echo -e "    2. In Seqera Compute Environment, use execution role:"
echo -e "       arn:aws:iam::${ACCOUNT_ID}:role/${PREFIX}-batch-execution"
echo ""
