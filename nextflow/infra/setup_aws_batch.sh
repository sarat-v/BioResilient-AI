#!/usr/bin/env bash
#
# BioResilient v2 — AWS Batch Setup for Nextflow
#
# Creates:
#   1. ECR repositories for all Docker images
#   2. Three Batch compute environments: spot (CPU), large (high-CPU), gpu
#   3. Job queues linking to compute environments
#   4. IAM roles for Batch execution
#
# Prerequisites: aws CLI configured, Docker logged into ECR
#
# Usage:
#   bash nextflow/infra/setup_aws_batch.sh
set -euo pipefail

REGION="${AWS_REGION:-ap-southeast-2}"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text 2>/dev/null || echo "UNKNOWN")
ECR_REGISTRY="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"
PREFIX="bioresilient"

GREEN="\033[0;32m"; CYAN="\033[0;36m"; BOLD="\033[1m"; RESET="\033[0m"
banner() { echo -e "\n${BOLD}${CYAN}━━━  $1  ━━━${RESET}\n"; }
ok()     { echo -e "${GREEN}  ✓ $1${RESET}"; }

banner "BioResilient v2 — AWS Batch Setup"
echo -e "  Region   : ${BOLD}${REGION}${RESET}"
echo -e "  Account  : ${BOLD}${ACCOUNT_ID}${RESET}"
echo -e "  ECR      : ${BOLD}${ECR_REGISTRY}${RESET}"
echo ""

# ── 1. Create ECR repositories ──────────────────────────────────────────────
banner "Step 1: Creating ECR repositories"

IMAGES=(base orthofinder align phylo hyphy esm clinical)
for img in "${IMAGES[@]}"; do
    REPO="${PREFIX}/${img}"
    aws ecr describe-repositories --repository-names "$REPO" --region "$REGION" &>/dev/null \
        && ok "ECR repo ${REPO} already exists" \
        || {
            aws ecr create-repository --repository-name "$REPO" --region "$REGION" --output text &>/dev/null
            ok "Created ECR repo ${REPO}"
        }
done

# ── 2. IAM role for Batch ───────────────────────────────────────────────────
banner "Step 2: IAM roles"

BATCH_ROLE="${PREFIX}-batch-execution"
aws iam get-role --role-name "$BATCH_ROLE" &>/dev/null \
    && ok "IAM role ${BATCH_ROLE} exists" \
    || {
        aws iam create-role \
            --role-name "$BATCH_ROLE" \
            --assume-role-policy-document '{
                "Version": "2012-10-17",
                "Statement": [{
                    "Effect": "Allow",
                    "Principal": {"Service": ["ecs-tasks.amazonaws.com", "batch.amazonaws.com"]},
                    "Action": "sts:AssumeRole"
                }]
            }' --output text &>/dev/null

        aws iam attach-role-policy --role-name "$BATCH_ROLE" \
            --policy-arn arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy
        aws iam attach-role-policy --role-name "$BATCH_ROLE" \
            --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess
        aws iam attach-role-policy --role-name "$BATCH_ROLE" \
            --policy-arn arn:aws:iam::aws:policy/AmazonEC2ContainerRegistryReadOnly

        ok "Created IAM role ${BATCH_ROLE}"
    }

SERVICE_ROLE="${PREFIX}-batch-service"
aws iam get-role --role-name "$SERVICE_ROLE" &>/dev/null \
    && ok "IAM role ${SERVICE_ROLE} exists" \
    || {
        aws iam create-role \
            --role-name "$SERVICE_ROLE" \
            --assume-role-policy-document '{
                "Version": "2012-10-17",
                "Statement": [{
                    "Effect": "Allow",
                    "Principal": {"Service": "batch.amazonaws.com"},
                    "Action": "sts:AssumeRole"
                }]
            }' --output text &>/dev/null

        aws iam attach-role-policy --role-name "$SERVICE_ROLE" \
            --policy-arn arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole

        ok "Created IAM service role ${SERVICE_ROLE}"
    }

# ── 3. Compute Environments ────────────────────────────────────────────────
banner "Step 3: Batch Compute Environments"

# Get default VPC subnets and security groups
DEFAULT_VPC=$(aws ec2 describe-vpcs --filters Name=isDefault,Values=true --region "$REGION" --query "Vpcs[0].VpcId" --output text 2>/dev/null || echo "")
if [[ -z "$DEFAULT_VPC" || "$DEFAULT_VPC" == "None" ]]; then
    echo "  No default VPC found. Provide subnet IDs and security group IDs manually."
    echo "  Set SUBNET_IDS and SG_IDS environment variables and re-run."
    exit 1
fi

SUBNET_IDS=$(aws ec2 describe-subnets --filters Name=vpc-id,Values="$DEFAULT_VPC" --region "$REGION" --query "Subnets[*].SubnetId" --output text 2>/dev/null | tr '\t' ',')
SG_IDS=$(aws ec2 describe-security-groups --filters Name=vpc-id,Values="$DEFAULT_VPC" Name=group-name,Values=default --region "$REGION" --query "SecurityGroups[0].GroupId" --output text 2>/dev/null)

echo -e "  VPC     : $DEFAULT_VPC"
echo -e "  Subnets : $SUBNET_IDS"
echo -e "  SG      : $SG_IDS"
echo ""

# Spot compute environment (for HyPhy per-OG jobs)
aws batch describe-compute-environments --compute-environments "${PREFIX}-spot" --region "$REGION" --query "computeEnvironments[0].status" --output text 2>/dev/null | grep -q VALID \
    && ok "Compute env ${PREFIX}-spot exists" \
    || {
        aws batch create-compute-environment \
            --compute-environment-name "${PREFIX}-spot" \
            --type MANAGED \
            --state ENABLED \
            --service-role "arn:aws:iam::${ACCOUNT_ID}:role/${SERVICE_ROLE}" \
            --compute-resources '{
                "type": "SPOT",
                "allocationStrategy": "SPOT_PRICE_CAPACITY_OPTIMIZED",
                "minvCpus": 0,
                "maxvCpus": 256,
                "desiredvCpus": 0,
                "instanceTypes": ["c6i.xlarge", "c6i.2xlarge", "c5.xlarge", "c5.2xlarge", "m6i.xlarge", "m6i.2xlarge"],
                "subnets": ["'"$(echo $SUBNET_IDS | tr ',' '","')"'"],
                "securityGroupIds": ["'"$SG_IDS"'"],
                "spotIamFleetRole": "arn:aws:iam::'"$ACCOUNT_ID"':role/aws-ec2-spot-fleet-tagging-role"
            }' \
            --region "$REGION" --output text &>/dev/null || true
        ok "Created spot compute environment"
    }

# Large compute environment (for OrthoFinder, alignment)
aws batch describe-compute-environments --compute-environments "${PREFIX}-large" --region "$REGION" --query "computeEnvironments[0].status" --output text 2>/dev/null | grep -q VALID \
    && ok "Compute env ${PREFIX}-large exists" \
    || {
        aws batch create-compute-environment \
            --compute-environment-name "${PREFIX}-large" \
            --type MANAGED \
            --state ENABLED \
            --service-role "arn:aws:iam::${ACCOUNT_ID}:role/${SERVICE_ROLE}" \
            --compute-resources '{
                "type": "SPOT",
                "allocationStrategy": "SPOT_PRICE_CAPACITY_OPTIMIZED",
                "minvCpus": 0,
                "maxvCpus": 64,
                "desiredvCpus": 0,
                "instanceTypes": ["c6i.4xlarge", "c6i.8xlarge", "m6i.4xlarge"],
                "subnets": ["'"$(echo $SUBNET_IDS | tr ',' '","')"'"],
                "securityGroupIds": ["'"$SG_IDS"'"]
            }' \
            --region "$REGION" --output text &>/dev/null || true
        ok "Created large compute environment"
    }

# GPU compute environment (for ESM-1v)
aws batch describe-compute-environments --compute-environments "${PREFIX}-gpu" --region "$REGION" --query "computeEnvironments[0].status" --output text 2>/dev/null | grep -q VALID \
    && ok "Compute env ${PREFIX}-gpu exists" \
    || {
        aws batch create-compute-environment \
            --compute-environment-name "${PREFIX}-gpu" \
            --type MANAGED \
            --state ENABLED \
            --service-role "arn:aws:iam::${ACCOUNT_ID}:role/${SERVICE_ROLE}" \
            --compute-resources '{
                "type": "SPOT",
                "allocationStrategy": "SPOT_PRICE_CAPACITY_OPTIMIZED",
                "minvCpus": 0,
                "maxvCpus": 16,
                "desiredvCpus": 0,
                "instanceTypes": ["g4dn.xlarge", "g4dn.2xlarge"],
                "subnets": ["'"$(echo $SUBNET_IDS | tr ',' '","')"'"],
                "securityGroupIds": ["'"$SG_IDS"'"]
            }' \
            --region "$REGION" --output text &>/dev/null || true
        ok "Created GPU compute environment"
    }

# ── 4. Job Queues ───────────────────────────────────────────────────────────
banner "Step 4: Job Queues"

for queue_name in spot large gpu; do
    CE="${PREFIX}-${queue_name}"
    QUEUE="${PREFIX}-${queue_name}"
    aws batch describe-job-queues --job-queues "$QUEUE" --region "$REGION" --query "jobQueues[0].status" --output text 2>/dev/null | grep -q VALID \
        && ok "Job queue ${QUEUE} exists" \
        || {
            aws batch create-job-queue \
                --job-queue-name "$QUEUE" \
                --state ENABLED \
                --priority 1 \
                --compute-environment-order "order=1,computeEnvironment=${CE}" \
                --region "$REGION" --output text &>/dev/null || true
            ok "Created job queue ${QUEUE}"
        }
done

# ── Summary ─────────────────────────────────────────────────────────────────
banner "Setup Complete"
echo -e "  ${BOLD}ECR Registry:${RESET} ${ECR_REGISTRY}"
echo -e "  ${BOLD}Compute Envs:${RESET} ${PREFIX}-spot, ${PREFIX}-large, ${PREFIX}-gpu"
echo -e "  ${BOLD}Job Queues:${RESET}   ${PREFIX}-spot, ${PREFIX}-large, ${PREFIX}-gpu"
echo ""
echo -e "  ${BOLD}Next steps:${RESET}"
echo -e "    1. Build & push Docker images:"
echo -e "       REGISTRY=${ECR_REGISTRY} TAG=latest bash docker/build_all.sh --push"
echo -e "    2. Run the pipeline:"
echo -e "       nextflow run nextflow/main.nf -profile aws"
echo ""
