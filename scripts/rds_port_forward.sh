#!/usr/bin/env bash
# Open a local TCP port that forwards to RDS through an EC2 in the same VPC.
#
# RDS is not on the public internet — your Mac must reach it via:
#   (A) This script — SSM Session Manager port forwarding to RDS (recommended), or
#   (B) scripts/sync_rds_to_local.sh — temporarily allow your IP on the RDS SG
#       (only works if the RDS subnet/route is reachable from the internet), or
#   (C) Classic SSH -L 5433:$RDS_HOST:5432 user@bastion-public-ip
#
# Prerequisites (option A):
#   - AWS CLI + session-manager-plugin: https://docs.aws.amazon.com/systems-manager/latest/userguide/session-manager-working-with-install-plugin.html
#   - IAM: ssm:StartSession on the EC2 instance; SSM agent running on the instance
#   - .env: AWS_PROFILE, AWS_REGION, RDS_HOST
#
# Usage:
#   source .env   # or set -a && source .env && set +a
#   export RDS_TUNNEL_EC2_ID=i-xxxxxxxx   # VPC instance that can reach RDS:5432
#   bash scripts/rds_port_forward.sh
#
# Then in a second terminal, point the client at localhost:
#   export DATABASE_URL="postgresql://bioresilient:${RDS_PASSWORD}@127.0.0.1:5433/bioresilient?sslmode=require"
#   psql "$DATABASE_URL" -c 'SELECT 1'
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
if [[ -f "$ROOT/.env" ]]; then
  set -a
  # shellcheck source=/dev/null
  source "$ROOT/.env"
  set +a
fi

AWS_PROFILE="${AWS_PROFILE:-seqera-runner}"
AWS_REGION="${AWS_REGION:-ap-south-1}"
RDS_HOST="${RDS_HOST:?Set RDS_HOST in .env}"
LOCAL_PORT="${RDS_LOCAL_PORT:-5433}"
REMOTE_PORT="${RDS_REMOTE_PORT:-5432}"

# EC2 in the same VPC as RDS; must reach RDS on REMOTE_PORT.
TUNNEL_EC2="${RDS_TUNNEL_EC2_ID:-}"

if [[ -z "$TUNNEL_EC2" ]]; then
  echo "Set RDS_TUNNEL_EC2_ID to a running EC2 instance id (same VPC as RDS), e.g.:"
  echo "  export RDS_TUNNEL_EC2_ID=i-0e9dfdd218ae18990"
  echo ""
  echo "List running instances:"
  echo "  aws ec2 describe-instances --profile $AWS_PROFILE --region $AWS_REGION \\"
  echo "    --filters Name=instance-state-name,Values=running \\"
  echo "    --query 'Reservations[*].Instances[*].[InstanceId,PrivateIpAddress,PublicIpAddress]' --output table"
  exit 1
fi

echo "Forwarding 127.0.0.1:${LOCAL_PORT} -> ${RDS_HOST}:${REMOTE_PORT} via SSM on ${TUNNEL_EC2}"
echo "Leave this terminal open. Use DATABASE_URL with host 127.0.0.1 port ${LOCAL_PORT}"
echo ""

exec aws ssm start-session \
  --profile "$AWS_PROFILE" \
  --region "$AWS_REGION" \
  --target "$TUNNEL_EC2" \
  --document-name AWS-StartPortForwardingSessionToRemoteHost \
  --parameters "{\"host\":[\"$RDS_HOST\"],\"portNumber\":[\"$REMOTE_PORT\"],\"localPortNumber\":[\"$LOCAL_PORT\"]}"
