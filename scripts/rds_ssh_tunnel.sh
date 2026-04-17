#!/usr/bin/env bash
# SSH local port forward: laptop -> EC2 (in VPC) -> RDS:5432
#
# OpenSSH -L forwards so the REMOTE host opens TCP to host:hostport — RDS DNS
# therefore resolves on EC2 (private IP works). Same pattern as scripts/pull_aws_data.sh.
#
# Prerequisites:
#   - EC2 in same VPC as RDS; SG allows EC2 -> RDS:5432
#   - Your IP allowed to EC2:22 (or use a VPN/bastion path you already use)
#   - Key: EC2_KEY path to .pem
#
# Usage (from repo root, after sourcing .env for RDS_HOST / RDS_PASSWORD):
#   export EC2_HOST=43.204.22.192                    # or ec2-...amazonaws.com
#   export EC2_KEY=~/.ssh/your-bioresilient.pem
#   export EC2_USER=ubuntu                         # optional; default ubuntu
#   bash scripts/rds_ssh_tunnel.sh
#
# Second terminal — talk to RDS through the tunnel:
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

RDS_HOST="${RDS_HOST:?Set RDS_HOST (e.g. in .env)}"
LOCAL_PORT="${RDS_LOCAL_PORT:-5433}"
RDS_PORT="${RDS_REMOTE_PORT:-5432}"

EC2_HOST="${EC2_HOST:?Set EC2_HOST — public DNS or IP of VPC EC2}"
EC2_KEY="${EC2_KEY:?Set EC2_KEY — path to .pem for SSH}"
EC2_KEY="${EC2_KEY/#\~/$HOME}"
EC2_USER="${EC2_USER:-ubuntu}"

[[ -f "$EC2_KEY" ]] || { echo "Key not found: $EC2_KEY"; exit 1; }

echo "Tunnel: 127.0.0.1:${LOCAL_PORT} -> ${RDS_HOST}:${RDS_PORT} (via ${EC2_USER}@${EC2_HOST})"
echo "Leave this terminal open. Press Ctrl+C to stop."
echo "Use DATABASE_URL with host 127.0.0.1 port ${LOCAL_PORT}"
echo ""

exec ssh -i "$EC2_KEY" \
  -o StrictHostKeyChecking=accept-new \
  -o ServerAliveInterval=30 \
  -N -L "${LOCAL_PORT}:${RDS_HOST}:${RDS_PORT}" \
  "${EC2_USER}@${EC2_HOST}"
