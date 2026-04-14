#!/bin/bash
#
# Sync BioResilient database from AWS RDS to local PostgreSQL.
#
# This script:
#   1. Gets your current public IP
#   2. Temporarily adds it to the RDS security group
#   3. Dumps the RDS database
#   4. Removes the security group rule
#   5. Restores the dump to local Postgres
#
# Usage: ./scripts/sync_rds_to_local.sh
#
# Prerequisites:
#   - AWS profile 'seqera-runner' configured with appropriate permissions
#   - psql and pg_dump installed
#   - Local PostgreSQL running on localhost:5432
#   - .env file sourced (AWS_PROFILE, RDS_HOST, RDS_PASSWORD, DATABASE_URL)
#

set -e

# Load environment
if [ -f .env ]; then
    set -a
    source .env
    set +a
else
    echo "❌ .env file not found. Run from project root."
    exit 1
fi

# Configuration
AWS_PROFILE="${AWS_PROFILE:-seqera-runner}"
AWS_REGION="${AWS_REGION:-ap-south-1}"
RDS_HOST="${RDS_HOST:-bioresilient.ctso0my4ykwn.ap-south-1.rds.amazonaws.com}"
RDS_PORT="${RDS_PORT:-5432}"
RDS_USER="${RDS_USER:-bioresilient}"
RDS_DB="${RDS_DB:-bioresilient}"
SECURITY_GROUP_ID="sg-062aa1e65f84c50cf"
LOCAL_USER=$(whoami)
DUMP_FILE="/tmp/bioresilient_rds.dump"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log() {
    echo -e "${GREEN}[$(date +'%H:%M:%S')]${NC} $1"
}

err() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

cleanup_sg_rule() {
    if [ -n "$RULE_DESCRIPTION" ]; then
        log "Removing security group rule..."
        aws ec2 revoke-security-group-ingress \
            --profile "$AWS_PROFILE" --region "$AWS_REGION" \
            --group-id "$SECURITY_GROUP_ID" \
            --protocol tcp --port $RDS_PORT \
            --cidr "$CURRENT_IP/32" \
            2>/dev/null || warn "Could not remove rule (may already be gone)"
    fi
}

trap cleanup_sg_rule EXIT

# ============================================================================
# Step 1: Get current public IP
# ============================================================================
log "Getting current public IP..."
CURRENT_IP=$(curl -s https://checkip.amazonaws.com | tr -d '\n')
if [ -z "$CURRENT_IP" ]; then
    err "Could not determine public IP. Try: curl https://checkip.amazonaws.com"
    exit 1
fi
log "Your IP: $CURRENT_IP"

# ============================================================================
# Step 2: Add IP to RDS security group
# ============================================================================
log "Adding $CURRENT_IP to RDS security group..."
RULE_DESCRIPTION="Temporary access from $CURRENT_IP at $(date +'%Y-%m-%d %H:%M:%S')"

aws ec2 authorize-security-group-ingress \
    --profile "$AWS_PROFILE" --region "$AWS_REGION" \
    --group-id "$SECURITY_GROUP_ID" \
    --protocol tcp --port $RDS_PORT \
    --cidr "$CURRENT_IP/32" \
    --output json >/dev/null 2>&1 || warn "Rule may already exist"

# Wait for rule to be active
log "Waiting for security group rule to activate..."
sleep 3

# ============================================================================
# Step 3: Verify RDS connectivity
# ============================================================================
log "Testing RDS connectivity..."
if ! PGPASSWORD="$RDS_PASSWORD" psql \
    -h "$RDS_HOST" -U "$RDS_USER" -d "$RDS_DB" \
    -c "SELECT 1" &>/dev/null; then
    err "Cannot connect to RDS. Check credentials and security group."
    exit 1
fi
log "✓ RDS connection successful"

# ============================================================================
# Step 4: Dump RDS database
# ============================================================================
log "Dumping RDS database to $DUMP_FILE (this may take a few minutes)..."
PGPASSWORD="$RDS_PASSWORD" pg_dump \
    -h "$RDS_HOST" -U "$RDS_USER" -d "$RDS_DB" \
    --no-owner --no-acl \
    -Fc -f "$DUMP_FILE" \
    2>&1 | grep -v "^pg_dump: warning" || true

if [ ! -f "$DUMP_FILE" ]; then
    err "Dump file not created"
    exit 1
fi

DUMP_SIZE=$(du -h "$DUMP_FILE" | cut -f1)
log "✓ Dump complete ($DUMP_SIZE)"

# ============================================================================
# Step 5: Remove security group rule
# ============================================================================
log "Removing security group rule (cleanup)..."
cleanup_sg_rule
RULE_DESCRIPTION=""  # Mark as cleaned up

# ============================================================================
# Step 6: Restore to local database
# ============================================================================
log "Dropping local 'bioresilient' database..."
dropdb -h localhost -U "$LOCAL_USER" bioresilient 2>/dev/null || true

log "Creating fresh local 'bioresilient' database..."
createdb -h localhost -U "$LOCAL_USER" bioresilient

log "Restoring dump to local database (this may take a few minutes)..."
pg_restore -h localhost -U "$LOCAL_USER" -d bioresilient \
    --no-owner --no-acl --jobs=4 \
    "$DUMP_FILE" 2>&1 | grep -v "^pg_restore: warning" | head -20 || true

log "✓ Restore complete"

# ============================================================================
# Step 7: Verify data
# ============================================================================
log "Verifying data..."
psql -h localhost -U "$LOCAL_USER" bioresilient -c "
SELECT
  (SELECT count(*) FROM gene) as genes,
  (SELECT count(*) FROM candidate_score) as candidates,
  (SELECT count(*) FROM evolution_score) as evolution_scores,
  (SELECT count(*) FROM disease_annotation) as disease_annotations,
  (SELECT count(*) FROM drug_target) as drug_targets,
  (SELECT count(*) FROM species) as species,
  (SELECT count(*) FROM ortholog) as orthologs
;" 2>&1

log "✓ All done! Your local database is now synced with RDS."
log "Run the API and frontend to view the results:"
echo ""
echo "  Backend:   uvicorn api.main:app --host 0.0.0.0 --port 8000"
echo "  Frontend:  cd frontend && npm run dev"
echo ""
