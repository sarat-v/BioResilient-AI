#!/usr/bin/env bash
# BioResilient AI — Pull all AWS data to local Mac
#
# RDS is in a private VPC — this script dumps the DB remotely on the EC2
# via SSH, uploads the dump to S3, then downloads it here along with all S3 data.
#
# Usage:
#   EC2_HOST=ec2-13-xx-xx.ap-southeast-2.compute.amazonaws.com \
#   EC2_KEY=~/.ssh/my-key.pem \
#   BACKUP_DIR="/Volumes/My Passport/BioResilient-backup" \
#   bash scripts/pull_aws_data.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BACKUP_DIR="${BACKUP_DIR:-$REPO_ROOT/aws-backup}"
TIMESTAMP=$(date +"%Y%m%dT%H%M%S")

RDS_HOST="bioresilient.chakgoo8k0sx.ap-southeast-2.rds.amazonaws.com"
RDS_USER="bioresilient"
RDS_DB="bioresilient"
S3_BUCKET="bioresilient-data"
AWS_REGION="ap-southeast-2"

EC2_HOST="${EC2_HOST:-}"      # public DNS or IP of your EC2
EC2_KEY="${EC2_KEY:-}"        # path to .pem key file
EC2_USER="${EC2_USER:-ubuntu}"

RED="\033[0;31m"; GREEN="\033[0;32m"; YELLOW="\033[1;33m"
CYAN="\033[0;36m"; BOLD="\033[1m"; RESET="\033[0m"

banner() { echo -e "\n${BOLD}${CYAN}━━━  $1  ━━━${RESET}\n"; }
ok()     { echo -e "${GREEN}  ✓ $1${RESET}"; }
warn()   { echo -e "${YELLOW}  ⚠ $1${RESET}"; }
die()    { echo -e "${RED}  ✗ $1${RESET}"; exit 1; }

# ── Check prerequisites ────────────────────────────────────────────────────────
banner "BioResilient — AWS Data Pull"
echo -e "  Backup dir : ${BOLD}${BACKUP_DIR}${RESET}"
echo -e "  Timestamp  : ${TIMESTAMP}"
echo ""

command -v aws &>/dev/null || die "aws CLI not found. Install: brew install awscli"

# ── Collect EC2 connection details ─────────────────────────────────────────────
if [[ -z "$EC2_HOST" ]]; then
    read -r -p "  EC2 public hostname/IP (for SSH): " EC2_HOST
fi
if [[ -z "$EC2_KEY" ]]; then
    read -r -p "  Path to EC2 .pem key file: " EC2_KEY
fi
EC2_KEY="${EC2_KEY/#\~/$HOME}"   # expand ~ if present

[[ -f "$EC2_KEY" ]] || die "Key file not found: $EC2_KEY"

SSH="ssh -i $EC2_KEY -o StrictHostKeyChecking=no -o ConnectTimeout=10 $EC2_USER@$EC2_HOST"

# ── Prompt for RDS password ────────────────────────────────────────────────────
if [[ -z "${RDS_PASSWORD:-}" ]]; then
    read -r -s -p "  Enter RDS password: " RDS_PASSWORD
    echo ""
fi

# ── Verify EC2 SSH connectivity ────────────────────────────────────────────────
banner "Step 1/3: Verifying EC2 SSH connectivity"
$SSH "echo ok" &>/dev/null || die "Cannot SSH to $EC2_HOST. Check that EC2 is running and key is correct."
ok "SSH to EC2 works"

# ── Dump RDS via EC2 and upload to S3 ─────────────────────────────────────────
banner "Step 2/3: Dumping RDS database via EC2"
REMOTE_DUMP="/tmp/bioresilient_${TIMESTAMP}.dump"
REMOTE_SQL="/tmp/bioresilient_${TIMESTAMP}.sql"
S3_DUMP_KEY="db-backups/bioresilient_${TIMESTAMP}.dump"
S3_SQL_KEY="db-backups/bioresilient_${TIMESTAMP}.sql"

echo -e "  Dumping on EC2 → ${REMOTE_DUMP} …"
# Pass password via environment using SSH SendEnv / inline export to avoid quoting issues
$SSH "bash -s" <<REMOTESCRIPT
export PGPASSWORD=$(printf '%q' "$RDS_PASSWORD")
# Use pg_dump 15 to match RDS server version (15)
PG_DUMP=\$(ls /usr/lib/postgresql/15/bin/pg_dump /usr/bin/pg_dump 2>/dev/null | head -1)
\$PG_DUMP "host=${RDS_HOST} dbname=${RDS_DB} user=${RDS_USER} sslmode=require" \
    --format=custom --file="${REMOTE_DUMP}" 2>&1 || { echo "DUMP_FAILED"; exit 1; }
\$PG_DUMP "host=${RDS_HOST} dbname=${RDS_DB} user=${RDS_USER} sslmode=require" \
    --format=plain  --file="${REMOTE_SQL}"  2>/dev/null || true
echo "dump_ok"
REMOTESCRIPT
[[ $? -eq 0 ]] || die "pg_dump on EC2 failed. Check RDS password and connectivity from EC2."

ok "DB dumped on EC2"

echo -e "  Uploading dump to s3://${S3_BUCKET}/${S3_DUMP_KEY} …"
$SSH "bash -s" <<REMOTESCRIPT2
aws s3 cp "${REMOTE_DUMP}" "s3://${S3_BUCKET}/${S3_DUMP_KEY}" --region ${AWS_REGION} --no-progress
aws s3 cp "${REMOTE_SQL}"  "s3://${S3_BUCKET}/${S3_SQL_KEY}"  --region ${AWS_REGION} --no-progress
echo "upload_ok"
REMOTESCRIPT2
[[ $? -eq 0 ]] || die "S3 upload from EC2 failed."

ok "Dump uploaded to S3"

# ── Download dump from S3 to local ────────────────────────────────────────────
mkdir -p "$BACKUP_DIR"
DUMP_FILE="$BACKUP_DIR/bioresilient_${TIMESTAMP}.dump"
DUMP_SQL="$BACKUP_DIR/bioresilient_${TIMESTAMP}.sql"

echo -e "  Downloading dump to local: ${DUMP_FILE} …"
aws s3 cp "s3://${S3_BUCKET}/${S3_DUMP_KEY}" "$DUMP_FILE" --region "$AWS_REGION" --no-progress
aws s3 cp "s3://${S3_BUCKET}/${S3_SQL_KEY}"  "$DUMP_SQL"  --region "$AWS_REGION" --no-progress

DUMP_SIZE=$(du -sh "$DUMP_FILE" | cut -f1)
SQL_SIZE=$(du -sh "$DUMP_SQL" | cut -f1)
ok "Database dump downloaded: $DUMP_SIZE (binary) + $SQL_SIZE (SQL)"

# Clean up temp files on EC2
$SSH "rm -f '${REMOTE_DUMP}' '${REMOTE_SQL}'" 2>/dev/null || true

# ── Sync S3 ───────────────────────────────────────────────────────────────────
banner "Step 3/3: Syncing S3 bucket"
S3_LOCAL="$BACKUP_DIR/s3"
mkdir -p "$S3_LOCAL"

echo -e "  Syncing s3://${S3_BUCKET}/ → ${S3_LOCAL}/"
echo -e "  ${YELLOW}(This may take a while for large proteome/genome files)${RESET}"
echo ""

aws s3 sync "s3://${S3_BUCKET}/" "$S3_LOCAL/" \
    --region "$AWS_REGION" \
    --no-progress \
    --exclude "*.tmp" \
    --exclude "db-backups/*" \
    2>&1 | grep -v "^$" | tail -10 || warn "S3 sync had warnings (check above)"

# Copy the dump files into s3/ as well so everything is in one place
cp "$DUMP_FILE" "$S3_LOCAL/db-backups/" 2>/dev/null || { mkdir -p "$S3_LOCAL/db-backups" && cp "$DUMP_FILE" "$S3_LOCAL/db-backups/"; }
cp "$DUMP_SQL"  "$S3_LOCAL/db-backups/" 2>/dev/null || true

S3_SIZE=$(du -sh "$S3_LOCAL" 2>/dev/null | cut -f1 || echo "unknown")
ok "S3 sync complete: $S3_SIZE downloaded"

# ── Pull step cache files into local step_cache/ dir ─────────────────────────
banner "Bonus: Pulling step cache reports"
STEP_CACHE_DIR="$REPO_ROOT/step_cache"
mkdir -p "$STEP_CACHE_DIR"

echo -e "  Syncing step reports → ${STEP_CACHE_DIR}/"
aws s3 sync "s3://${S3_BUCKET}/step_cache/cancer_resistance/" "$STEP_CACHE_DIR/" \
    --region "$AWS_REGION" \
    --no-progress 2>/dev/null \
    && ok "Step cache reports downloaded to step_cache/" \
    || warn "Step cache sync had warnings (non-fatal)"

RESULTS_DIR="$REPO_ROOT/BioResilient-Results"
mkdir -p "$RESULTS_DIR"
for f in "$STEP_CACHE_DIR"/*.json "$STEP_CACHE_DIR"/*.md; do
    [[ -f "$f" ]] && cp "$f" "$RESULTS_DIR/" 2>/dev/null || true
done
ok "Results also copied to BioResilient-Results/"

# ── Summary ───────────────────────────────────────────────────────────────────
banner "Download Complete"
echo -e "  ${BOLD}Files saved to: ${BACKUP_DIR}${RESET}"
echo ""
echo -e "  DB dump (binary): ${BOLD}${DUMP_FILE}${RESET}"
echo -e "  DB dump (SQL)   : ${BOLD}${DUMP_SQL}${RESET}"
echo -e "  S3 files        : ${BOLD}${S3_LOCAL}/${RESET}"
echo -e "  Step reports    : ${BOLD}${STEP_CACHE_DIR}/${RESET}"
echo -e "  Results folder  : ${BOLD}${RESULTS_DIR}/${RESET}"
echo ""
echo -e "  ${BOLD}Total backup size:${RESET}"
du -sh "$BACKUP_DIR"
echo ""

# ── Instructions for Dell PC restore ──────────────────────────────────────────
echo -e "${BOLD}${CYAN}━━━  Next: Restore on Dell PC  ━━━${RESET}"
echo ""
echo -e "  ${BOLD}1. Copy backup to Dell PC (from My Passport):${RESET}"
echo -e "     rsync -avz '/Volumes/My Passport/BioResilient-backup/' user@dell-pc-ip:~/BioResilient-backup/"
echo ""
echo -e "  ${BOLD}2. On Dell PC — restore DB:${RESET}"
echo -e "     createdb bioresilient"
echo -e "     pg_restore -d bioresilient -Fc ~/BioResilient-backup/bioresilient_${TIMESTAMP}.dump"
echo ""
echo -e "  ${BOLD}3. On Dell PC — place S3 files:${RESET}"
echo -e "     cp -r ~/BioResilient-backup/s3/proteomes/ ~/BioResilient-AI/data/proteomes/"
echo -e "     mkdir -p /tmp/bioresilient/phylo"
echo -e "     cp ~/BioResilient-backup/s3/cache/aligned_orthogroups.pkl /tmp/bioresilient/"
echo -e "     cp ~/BioResilient-backup/s3/cache/species.treefile /tmp/bioresilient/phylo/"
echo -e "     cp -r ~/BioResilient-backup/s3/step_cache/ ~/BioResilient-AI/step_cache/"
echo ""
echo -e "  ${BOLD}4. Set deployment to local in config/environment.yml:${RESET}"
echo -e "     sed -i 's/deployment: cloud/deployment: local/' config/environment.yml"
echo ""
echo -e "  ${BOLD}5. Run pipeline from step 6:${RESET}"
echo -e "     ./run_cancer_resistance_stepwise.sh --from step6"
echo ""
