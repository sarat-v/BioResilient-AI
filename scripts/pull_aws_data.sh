#!/usr/bin/env bash
# BioResilient AI — Pull all AWS data to local Mac
# Downloads: RDS database dump + all S3 files
# Run from the repo root: bash scripts/pull_aws_data.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BACKUP_DIR="$REPO_ROOT/aws-backup"
TIMESTAMP=$(date +"%Y%m%dT%H%M%S")

RDS_HOST="bioresilient.chakgoo8k0sx.ap-southeast-2.rds.amazonaws.com"
RDS_USER="bioresilient"
RDS_DB="bioresilient"
S3_BUCKET="bioresilient-data"
AWS_REGION="ap-southeast-2"

RED="\033[0;31m"; GREEN="\033[0;32m"; YELLOW="\033[1;33m"
CYAN="\033[0;36m"; BOLD="\033[1m"; RESET="\033[0m"

banner() { echo -e "\n${BOLD}${CYAN}━━━  $1  ━━━${RESET}\n"; }
ok()     { echo -e "${GREEN}  ✓ $1${RESET}"; }
warn()   { echo -e "${YELLOW}  ⚠ $1${RESET}"; }
die()    { echo -e "${RED}  ✗ $1${RESET}"; exit 1; }

# ── Check prerequisites ───────────────────────────────────────────────────────
banner "BioResilient — AWS Data Pull"
echo -e "  Backup dir : ${BOLD}${BACKUP_DIR}${RESET}"
echo -e "  Timestamp  : ${TIMESTAMP}"
echo ""

command -v pg_dump  &>/dev/null || die "pg_dump not found. Install: brew install libpq && brew link libpq --force"
command -v aws      &>/dev/null || die "aws CLI not found. Install: brew install awscli"

# ── Prompt for RDS password ───────────────────────────────────────────────────
if [[ -z "${RDS_PASSWORD:-}" ]]; then
    read -r -s -p "  Enter RDS password: " RDS_PASSWORD
    echo ""
fi
export PGPASSWORD="$RDS_PASSWORD"

# ── Verify RDS connectivity ───────────────────────────────────────────────────
banner "Step 1/3: Verifying RDS connectivity"
pg_isready -h "$RDS_HOST" -U "$RDS_USER" -d "$RDS_DB" -q 2>/dev/null \
    || die "Cannot reach RDS at $RDS_HOST. Check VPN/network and that RDS is started."
ok "RDS is reachable"

# Row counts for reference
echo ""
echo -e "  ${BOLD}Current DB row counts:${RESET}"
psql "host=$RDS_HOST dbname=$RDS_DB user=$RDS_USER sslmode=require" \
    -t -c "
SELECT '  ' || table_name || ': ' || to_char(n_live_tup, 'FM999,999,999') || ' rows'
FROM pg_stat_user_tables
ORDER BY n_live_tup DESC;" 2>/dev/null || warn "Could not get row counts (non-fatal)"

# ── Dump RDS ──────────────────────────────────────────────────────────────────
banner "Step 2/3: Dumping RDS database"
mkdir -p "$BACKUP_DIR"
DUMP_FILE="$BACKUP_DIR/bioresilient_${TIMESTAMP}.dump"
DUMP_SQL="$BACKUP_DIR/bioresilient_${TIMESTAMP}.sql"

echo -e "  Dumping to: ${BOLD}${DUMP_FILE}${RESET}"
pg_dump \
    "host=$RDS_HOST dbname=$RDS_DB user=$RDS_USER sslmode=require" \
    --format=custom \
    --file="$DUMP_FILE" \
    --verbose 2>&1 | grep -E "dumping|creating|setting" | head -20 || true

# Also dump plain SQL for easy inspection
pg_dump \
    "host=$RDS_HOST dbname=$RDS_DB user=$RDS_USER sslmode=require" \
    --format=plain \
    --file="$DUMP_SQL" 2>/dev/null

DUMP_SIZE=$(du -sh "$DUMP_FILE" | cut -f1)
SQL_SIZE=$(du -sh "$DUMP_SQL" | cut -f1)
ok "Database dump complete: $DUMP_SIZE (binary) + $SQL_SIZE (SQL)"

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
    2>&1 | grep -v "^$" | tail -5 || warn "S3 sync had warnings (check above)"

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

# Also copy to BioResilient-Results/ for the local results folder
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

# ── Instructions for Dell PC restore ─────────────────────────────────────────
echo -e "${BOLD}${CYAN}━━━  Next: Restore on Dell PC  ━━━${RESET}"
echo ""
echo -e "  ${BOLD}1. Copy backup to Dell PC:${RESET}"
echo -e "     scp -r ${BACKUP_DIR} user@dell-pc-ip:~/BioResilient/"
echo ""
echo -e "  ${BOLD}2. On Dell PC — restore DB:${RESET}"
echo -e "     createdb bioresilient"
echo -e "     pg_restore -d bioresilient -Fc bioresilient_${TIMESTAMP}.dump"
echo ""
  echo -e "  ${BOLD}3. On Dell PC — place S3 files:${RESET}"
  echo -e "     cp -r s3/proteomes/ data/proteomes/"
  echo -e "     mkdir -p /tmp/bioresilient/phylo"
  echo -e "     cp s3/cache/aligned_orthogroups.pkl /tmp/bioresilient/"
  echo -e "     cp s3/cache/species.treefile /tmp/bioresilient/phylo/"
  echo -e "     cp -r step_cache/ ./step_cache/"
  echo -e "     cp -r BioResilient-Results/ ./BioResilient-Results/"
echo ""
echo -e "  ${BOLD}4. On Dell PC — set deployment to local in config/environment.yml:${RESET}"
echo -e "     sed -i 's/deployment: cloud/deployment: local/' config/environment.yml"
echo ""
echo -e "  ${BOLD}5. Run pipeline from step 4c:${RESET}"
echo -e "     ./run_cancer_resistance_stepwise.sh --from step4c"
echo ""
