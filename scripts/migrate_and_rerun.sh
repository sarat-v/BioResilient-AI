#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# migrate_and_rerun.sh — Apply accuracy-fix migrations and re-run steps 7→9
#
# Prerequisites
# -------------
# Option A (AWS SSM port-forward):
#   1. aws ssm start-session --profile seqera-runner --region ap-south-1 \
#        --target <INSTANCE_ID> \
#        --document-name AWS-StartPortForwardingSessionToRemoteHost \
#        --parameters '{"portNumber":["5432"],"localPortNumber":["5433"],
#                       "host":["bioresilient.ctso0my4ykwn.ap-south-1.rds.amazonaws.com"]}'
#   2. Update DATABASE_URL in .env: replace RDS_HOST with localhost and 5432 with 5433
#   3. Run this script
#
# Option B (temporarily open security group via AWS Console):
#   1. EC2 → Security Groups → bioresilient-rds-sg → Inbound rules → Add rule:
#      PostgreSQL | TCP | 5432 | My IP
#   2. Run this script
#   3. Delete that inbound rule afterwards
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_DIR"

# Load environment
set -a && source .env && set +a

echo "=== [1/5] Running Alembic migration 0023 (convergence_weight + convergence_pval) ==="
python3 -m alembic upgrade head

echo ""
echo "=== [2/5] Re-running Step 7 (convergence detection + permutation null model) ==="
echo "    This will recompute convergence_weight into the correct column and run"
echo "    200-iteration permutation test to assign convergence_pval."
python3 -c "
from pipeline.orchestrator import step7_convergence
# Force permutation re-run by NULLing existing pvals (migration back-fill preserved old values)
from db.models import EvolutionScore
from db.session import get_session
with get_session() as session:
    session.query(EvolutionScore).update({'convergence_pval': None})
    session.commit()
step7_convergence()
"

echo ""
echo "=== [3/5] Re-running Step 7b (convergent AA with trusted TimeTree topology) ==="
python3 -c "from pipeline.orchestrator import step7b_convergent_aa; step7b_convergent_aa()"

echo ""
echo "=== [4/5] Re-running Step 9 (composite scoring with fixed convergence layers) ==="
python3 -c "from pipeline.orchestrator import step9_scoring; step9_scoring()"

echo ""
echo "=== [5/5] Tier distribution after fixes ==="
python3 -c "
from db.models import CandidateScore
from db.session import get_session
from collections import Counter
with get_session() as session:
    tiers = [cs.tier for cs in session.query(CandidateScore).all()]
counts = Counter(tiers)
for tier in ['Tier1', 'Tier2', 'Tier3']:
    print(f'  {tier}: {counts.get(tier, 0)} genes')
"

echo ""
echo "Done. Review the tier distribution above."
echo "Compare against the original: Tier1 ~120, Tier2 ~350, Tier3 ~11500"
