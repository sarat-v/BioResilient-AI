#!/bin/bash
# BioResilient Pipeline — Validation Test Runner
#
# Two modes:
#
#   --quick   (default) Runs only the key validation gates through step9
#             using the stepwise runner in non-interactive mode.
#             Exits 1 if any FAIL gate is detected.
#             Expected runtime: 4-8h (Phase 1 only).
#
#   --full    Runs all 21 step groups end-to-end, non-interactive.
#             Expected runtime: 8-14h.
#
# Test Panel (8 species, from species_registry_test.json):
#   Resilient (5):  naked_mole_rat, bowhead_whale, african_elephant,
#                   little_brown_bat, greenland_shark
#   Controls (2):   rat, macaque
#   Reference (1):  human
#
# Expected Tier1 recoveries: TP53, ATM, ERCC1

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

MODE="quick"
if [[ "${1:-}" == "--full" ]]; then
    MODE="full"
fi

echo "============================================================"
echo "  BioResilient Pipeline — Validation Test ($MODE)"
echo "============================================================"
echo ""
echo "  Test panel: 5 resilient × 5 independent lineages + 3 controls"
echo "  Mode: $MODE (non-interactive, exit 1 on FAIL)"
echo "  Log:  pipeline_test.log"
echo "============================================================"
echo ""

# ── Conda activation ──────────────────────────────────────────────────────────
if [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "$CONDA_DEFAULT_ENV" != "bioresilient" ]]; then
    echo ">>> Activating bioresilient conda environment..."
    if command -v conda &>/dev/null; then
        eval "$(conda shell.bash hook 2>/dev/null || conda shell.zsh hook 2>/dev/null)"
    elif [ -f "$HOME/miniconda3/bin/conda" ]; then
        eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
    fi
    conda activate bioresilient 2>/dev/null || {
        echo "ERROR: conda environment 'bioresilient' not found."
        echo "Run: bash scripts/setup_cloud.sh"
        exit 1
    }
fi

# ── Switch to test species registry ───────────────────────────────────────────
BACKUP="config/species_registry.json.validation_backup"
echo ">>> Switching to test species registry (8 species)..."
cp config/species_registry.json "$BACKUP"
cp config/species_registry_test.json config/species_registry.json

restore_registry() {
    if [[ -f "$BACKUP" ]]; then
        mv "$BACKUP" config/species_registry.json
        echo ">>> Original species registry restored."
    fi
}
trap restore_registry EXIT

# ── Clear prior state ─────────────────────────────────────────────────────────
echo ">>> Clearing prior pipeline state..."
rm -f pipeline_state.json pipeline_cache.json pipeline.log pipeline_test.log
rm -f data/aligned_orthogroups.pkl data/motifs_by_og.pkl 2>/dev/null || true
rm -rf data/alignments/ data/hyphy/ data/phylo/ data/orthofinder_out/ 2>/dev/null || true
rm -rf step_cache/ 2>/dev/null || true
echo "    Cleared. Proteomes in data/proteomes/ are preserved."
echo ""

# ── Reset database ────────────────────────────────────────────────────────────
echo ">>> Resetting database..."
python3 - <<'PYEOF'
import sys
sys.path.insert(0, '.')
from db.session import get_engine
from db.models import Base

engine = get_engine()
try:
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)
    print("    Database reset complete.")
except Exception as e:
    print(f"    ERROR: {e}")
    sys.exit(1)
PYEOF

# Seed test species
echo ">>> Seeding database with test species..."
python db/seed.py

echo ""
echo ">>> Starting stepwise runner (non-interactive, mode=$MODE)..."
echo ""

START_TS=$(date +%s)

if [[ "$MODE" == "quick" ]]; then
    # Quick validation: only steps needed for Phase 1 checkpoint
    bash run_cancer_resistance_stepwise.sh \
        --non-interactive \
        --only "step1,step2,step3,step3b,step4,step4b,step4c,step4d,step5,step6,step6b,step6c,step7,step7b,step8,step8b,step9" \
        2>&1 | tee pipeline_test.log
    EXIT_CODE=${PIPESTATUS[0]}
else
    # Full run: all step groups
    bash run_cancer_resistance_stepwise.sh \
        --non-interactive \
        2>&1 | tee pipeline_test.log
    EXIT_CODE=${PIPESTATUS[0]}
fi

END_TS=$(date +%s)
ELAPSED=$(( (END_TS - START_TS) / 60 ))

echo ""
echo "============================================================"
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "  RESULT: SUCCESS  (${ELAPSED} min)"
    echo "============================================================"
    echo ""
    echo "  All validation gates passed."
    echo ""
    echo "  Run benchmark recall:"
    echo "    python scripts/benchmark_recall.py --trait cancer_resistance"
    echo ""
    echo "  Inspect cached reports:"
    echo "    ls step_cache/"
    echo "    cat step_cache/step9.md"
    echo ""
else
    echo "  RESULT: FAILED  (exit code: $EXIT_CODE, ${ELAPSED} min)"
    echo "============================================================"
    echo ""
    echo "  One or more validation gates failed. Investigate:"
    echo "    grep -E 'FAIL|ERROR|Traceback' pipeline_test.log | tail -20"
    echo ""
    echo "  Step reports are in step_cache/ — check the last completed step."
    echo ""
fi

exit $EXIT_CODE
