#!/bin/bash
# BioResilient Pipeline — Validation Test Runner
#
# Purpose: End-to-end pipeline validation on a fast subset of species.
# This test is designed to give MEANINGFUL results, not just check that
# the code runs. Species were chosen to cover 5 independent lineages across
# two phenotypes (cancer resistance + longevity), which is enough for the
# convergence scoring to produce real Tier1/Tier2 signals on known genes.
#
# Test Panel (8 species):
#   Resilient (5 species, 5 independent lineages):
#     - Naked mole rat    (Rodents)      cancer resistance + longevity
#     - Bowhead whale     (Cetaceans)    longevity + cancer resistance
#     - African elephant  (Proboscideans) cancer resistance (TP53 copies)
#     - Little brown bat  (Bats)         viral tolerance + longevity
#     - Greenland shark   (Sharks)       longevity (400+ yr)
#   Controls (2):
#     - Rat (Rodents control)
#     - Macaque (Primates control)
#   Reference:
#     - Human
#
# Expected Tier1 recoveries: TP53, ATM, ERCC1, FOXO3 pathway members
# Expected runtime: 4-8 hours on itpl-srv3 (10-core Xeon)

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "============================================================"
echo "  BioResilient Pipeline — Validation Test"
echo "============================================================"
echo ""
echo "  Test panel: 5 resilient species × 5 independent lineages"
echo "    · Naked mole rat   (Rodents)"
echo "    · Bowhead whale    (Cetaceans)"
echo "    · African elephant (Proboscideans)"
echo "    · Little brown bat (Bats)"
echo "    · Greenland shark  (Sharks)"
echo "  Controls: Rat, Macaque"
echo "  Reference: Human"
echo ""
  echo "  Expected runtime: 4-8 hrs (10-core server) / 8-14 hrs (8-core Mac M-series)"
echo "  Expected Tier1 genes: TP53, ATM, ERCC1, FOXO3 pathway"
echo "============================================================"
echo ""

# Check conda environment — activate if not already active
if [[ -z "$CONDA_DEFAULT_ENV" ]] || [[ "$CONDA_DEFAULT_ENV" != "bioresillient" ]]; then
    echo ">>> Activating bioresillient conda environment..."
    # Works on both bash (Linux) and zsh (Mac)
    if command -v conda &>/dev/null; then
        eval "$(conda shell.bash hook 2>/dev/null || conda shell.zsh hook 2>/dev/null)"
    elif [ -f "$HOME/miniconda3/bin/conda" ]; then
        eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
    fi
    conda activate bioresillient 2>/dev/null || {
        echo "ERROR: conda environment 'bioresillient' not found."
        echo "Run: bash scripts/setup_mac.sh   (Mac)"
        echo "  or bash scripts/setup_local.sh  (Linux)"
        exit 1
    }
fi

# Backup original species registry
BACKUP="config/species_registry.json.backup"
if [ ! -f "$BACKUP" ]; then
    echo ">>> Backing up config/species_registry.json..."
    cp config/species_registry.json "$BACKUP"
fi

# Activate test registry
echo ">>> Switching to test species registry (8 species)..."
cp config/species_registry_test.json config/species_registry.json

# Wipe any prior test run state so steps don't appear pre-complete
echo ">>> Clearing prior pipeline state and data artifacts..."
rm -f pipeline_state.json pipeline_cache.json pipeline.log pipeline_test.log
# Remove computed data (but keep downloaded proteomes if present to save NCBI bandwidth)
rm -f data/aligned_orthogroups.pkl data/motifs_by_og.pkl 2>/dev/null || true
rm -rf data/alignments/ data/hyphy/ data/phylo/ data/orthofinder_out/ 2>/dev/null || true
echo "    State files cleared."
echo "    NOTE: data/proteomes/ kept intact (re-download skipped if already present)"
echo ""

# Drop and recreate database tables (full clean slate)
echo ">>> Resetting database to clean state..."
python3 - <<'PYEOF'
import sys
sys.path.insert(0, '.')
from db.session import get_engine
from db.models import Base
from sqlalchemy import text

engine = get_engine()
try:
    Base.metadata.drop_all(engine)
    print("    Dropped all tables.")
    Base.metadata.create_all(engine)
    print("    Recreated all tables.")
except Exception as e:
    print(f"    ERROR: {e}")
    sys.exit(1)
PYEOF

# Seed the database with test species
echo ">>> Seeding database with test species..."
python db/seed.py

# Auto-detect available CPU cores (works on both Mac and Linux)
if [[ "$(uname)" == "Darwin" ]]; then
    NCORES=$(sysctl -n hw.physicalcpu)
else
    NCORES=$(nproc)
fi
echo ">>> Detected $NCORES CPU cores"
export ORTHOFINDER_THREADS=$NCORES
export MAFFT_THREADS=$NCORES
export IQTREE_THREADS=$NCORES
export HYPHY_THREADS=$NCORES

echo ""
echo ">>> Starting full pipeline (all 30 steps)..."
echo "    Log: pipeline_test.log"
echo "    Monitor with: tail -f pipeline_test.log"
echo ""

START_TS=$(date +%s)

# Run pipeline — log both to file and stdout so you can see it via SSH
python pipeline/orchestrator.py 2>&1 | tee pipeline_test.log
EXIT_CODE=${PIPESTATUS[0]}

END_TS=$(date +%s)
ELAPSED=$(( (END_TS - START_TS) / 60 ))

# Restore original registry
echo ""
echo ">>> Restoring original species registry..."
mv "$BACKUP" config/species_registry.json

echo ""
echo "============================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "  RESULT: SUCCESS  (${ELAPSED} min)"
    echo "============================================================"
    echo ""
    echo "  Pipeline ran all 30 steps cleanly."
    echo ""
    echo "  Next — run the benchmark:"
    echo "    python scripts/benchmark_recall.py --trait cancer_resistance"
    echo ""
    echo "  Then start the API to inspect results:"
    echo "    uvicorn api.main:app --host 0.0.0.0 --port 8000"
    echo ""
    echo "  Access the dashboard from your Mac:"
    echo "    ssh -L 8000:localhost:8000 <user>@itpl-srv3"
    echo "    Open: http://localhost:8000"
    echo ""
else
    echo "  RESULT: FAILED  (exit code: $EXIT_CODE, ${ELAPSED} min)"
    echo "============================================================"
    echo ""
    echo "  Check log for the failing step:"
    echo "    grep -E 'ERROR|FAILED|Traceback' pipeline_test.log | tail -20"
    echo ""
    echo "  Resume from the last completed step (e.g. step6):"
    echo "    python pipeline/orchestrator.py --resume-from step6b"
    echo ""
fi

exit $EXIT_CODE
