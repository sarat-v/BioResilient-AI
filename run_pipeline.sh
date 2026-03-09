#!/usr/bin/env bash
# BioResilient AI — Full Production Pipeline Runner
# Usage:
#   bash run_pipeline.sh                  # Run all 12 species (full pipeline)
#   bash run_pipeline.sh --resume         # Resume from last completed step
#   bash run_pipeline.sh --dry-run        # Validate config without running
#
# Prerequisites:
#   1. bash scripts/setup_local.sh        # One-time setup (Linux/WSL2)
#   2. conda activate bioresillient       # Activate environment
#   3. Edit config/environment.yml        # Add NCBI API key and DB settings

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="bioresillient"
LOG_FILE="$REPO_ROOT/pipeline.log"
STATE_FILE="$REPO_ROOT/pipeline_cache.json"

# ── Parse arguments ──────────────────────────────────────────────────────────
RESUME=false
DRY_RUN=false
RESET=false
for arg in "$@"; do
    case $arg in
        --resume)   RESUME=true ;;
        --dry-run)  DRY_RUN=true ;;
        --reset)    RESET=true ;;
        *) echo "Unknown option: $arg"; echo "Usage: $0 [--resume] [--dry-run] [--reset]"; exit 1 ;;
    esac
done

echo "============================================================"
echo " BioResilient AI — Production Pipeline"
echo " Repo:   $REPO_ROOT"
echo " Log:    $LOG_FILE"
echo " Dry:    $DRY_RUN"
echo " Resume: $RESUME"
echo "============================================================"
echo ""

# ── Activate conda ────────────────────────────────────────────────────────────
if command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate "$ENV_NAME" 2>/dev/null || {
        echo "ERROR: conda env '$ENV_NAME' not found. Run: bash scripts/setup_local.sh"
        exit 1
    }
elif [[ -n "${CONDA_PREFIX:-}" ]]; then
    echo "Using active conda env: $CONDA_PREFIX"
else
    echo "WARNING: conda not found — assuming tools are on PATH"
fi

# ── Verify PostgreSQL is running ──────────────────────────────────────────────
if ! pg_isready -q 2>/dev/null; then
    echo "Starting PostgreSQL..."
    # Try systemd first (Ubuntu 22.04+), fall back to service
    sudo systemctl start postgresql 2>/dev/null || sudo service postgresql start 2>/dev/null || {
        echo "ERROR: Could not start PostgreSQL. Start it manually and re-run."
        exit 1
    }
fi

# ── Check config file exists ─────────────────────────────────────────────────
if [[ ! -f "$REPO_ROOT/config/environment.yml" ]]; then
    echo "ERROR: config/environment.yml not found."
    echo "  Run: cp config/environment.example.yml config/environment.yml"
    echo "  Then edit it to add your NCBI API key."
    exit 1
fi

# ── Reset cache if requested ──────────────────────────────────────────────────
if $RESET; then
    echo "Resetting pipeline cache (will re-run all steps)..."
    rm -f "$STATE_FILE"
    rm -f "$REPO_ROOT/pipeline_state.json"
fi

# ── Run pipeline ──────────────────────────────────────────────────────────────
cd "$REPO_ROOT"

PYTHON_CMD="python"
if $DRY_RUN; then
    PYTHON_CMD="python -c \"
import os, sys
sys.path.insert(0, '.')
from pipeline.orchestrator import run_pipeline
from pipeline.config import _load_species_registry
species = _load_species_registry()
print(f'Loaded {len(species)} species from registry:')
for s in species:
    print(f'  {s[\\\"id\\\"]}: {s[\\\"name\\\"]} ({s[\\\"lineage_group\\\"]})') 
print('Dry run — not executing pipeline.')
\""
    eval $PYTHON_CMD
    exit 0
fi

echo "Starting pipeline at $(date)"
echo "Logs streaming to: $LOG_FILE"
echo ""

python - <<'PYEOF'
import os
import sys
import json
from pathlib import Path

sys.path.insert(0, ".")

# ── Load species from registry ────────────────────────────────────────────────
from pipeline.orchestrator import _load_species_registry, run_pipeline

species_list = _load_species_registry()
print(f"Running pipeline with {len(species_list)} species:")
for s in species_list:
    print(f"  {s['id']:25s} {s['name']} ({s['lineage_group']})")
print()

# ── Determine resume point ────────────────────────────────────────────────────
step_order = [
    'step1','step2','step3','step3b','step4','step5','step6',
    'step7','step8','step9','step10','step10b','step11',
    'step12','step13','step14','step15','step16'
]

cache_file = Path('pipeline_cache.json')
completed_steps = []
if cache_file.exists():
    completed_steps = json.loads(cache_file.read_text()).get('completed', [])

# Treat step1 and step2 as done if proteomes exist for ALL species
proteomes_dir = Path('data/proteomes')
proteomes_done = proteomes_dir.exists() and all(
    (proteomes_dir / f"{s['id']}.reheadered.faa").exists()
    for s in species_list
)
if proteomes_done:
    for step in ('step1', 'step2'):
        if step not in completed_steps:
            completed_steps.append(step)

# Find first incomplete step
resume_from = step_order[0]
for step in step_order:
    if step not in completed_steps:
        resume_from = step
        break

if completed_steps:
    print(f"Resuming from: {resume_from}")
    print(f"Already done:  {completed_steps}")
else:
    print("Starting fresh from step1")
print()

run_pipeline(species_list=species_list, resume_from=resume_from)

print()
print("============================================================")
print("Pipeline complete!")
print("To view results:")
print("  uvicorn api.main:app --host 0.0.0.0 --port 8000")
print("============================================================")
PYEOF

echo ""
echo "Pipeline finished at $(date)"
