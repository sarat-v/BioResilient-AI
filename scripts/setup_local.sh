#!/usr/bin/env bash
# BioResilient AI — Local environment bootstrap (Ubuntu 22.04 / WSL2)
# Run once on a fresh machine. Idempotent — safe to re-run.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="bioresilient"

echo "============================================"
echo " BioResilient AI — Local Setup"
echo " Repo: $REPO_ROOT"
echo "============================================"

# ── 1. Conda ─────────────────────────────────────────────────────────────────
if ! command -v conda &>/dev/null; then
    echo "[1/6] Installing Miniconda..."
    MINICONDA_INSTALLER="/tmp/miniconda.sh"
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
         -O "$MINICONDA_INSTALLER"
    bash "$MINICONDA_INSTALLER" -b -p "$HOME/miniconda3"
    eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
    conda init bash
    rm "$MINICONDA_INSTALLER"
else
    echo "[1/6] Conda already installed — skipping."
fi

# ── 2. Conda environment ─────────────────────────────────────────────────────
echo "[2/6] Creating/updating conda environment '$ENV_NAME'..."
conda env update --file "$REPO_ROOT/environment.yml" --prune
# Activate for the rest of this script
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

# ── 3. App config ─────────────────────────────────────────────────────────────
CONFIG_FILE="$REPO_ROOT/config/environment.yml"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "[3/6] Creating config/environment.yml from example..."
    cp "$REPO_ROOT/config/environment.example.yml" "$CONFIG_FILE"
    echo "  ➜  Edit $CONFIG_FILE and add your NCBI_API_KEY before running the pipeline."
else
    echo "[3/6] config/environment.yml already exists — skipping."
fi

# ── 4. PostgreSQL ─────────────────────────────────────────────────────────────
echo "[4/6] Setting up PostgreSQL..."
if ! command -v pg_isready &>/dev/null; then
    sudo apt-get update -qq && sudo apt-get install -y postgresql postgresql-contrib
fi

if ! pg_isready -q; then
    sudo service postgresql start
fi

DB_NAME="bioresilient"
if ! psql -U postgres -lqt 2>/dev/null | cut -d '|' -f1 | grep -qw "$DB_NAME"; then
    echo "  Creating database '$DB_NAME'..."
    sudo -u postgres createdb "$DB_NAME"
    sudo -u postgres psql -c "CREATE USER $USER WITH SUPERUSER;" 2>/dev/null || true
else
    echo "  Database '$DB_NAME' already exists — skipping."
fi

# ── 5. Alembic migrations + seed ─────────────────────────────────────────────
echo "[5/6] Running Alembic migrations..."
cd "$REPO_ROOT"
alembic upgrade head

echo "  Seeding species registry..."
python db/seed.py

# ── 6. Tool validation ────────────────────────────────────────────────────────
echo "[6/6] Validating installed tools..."

check_tool() {
    local tool="$1"
    local flag="${2:---version}"
    if command -v "$tool" &>/dev/null; then
        version=$("$tool" $flag 2>&1 | head -1)
        echo "  ✓  $tool — $version"
    else
        echo "  ✗  $tool NOT FOUND — check environment.yml"
    fi
}

check_tool orthofinder
check_tool diamond --version
check_tool mafft --version
check_tool iqtree2 --version
check_tool hyphy --version
check_tool fpocket --version
check_tool python --version

echo ""
python -c "import torch; print('  GPU available:', torch.cuda.is_available(), '|', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU only')"

echo ""
echo "============================================"
echo " Setup complete."
echo " Activate with: conda activate $ENV_NAME"
echo " Run pipeline:  python pipeline/orchestrator.py"
echo "============================================"
