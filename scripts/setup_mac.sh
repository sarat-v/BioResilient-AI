#!/usr/bin/env bash
# BioResilient AI — Mac Setup Script (Apple Silicon + Intel)
# Run once on a fresh Mac. Idempotent — safe to re-run.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="bioresillient"

echo "============================================"
echo " BioResilient AI — Mac Setup"
echo " Repo: $REPO_ROOT"
echo " Arch: $(uname -m)"
echo "============================================"
echo ""

# ── 1. Homebrew ───────────────────────────────────────────────────────────────
echo "[1/7] Checking Homebrew..."
if ! command -v brew &>/dev/null; then
    echo "  Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    # Add brew to PATH for Apple Silicon
    if [[ "$(uname -m)" == "arm64" ]]; then
        eval "$(/opt/homebrew/bin/brew shellenv)"
        echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
    fi
else
    echo "  Homebrew already installed — skipping."
fi

# ── 2. Miniconda (Apple Silicon native) ──────────────────────────────────────
echo "[2/7] Checking Conda..."
if ! command -v conda &>/dev/null; then
    echo "  Installing Miniconda for $(uname -m)..."
    if [[ "$(uname -m)" == "arm64" ]]; then
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
    else
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    fi
    curl -fsSL "$MINICONDA_URL" -o /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
    eval "$("$HOME/miniconda3/bin/conda" shell.zsh hook)"
    conda init zsh
    rm /tmp/miniconda.sh
    echo "  Miniconda installed at ~/miniconda3"
else
    echo "  Conda already installed — skipping."
    eval "$(conda shell.zsh hook 2>/dev/null || conda shell.bash hook)"
fi

# ── 3. Conda environment ──────────────────────────────────────────────────────
echo "[3/7] Creating/updating conda environment '$ENV_NAME'..."
# Note: fpocket is not on conda-forge for macOS — installed via brew below
conda env update --file "$REPO_ROOT/environment.yml" --prune -n "$ENV_NAME" || \
    conda env create --file "$REPO_ROOT/environment.yml"

eval "$(conda shell.zsh hook 2>/dev/null || conda shell.bash hook)"
conda activate "$ENV_NAME"
echo "  Environment '$ENV_NAME' active."

# ── 4. Extra tools via Homebrew (not in conda-forge for Apple Silicon) ────────
echo "[4/7] Installing bioinformatics tools via Homebrew..."

brew_install() {
    local pkg="$1"
    local cmd="${2:-$1}"
    if command -v "$cmd" &>/dev/null; then
        echo "  ✓  $cmd already installed"
    else
        echo "  Installing $pkg..."
        brew install "$pkg"
    fi
}

brew_install mafft mafft
# fpocket — build from source (no bottle for arm64)
if ! command -v fpocket &>/dev/null; then
    echo "  Building fpocket from source..."
    FPOCKET_TMP=$(mktemp -d)
    git clone --depth=1 https://github.com/Discngine/fpocket.git "$FPOCKET_TMP/fpocket"
    make -C "$FPOCKET_TMP/fpocket" -j4
    sudo make -C "$FPOCKET_TMP/fpocket" install
    rm -rf "$FPOCKET_TMP"
    echo "  ✓  fpocket installed"
else
    echo "  ✓  fpocket already installed"
fi

# ── 5. PostgreSQL ─────────────────────────────────────────────────────────────
echo "[5/7] Setting up PostgreSQL..."
if ! command -v psql &>/dev/null; then
    echo "  Installing PostgreSQL via Homebrew..."
    brew install postgresql@16
    brew link postgresql@16 --force
fi

# Start PostgreSQL service
if ! pg_isready -q 2>/dev/null; then
    echo "  Starting PostgreSQL..."
    brew services start postgresql@16 2>/dev/null || brew services start postgresql 2>/dev/null
    sleep 3
fi

DB_NAME="bioresillient"
if psql -d "$DB_NAME" -c "" 2>/dev/null; then
    echo "  Database '$DB_NAME' already exists — skipping."
else
    echo "  Creating database '$DB_NAME'..."
    createdb "$DB_NAME" 2>/dev/null || true
fi

# ── 6. App config ─────────────────────────────────────────────────────────────
echo "[6/7] Checking app config..."
CONFIG_FILE="$REPO_ROOT/config/environment.yml"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "  Creating config/environment.yml from example..."
    cp "$REPO_ROOT/config/environment.example.yml" "$CONFIG_FILE"
    echo ""
    echo "  *** ACTION REQUIRED ***"
    echo "  Edit $CONFIG_FILE"
    echo "  Add your NCBI API key (free at https://www.ncbi.nlm.nih.gov/account/)"
    echo ""
else
    echo "  config/environment.yml already exists — skipping."
fi

# ── 7. Alembic migrations + validation ───────────────────────────────────────
echo "[7/7] Running database migrations and validating tools..."
cd "$REPO_ROOT"
alembic upgrade head

echo ""
echo "  Tool check:"
check_tool() {
    local tool="$1"
    if command -v "$tool" &>/dev/null; then
        echo "    ✓  $tool"
    else
        echo "    ✗  $tool — NOT FOUND (check conda env or brew)"
    fi
}
check_tool orthofinder
check_tool diamond
check_tool mafft
check_tool iqtree2
check_tool hyphy
check_tool fpocket
check_tool python

echo ""
echo "============================================"
echo " Setup complete!"
echo ""
echo " Each new terminal session, run:"
echo "   conda activate $ENV_NAME"
echo ""
echo " Then start the test with:"
echo "   bash scripts/run_validation_test.sh"
echo "============================================"
