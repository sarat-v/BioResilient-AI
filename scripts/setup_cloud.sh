#!/usr/bin/env bash
# BioResilient AI — AWS cloud environment bootstrap
# Tested on Ubuntu 22.04 AMI (ami-0c7217cdde317cfec)
# Run on: c6i.4xlarge (CPU tasks) or g4dn.xlarge (GPU tasks)
# Expected cost: ~$0.27/hr (CPU spot) | ~$0.16/hr (GPU spot)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="bioresillient"

echo "============================================"
echo " BioResilient AI — AWS Cloud Setup"
echo "============================================"

# ── Required environment variables ───────────────────────────────────────────
: "${NCBI_API_KEY:?NCBI_API_KEY environment variable not set}"
: "${NCBI_EMAIL:?NCBI_EMAIL environment variable not set}"
: "${RDS_HOST:?RDS_HOST environment variable not set (RDS endpoint)}"
: "${S3_BUCKET:=bioresillient-data}"
# Optional API keys (pipeline degrades gracefully without them)
ALPHAGENOME_KEY="${ALPHAGENOME_API_KEY:-}"
ANTHROPIC_KEY="${ANTHROPIC_API_KEY:-}"

echo "  RDS host : $RDS_HOST"
echo "  S3 bucket: $S3_BUCKET"

# ── 1. System packages ────────────────────────────────────────────────────────
echo "[1/7] Installing system packages..."
sudo apt-get update -qq
sudo apt-get install -y wget curl git postgresql-client build-essential

# ── 2. Conda ─────────────────────────────────────────────────────────────────
if ! command -v conda &>/dev/null; then
    echo "[2/7] Installing Miniconda..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
         -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
    eval "$("$HOME/miniconda3/bin/conda" shell.bash hook)"
    conda init bash
else
    echo "[2/7] Conda already installed — skipping."
fi

eval "$(conda shell.bash hook)"

# ── 3. Conda environment ─────────────────────────────────────────────────────
echo "[3/7] Creating conda environment (includes minimap2 + PHAST)..."
conda env update --file "$REPO_ROOT/environment.yml" --prune
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
conda activate "$ENV_NAME"

# Install tools that need bioconda and may not be in environment.yml
conda install -y -c bioconda minimap2   2>/dev/null || true
conda install -y -c bioconda phast      2>/dev/null || true   # provides phyloP + phastCons
conda install -y -c bioconda lastz      2>/dev/null || true   # fallback aligner for step3c

# ── 4. App config ─────────────────────────────────────────────────────────────
echo "[4/7] Writing cloud config..."
# This mirrors the exact YAML structure that pipeline/config.py expects.
# Keys must stay in sync with config.py:get_tool_config() and get_db_url().
cat > "$REPO_ROOT/config/environment.yml" <<EOF
# BioResilient AI — Application Configuration (LIVE SECRETS FILE)
# DO NOT COMMIT THIS FILE — it is in .gitignore
deployment: cloud

database:
  local:  postgresql://localhost/bioresillient
  cloud:  postgresql://${RDS_HOST}/bioresillient

storage:
  local:  ./data/
  cloud:  s3://${S3_BUCKET}/

gpu:
  device: auto
  esm1v_model: esm1v_t33_650M_UR90S_1
  esm_chunk_size: 128

target_phenotype: cancer_resistance

ncbi:
  api_key: ${NCBI_API_KEY}
  email: ${NCBI_EMAIL}

alphagenome:
  api_key: ${ALPHAGENOME_KEY}

anthropic:
  api_key: ${ANTHROPIC_KEY}

tools:
  local:
    orthofinder_threads: 6
    orthofinder_align_threads: 3
    mafft_threads: 6
    iqtree_threads: AUTO
    iqtree_bootstrap: 1000
    hyphy_threads: 4
    minimap2_threads: 6
  cloud:
    orthofinder_threads: 14
    orthofinder_align_threads: 7
    mafft_threads: 14
    iqtree_threads: AUTO
    iqtree_bootstrap: 1000
    hyphy_threads: 14
    minimap2_threads: 14

thresholds:
  divergence_identity_max: 0.85
  divergence_min_species: 2
  convergence_min_lineages: 3
  expression_log2fc_min: 1.0
  expression_padj_max: 0.05
  tier1_composite_min: 0.70
  tier2_composite_min: 0.40
EOF

# ── 5. Alembic migrations + seed ─────────────────────────────────────────────
echo "[5/7] Running Alembic migrations against RDS..."
cd "$REPO_ROOT"
alembic upgrade head
python db/seed.py

# ── 6. S3 data directory structure ────────────────────────────────────────────
echo "[6/7] Creating S3 data structure..."
aws s3api put-object --bucket "$S3_BUCKET" --key "proteomes/"
aws s3api put-object --bucket "$S3_BUCKET" --key "orthofinder_out/"
aws s3api put-object --bucket "$S3_BUCKET" --key "alignments/"
aws s3api put-object --bucket "$S3_BUCKET" --key "geo_data/"
aws s3api put-object --bucket "$S3_BUCKET" --key "genomes/"            # step3c: genomic FASTA cache
aws s3api put-object --bucket "$S3_BUCKET" --key "nucleotide_regions/" # step3c: extracted region sequences
aws s3api put-object --bucket "$S3_BUCKET" --key "cache/"              # step3/4/5: pkl/treefile cache

# ── 7. Tool validation ────────────────────────────────────────────────────────
echo "[7/7] Validating tools..."
all_ok=true
for tool in orthofinder diamond mafft iqtree2 hyphy fpocket minimap2 phyloP phastCons; do
    if command -v "$tool" &>/dev/null; then
        echo "  ✓  $tool"
    else
        echo "  ✗  $tool NOT FOUND"
        # minimap2 / phyloP / phastCons are optional (pipeline falls back gracefully)
        if [[ "$tool" != "minimap2" && "$tool" != "phyloP" && "$tool" != "phastCons" ]]; then
            all_ok=false
        fi
    fi
done

python -c "import torch; print('  GPU:', torch.cuda.is_available())"

if [[ "$all_ok" != "true" ]]; then
    echo ""
    echo "⚠  One or more required tools are missing. Check the output above."
    exit 1
fi

echo ""
echo "============================================"
echo " Cloud setup complete."
echo " Run: conda activate $ENV_NAME"
echo " Then: bash run_cancer_resistance_stepwise.sh"
echo "============================================"
