#!/usr/bin/env bash
# BioResilient AI — AWS cloud environment bootstrap
# Tested on Ubuntu 22.04 AMI (ami-0c7217cdde317cfec)
# Run on: c5.4xlarge (CPU tasks) or g4dn.xlarge (GPU tasks)
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
echo "[3/7] Creating conda environment..."
conda env update --file "$REPO_ROOT/environment.yml" --prune
conda activate "$ENV_NAME"

# ── 4. App config ─────────────────────────────────────────────────────────────
echo "[4/7] Writing cloud config..."
cat > "$REPO_ROOT/config/environment.yml" <<EOF
deployment: cloud

database:
  local:  postgresql://localhost/bioresillient
  cloud:  postgresql://${RDS_HOST}/bioresillient

storage:
  local:  ./data/
  cloud:  s3://${S3_BUCKET}/

gpu:
  device: auto
  esm_model: esm2_t33_650M_UR50D
  esm_chunk_size: 64

ncbi:
  api_key: ${NCBI_API_KEY}
  email: ${NCBI_EMAIL}

tools:
  orthofinder_threads: 16
  orthofinder_align_threads: 8
  mafft_threads: 8
  iqtree_threads: AUTO
  iqtree_bootstrap: 1000
  hyphy_threads: 8

scoring_weights:
  convergence:      0.25
  selection:        0.20
  disease:          0.20
  druggability:     0.15
  expression:       0.10
  safety:           0.10

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

# ── 7. Tool validation ────────────────────────────────────────────────────────
echo "[7/7] Validating tools..."
for tool in orthofinder diamond mafft iqtree2 hyphy fpocket; do
    if command -v "$tool" &>/dev/null; then
        echo "  ✓  $tool"
    else
        echo "  ✗  $tool NOT FOUND"
    fi
done

python -c "import torch; print('  GPU:', torch.cuda.is_available())"

echo ""
echo "============================================"
echo " Cloud setup complete."
echo " Run: conda activate $ENV_NAME"
echo " Then: python pipeline/orchestrator.py"
echo "============================================"
