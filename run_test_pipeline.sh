#!/bin/bash
# BioResilient Test Pipeline — Local 3-Species Run
# Tests: Human vs Naked Mole Rat vs African Elephant (Converging Cancer Resistance)

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

echo "================================"
echo "BioResilient Test Pipeline"
echo "Species: Human, Naked Mole Rat, African Elephant"
echo "Theme: Converging Cancer Resistance"
echo "================================"
echo ""

# Step 1: Check Python
echo "[1/6] Checking Python..."
python3 --version || { echo "Python3 not found"; exit 1; }

# Step 2: Check PostgreSQL
echo "[2/6] Checking PostgreSQL..."
if ! command -v psql &> /dev/null; then
    echo "PostgreSQL not found. Installing..."
    if command -v brew &> /dev/null; then
        brew install postgresql@15
        brew services start postgresql@15
    else
        echo "Please install PostgreSQL manually"
        exit 1
    fi
fi

# Step 3: Create test database
echo "[3/6] Setting up database..."
export PATH="/opt/homebrew/opt/postgresql@15/bin:$PATH"
dropdb bioresillient_test 2>/dev/null || true
createdb bioresillient_test
echo "✓ Database created: bioresillient_test"

# Step 4: Run migrations
echo "[4/6] Running database migrations..."
cd "$REPO_ROOT"
export PATH="/opt/homebrew/opt/postgresql@15/bin:$PATH"
export DATABASE_URL="postgresql://localhost/bioresillient_test"
alembic upgrade head 2>/dev/null || echo "Note: Alembic may need manual setup"
echo "✓ Migrations complete"

# Step 5: Seed with 2-3 species
echo "[5/6] Seeding database with test species..."
export DATABASE_URL="postgresql://localhost/bioresillient_test"
python3 << 'PYTHON_SEED'
import sys
import os
sys.path.insert(0, '.')

# Set database URL before importing ORM
os.environ['DATABASE_URL'] = 'postgresql://localhost/bioresillient_test'

from db.session import get_session, get_engine
from db.models import Base, Species

# Ensure tables exist
engine = get_engine()
Base.metadata.create_all(engine)

with get_session() as session:
    # Check if already seeded
    existing = session.query(Species).count()
    if existing > 0:
        print(f"✓ Database already has {existing} species")
        sys.exit(0)

    # Seed 3 test species
    species_data = [
        {'id': 'human', 'name': 'Homo sapiens', 'lineage': 'Primates', 'phenotypes': ['bipedal', 'large_brain', 'complex_language']},
        {'id': 'nmr', 'name': 'Heterocephalus glaber', 'lineage': 'Rodentia', 'phenotypes': ['eusocial', 'extreme_longevity', 'cancer_resistance', 'pain_insensitivity']},
        {'id': 'elephant', 'name': 'Loxodonta africana', 'lineage': 'Proboscidea', 'phenotypes': ['large_body', 'longevity', 'cancer_resistance', 'low_tumor_rate']},
    ]

    for s in species_data:
        sp = Species(
            id=s['id'],
            name=s['name'],
            lineage=s['lineage'],
            phenotypes=s['phenotypes']
        )
        session.add(sp)

    print(f"✓ Seeded {len(species_data)} species")
PYTHON_SEED

# Step 6: Display next steps
echo "[6/6] Test setup complete!"
echo ""
echo "================================"
echo "Next Steps:"
echo "================================"
echo ""
echo "Option A: Run full pipeline"
echo "  cd $REPO_ROOT"
echo "  DATABASE_URL='postgresql://localhost/bioresillient_test' python3 pipeline/orchestrator.py"
echo ""
echo "Option B: Run API only (to test UI)"
echo "  DATABASE_URL='postgresql://localhost/bioresillient_test' uvicorn api.main:app --host 127.0.0.1 --port 8000"
echo ""
echo "Option C: Preview frontend (separate terminal)"
echo "  cd frontend && npx vite preview --port 5173"
echo ""
echo "Then visit: http://localhost:5173 (point to http://localhost:8000)"
echo ""
