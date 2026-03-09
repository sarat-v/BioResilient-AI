#!/usr/bin/env python3
"""
BioResilient Test Runner — 2 Species (Human, Naked Mole Rat, Axolotl)
Run this to test the full pipeline end-to-end locally
"""

import os
import sys
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

def main():
    repo_root = Path(__file__).parent
    os.chdir(repo_root)
    
    # Set test database
    os.environ['DATABASE_URL'] = 'postgresql://localhost/bioresillient_test'
    os.environ['NCBI_API_KEY'] = os.environ.get('NCBI_API_KEY', '')
    os.environ['ALPHAGENOME_API_KEY'] = os.environ.get('ALPHAGENOME_API_KEY', '')
    
    print("\n" + "="*60)
    print("BioResilient Test Pipeline")
    print("Species: Human, Naked Mole Rat, African Elephant")
    print("Theme: Converging Cancer Resistance")
    print("="*60 + "\n")
    
    # Import after path setup
    from pipeline.orchestrator import run_pipeline
    from db.session import SessionLocal
    from db.models import Species
    
    # Verify species exist
    session = SessionLocal()
    species = session.query(Species).all()
    session.close()
    
    if not species:
        print("❌ No species in database. Run setup script first:")
        print("   bash run_test_pipeline.sh")
        sys.exit(1)
    
    print(f"✓ Found {len(species)} species: {', '.join([s.name for s in species])}")
    print()
    
    # Run pipeline
    try:
        run_pipeline()
        print("\n" + "="*60)
        print("✅ Pipeline Complete!")
        print("="*60)
        print("\nTo view results:")
        print("  1. Start API: uvicorn api.main:app --host 127.0.0.1 --port 8000")
        print("  2. Start Frontend: cd frontend && npx vite preview --port 5173")
        print("  3. Visit: http://localhost:5173")
        print()
    except Exception as e:
        print(f"\n❌ Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
