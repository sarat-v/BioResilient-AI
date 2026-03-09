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
    from db.session import get_session
    from db.models import Species
    
    # Verify species exist
    with get_session() as session:
        species = session.query(Species).all()
    
    if not species:
        print("❌ No species in database. Run setup script first:")
        print("   bash run_test_pipeline.sh")
        sys.exit(1)
    
    print(f"✓ Found {len(species)} species: {', '.join([s.scientific_name for s in species])}")
    print()

    # Build species_list dict as expected by run_pipeline
    species_list = [
        {
            "id": s.id,
            "taxid": s.taxid,
            "scientific_name": s.scientific_name,
            "lineage_group": s.lineage_group,
            "phenotypes": s.phenotypes or [],
        }
        for s in species
    ]
    
    # Run pipeline — auto-detect furthest completed step and resume from there
    try:
        import json
        cache_file = Path('pipeline_cache.json')
        completed_steps = []
        if cache_file.exists():
            completed_steps = json.loads(cache_file.read_text()).get('completed', [])

        step_order = ['step1','step2','step3','step3b','step4','step5','step6',
                      'step7','step8','step9','step10','step10b','step11',
                      'step12','step13','step14','step15','step16']

        # Treat step1 and step2 as done if proteomes already exist on disk
        proteomes_done = (Path('data/proteomes/human.reheadered.faa').exists() and
                         Path('data/proteomes/nmr.reheadered.faa').exists() and
                         Path('data/proteomes/elephant.reheadered.faa').exists())
        if proteomes_done:
            for s in ('step1', 'step2'):
                if s not in completed_steps:
                    completed_steps.append(s)

        # Find first step not yet completed
        resume_from = step_order[-1]
        for s in step_order:
            if s not in completed_steps:
                resume_from = s
                break

        print(f"▶ Resuming from {resume_from} (completed: {completed_steps or 'none'})")
        run_pipeline(species_list=species_list, resume_from=resume_from)
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
