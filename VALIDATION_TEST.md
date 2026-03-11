# Pipeline Validation Test Guide

## Quick Test Configuration

**Purpose**: Validate that all 30 pipeline steps execute correctly without waiting 2-4 days

**Test Panel**: 5 species
- 3 resilient species (3 independent lineages):
  - Naked mole rat (Rodents) - cancer resistance, longevity
  - Bowhead whale (Cetaceans) - longevity, cancer resistance  
  - African elephant (Proboscideans) - cancer resistance (TP53)
- 1 control species: Rat
- 1 reference: Human

**Expected Runtime on Your PC**: 2-4 hours

## How to Run

### Option 1: Automated Test Script (Recommended)
```bash
cd /Users/saratvakkalanka/Desktop/CoworkFiles/BioResilient/claude-code
bash scripts/run_validation_test.sh
```

This script will:
1. Backup your original `species_registry.json`
2. Switch to the 5-species test configuration
3. Run the full pipeline with all 30 steps
4. Restore your original configuration
5. Report success/failure

### Option 2: Manual Test
```bash
# Activate environment
conda activate bioresillient

# Backup original registry
cp config/species_registry.json config/species_registry.json.backup

# Switch to test registry
cp config/species_registry_test.json config/species_registry.json

# Seed database
python db/seed.py

# Set thread counts for your CPU
export ORTHOFINDER_THREADS=10
export MAFFT_THREADS=10
export IQTREE_THREADS=10
export HYPHY_THREADS=10

# Run pipeline
python pipeline/orchestrator.py

# Restore original registry
mv config/species_registry.json.backup config/species_registry.json
```

## What to Check After the Test

### 1. Verify All Steps Completed
```bash
# Check pipeline state
cat pipeline_state.json | jq '.steps | to_entries[] | select(.value.status != "complete")'

# Should return empty if all steps passed
```

### 2. Check Candidate Results
```bash
# Start API server
uvicorn api.main:app --host 0.0.0.0 --port 8000

# Query candidates
curl http://localhost:8000/candidates | jq '.[0:5]'
```

Expected: You should see some candidates with scores, even if low (convergence requires ≥2 lineages by default)

### 3. Verify Key Outputs

**Phase 1 outputs to check:**
- `data/proteomes/` - 5 FASTA files downloaded
- `data/orthofinder/Results_*/` - OrthoFinder completed
- `data/phylo/species.treefile` - IQ-TREE phylogenetic tree
- Database has populated tables:
  ```sql
  SELECT COUNT(*) FROM gene;           -- Should have genes
  SELECT COUNT(*) FROM ortholog;       -- Should have orthologs
  SELECT COUNT(*) FROM divergent_motif; -- Should have motifs
  SELECT COUNT(*) FROM candidate_score; -- Should have scores
  ```

**Known limitations with 3 species:**
- Convergence scores will be lower (requires ≥2 independent lineages)
- Many genes may fall into Tier3 (that's expected)
- You're testing the **pipeline mechanics**, not the scientific results

### 4. Check for the Fixed Bugs

Verify the 12 previously-skipped steps now execute:
```bash
grep -E "Step (4b|4c|4d|6b|6c|7b|8b|11b|11c|11d|12b|14b)" pipeline_test.log
```

Each should show "complete" status.

## Estimated Resource Usage

| Resource | Usage | Your PC | Status |
|----------|-------|---------|--------|
| CPU | 100% of 10 cores | 10 cores | ✅ OK |
| RAM | 8-12 GB peak | 32 GB | ✅ OK |
| Disk | ~50 GB | 3.1 TB | ✅ OK |
| Runtime | 2-4 hours | - | ✅ OK |

## Timeline Breakdown

| Step | Time | Cumulative |
|------|------|------------|
| Download proteomes | 5-10 min | 10 min |
| OrthoFinder | 30-60 min | 1h 10m |
| Alignment & divergence | 20-30 min | 1h 40m |
| Phylogenetic tree | 10-15 min | 2h |
| MEME selection | 30-60 min | 3h |
| FEL + BUSTED | 20-30 min | 3h 30m |
| RELAX | 15-20 min | 4h |
| Convergence & expression | 10-15 min | 4h 15m |
| Phase 1 scoring | 1-2 min | 4h 17m |
| Phase 2 steps | 15-30 min | 4h 45m |

**Total: ~2-5 hours** (depends on network speed for NCBI downloads)

## Troubleshooting

### If OrthoFinder fails
- Check NCBI proteome downloads completed: `ls -lh data/proteomes/`
- Verify no empty FASTA files: `wc -l data/proteomes/*.faa`

### If HyPhy fails
- Check that codon alignments exist
- Verify IQ-TREE species tree was created: `ls -lh data/phylo/species.treefile`

### If memory issues occur
Reduce parallelism:
```bash
export ORTHOFINDER_THREADS=5
export HYPHY_THREADS=5
```

### If steps are still being skipped
Check the orchestrator was updated:
```bash
grep -c "elif step_name ==" pipeline/orchestrator.py
# Should output: 29
```

## After Validation Passes

1. **Review the fixes**: Verify all sections now appear in the frontend
2. **Run benchmark** (will likely fail with 3 species, but check it runs):
   ```bash
   python scripts/benchmark_recall.py
   ```
3. **Decide next steps**:
   - Run full 19-species panel on your PC (48-96 hours)
   - Switch to AWS for production runs
   - Iterate on specific components

## Notes

- This test uses **real NCBI data** (not mocked), so results are scientifically valid for these 3 species
- The test proves the **pipeline mechanics work end-to-end**
- With only 3 species, you won't see strong convergence signals (need ≥5-8 for that)
- Perfect for validating that the 12 bug fixes we made actually work

---

**Ready to run?**
```bash
bash scripts/run_validation_test.sh
```

Monitor progress in another terminal:
```bash
tail -f pipeline_test.log
```
