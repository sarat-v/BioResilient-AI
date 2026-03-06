# BioResilient AI

A bioinformatics pipeline that identifies therapeutic target candidates by mining the genomes of naturally disease-resistant animals. The platform identifies genes under positive evolutionary selection across multiple independent resilient lineages, prioritising candidates where the same functional motif has converged in geographically and phylogenetically distant species.

## Architecture

Six-layer scoring system:

| Layer | Description | Phase |
|---|---|---|
| Layer 1 | Sequence divergence — ortholog finding, alignment, motif scoring | Phase 1 |
| Layer 2 | Evolutionary selection — dN/dS, convergence detection | Phase 1 |
| Layer 3 | Disease annotation — OpenTargets, GWAS, gnomAD | Phase 2 |
| Layer 4 | Druggability — structure, pockets, ChEMBL | Phase 2 |
| Layer 5 | Gene therapy feasibility — AAV compatibility, CRISPR sites | Phase 2 |
| Layer 6 | Safety pre-screen — PheWAS, network centrality | Phase 2 |

## Phase 1 (current)

Runs end-to-end on 12 species (11 resilient + human baseline), producing Layer 1 + Layer 2 scores. Outputs a ranked candidate list with composite scores and tier assignments.

## Quick Start

### Local (Ubuntu / WSL2)

```bash
bash scripts/setup_local.sh
conda activate bioresillient
createdb bioresillient
alembic upgrade head
python db/seed.py
python pipeline/orchestrator.py
```

### Cloud (AWS)

```bash
bash scripts/setup_cloud.sh
# Remaining steps identical to local after bootstrap
```

### Run tests

```bash
pytest tests/ -v
pytest tests/ -k "not slow"   # skip full-proteome integration tests
```

## Species (Phase 1 seed — 12 species)

| Species | Taxon ID | Key Phenotype |
|---|---|---|
| Naked mole rat | 10181 | Cancer resistance, longevity |
| Bowhead whale | 13397 | Longevity, cancer resistance |
| Axolotl | 8296 | Regeneration |
| 13-lined ground squirrel | 43179 | Hibernation, ischaemia resistance |
| Little brown bat | 59463 | Viral resistance, longevity |
| Damaraland mole rat | 885580 | Longevity, cancer resistance |
| Greenland shark | 58051 | Longevity |
| African elephant | 9785 | Cancer resistance |
| Mouse lemur | 30608 | Hibernation |
| Spiny mouse | 10103 | Regeneration |
| Rougheye rockfish | 161559 | Longevity |
| Human | 9606 | Baseline |

## API Endpoints (Phase 1)

```
GET /candidates?tier=1          Ranked candidates filtered by tier
GET /candidates/{gene_id}       Full detail: scores, motifs, orthologs
GET /species                    Species registry
```

## Environment Requirements

- Python 3.11 via Conda
- PostgreSQL 14+
- OrthoFinder 2.5.5 (with DIAMOND)
- MAFFT 7.520
- IQ-TREE2 2.2.6
- HyPhy 2.5.62
- GPU optional for Phase 1 (required for Phase 2 ESMFold / ESM-2)

## Configuration

Copy `config/environment.example.yml` to `config/environment.yml` and set:

```yaml
deployment: local   # or "cloud"
ncbi:
  api_key: YOUR_NCBI_KEY
  email: your@email.com
```

Set environment variables before running:

```bash
export NCBI_API_KEY=your_key
export NCBI_EMAIL=your@email.com
# Cloud only:
export RDS_HOST=your-rds-endpoint.region.rds.amazonaws.com
```
