# BioResilient AI — Technical Build Plan
> Claude Code engineering reference. Read alongside `BioResilient_Master_Build_Plan.docx` for scientific context.

---

## 0. Hardware Strategy

The pipeline is hardware-agnostic by design. Two valid deployment targets:

| | Local (GPU Rig) | Cloud (AWS) |
|---|---|---|
| GPU tasks (ESMFold, ESM-2) | RTX 3080 Ti / 4070 Ti (12GB+ VRAM) | `g4dn.xlarge` (T4, 16GB VRAM) |
| CPU tasks (OrthoFinder, HyPhy, IQ-TREE) | 8+ core CPU, 32GB RAM | `c5.4xlarge` (16 vCPU, 32GB RAM) |
| Database | PostgreSQL local | AWS RDS `db.t3.medium` |
| Storage | Local NVMe (500GB+) | S3 + EBS |
| OS | Ubuntu 22.04 or WSL2 (Windows) | Ubuntu 22.04 AMI |

**Decision rule:** Use local if GPU rig is available. Use cloud if not. The code is identical either way — only `config/environment.yml` and connection strings change. All expensive cloud steps use **spot instances**.

---

## 1. Repository Structure

```
bioresilient/
├── README.md
├── environment.yml              # Conda env (all tools)
├── config/
│   ├── environment.yml          # LOCAL vs CLOUD toggle + secrets
│   ├── species_registry.json    # 20-species seed list
│   └── scoring_weights.json     # Layer composite weights
├── pipeline/
│   ├── orchestrator.py          # Runs full pipeline, phase-gated
│   ├── layer1_sequence/
│   │   ├── download.py          # NCBI genome + proteome fetch
│   │   ├── orthofinder.py       # OrthoFinder wrapper (--diamond flag)
│   │   ├── alignment.py         # MAFFT wrapper
│   │   ├── divergence.py        # Sliding window, sequence identity
│   │   └── expression.py        # GEO dataset fetch + DESeq2 via rpy2
│   ├── layer2_evolution/
│   │   ├── phylo_tree.py        # IQ-TREE2 wrapper
│   │   ├── selection.py         # HyPhy aBSREL wrapper
│   │   └── convergence.py       # PhyloP scores + multi-lineage check
│   ├── layer3_disease/
│   │   ├── opentargets.py       # REST API calls
│   │   ├── gwas.py              # GWAS Catalog FTP + query
│   │   ├── gnomad.py            # gnomAD GraphQL API
│   │   ├── impc.py              # Mouse KO phenotype REST API
│   │   └── protein_atlas.py    # HPA REST API
│   ├── layer4_druggability/
│   │   ├── structure.py         # AlphaFold DB download + ESMFold API
│   │   ├── pockets.py           # fpocket wrapper
│   │   ├── chembl.py            # ChEMBL REST API
│   │   ├── cansar.py            # CanSAR API
│   │   └── peptide.py           # Tractability filters on motifs
│   ├── layer5_gene_therapy/
│   │   ├── aav.py               # Gene size check, tissue tropism map
│   │   └── crispr.py            # CRISPOR + Cas-OFFinder wrappers
│   ├── layer6_safety/
│   │   ├── phewas.py            # PheWAS association lookup
│   │   ├── network.py           # STRING centrality scoring
│   │   └── selectivity.py       # Protein family off-target risk
│   └── scoring.py               # Composite score assembly
├── db/
│   ├── models.py                # SQLAlchemy ORM (schema below)
│   ├── migrations/              # Alembic
│   └── seed.py                  # Loads species_registry.json
├── api/
│   ├── main.py                  # FastAPI app
│   └── routes/
│       ├── candidates.py        # GET /candidates, /candidates/{id}
│       ├── species.py           # GET /species
│       └── scores.py            # GET /scores/{gene_id}
├── scripts/
│   ├── setup_local.sh           # Local environment bootstrap
│   └── setup_cloud.sh           # AWS environment bootstrap
└── tests/
    ├── test_layer1.py
    ├── test_layer2.py
    └── test_scoring.py
```

---

## 2. Environment Setup

### Option A — Local (WSL2 or native Ubuntu)

```bash
# 1. Install Conda (if not present)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# 2. Create environment
conda env create -f environment.yml
conda activate bioresilient

# 3. Install CUDA (if using GPU locally — skip on cloud, AMI includes it)
# Follow: https://developer.nvidia.com/cuda-downloads (CUDA 12.x)

# 4. WSL2 only — set memory limits in %UserProfile%\.wslconfig
# [wsl2]
# memory=28GB
# processors=8

# 5. Init database
createdb bioresilient
alembic upgrade head
python db/seed.py
```

### Option B — AWS Cloud

```bash
# Recommended: use a single bootstrap script
bash scripts/setup_cloud.sh

# GPU jobs: g4dn.xlarge spot (~$0.16/hr)  → ESMFold, ESM-2
# CPU jobs: c5.4xlarge spot (~$0.27/hr)   → OrthoFinder, HyPhy, IQ-TREE
# DB:       RDS db.t3.medium (~$50/month) → PostgreSQL
# Storage:  S3 standard (~$0.023/GB/mo)   → genome files
```

### `environment.yml`

```yaml
name: bioresilient
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - postgresql
  - blast=2.15
  - mafft=7.520
  - orthofinder=2.5.5       # includes diamond
  - diamond=2.1.8
  - iqtree=2.2.6
  - hyphy=2.5.62
  - fpocket=4.0
  - r-base=4.3
  - bioconductor-deseq2
  - bioconductor-edger
  - pip
  - pip:
    - sqlalchemy==2.0
    - alembic==1.13
    - fastapi==0.110
    - uvicorn==0.29
    - biopython==1.83
    - geoparse==2.0
    - rpy2==3.5
    - torch==2.2          # GPU: auto-detects CUDA; CPU: runs slower
    - fair-esm==2.0       # ESM-2 embeddings
    - requests==2.31
    - pydantic==2.6
    - pandas==2.2
    - numpy==1.26
    - pytest==8.1
```

---

## 3. Configuration

### `config/environment.yml` (not the Conda one — app config)

```yaml
# Toggle: "local" or "cloud"
deployment: local

database:
  local:  postgresql://localhost/bioresilient
  cloud:  postgresql://${RDS_HOST}/bioresilient

storage:
  local:  ./data/
  cloud:  s3://bioresilient-data/

gpu:
  device: cuda          # "cuda" or "cpu" — auto-detected if omitted
  esm_model: esm2_t33_650M_UR50D   # 650M fits all VRAM configs
  # esm2_t36_3B_UR50D             # 3B: use only if 16GB+ VRAM available

ncbi:
  api_key: ${NCBI_API_KEY}         # Free key, increases rate limits
  email: ${NCBI_EMAIL}

scoring_weights:
  convergence:      0.25
  selection:        0.20
  disease:          0.20
  druggability:     0.15
  expression:       0.10
  safety:           0.10
```

### `config/species_registry.json` (Phase 1 seed)

```json
[
  {"id": "naked_mole_rat",     "taxid": 10181,  "name": "Heterocephalus glaber",       "phenotype": ["cancer_resistance", "longevity"],          "genome_assembly": "GCF_000247695.1"},
  {"id": "bowhead_whale",      "taxid": 13397,  "name": "Balaena mysticetus",           "phenotype": ["longevity", "cancer_resistance"],          "genome_assembly": "GCF_028564815.1"},
  {"id": "axolotl",            "taxid": 8296,   "name": "Ambystoma mexicanum",          "phenotype": ["regeneration"],                           "genome_assembly": "GCF_002915635.3"},
  {"id": "ground_squirrel",    "taxid": 43179,  "name": "Ictidomys tridecemlineatus",   "phenotype": ["hibernation", "ischaemia_resistance"],     "genome_assembly": "GCF_000236235.1"},
  {"id": "little_brown_bat",   "taxid": 59463,  "name": "Myotis lucifugus",             "phenotype": ["viral_resistance", "longevity"],           "genome_assembly": "GCF_000147115.1"},
  {"id": "damaraland_mole_rat","taxid": 885580, "name": "Fukomys damarensis",           "phenotype": ["longevity", "cancer_resistance"],          "genome_assembly": "GCF_000743615.1"},
  {"id": "greenland_shark",    "taxid": 58051,  "name": "Somniosus microcephalus",      "phenotype": ["longevity"],                              "genome_assembly": "GCF_963853765.1"},
  {"id": "african_elephant",   "taxid": 9785,   "name": "Loxodonta africana",           "phenotype": ["cancer_resistance"],                      "genome_assembly": "GCF_000001905.1"},
  {"id": "mouse_lemur",        "taxid": 30608,  "name": "Microcebus murinus",           "phenotype": ["hibernation"],                            "genome_assembly": "GCF_000165445.2"},
  {"id": "spiny_mouse",        "taxid": 10103,  "name": "Acomys cahirinus",             "phenotype": ["regeneration"],                           "genome_assembly": "GCF_009650055.2"},
  {"id": "bowhead_rockfish",   "taxid": 161559, "name": "Sebastes aleutianus",          "phenotype": ["longevity"],                              "genome_assembly": "GCF_021917145.1"},
  {"id": "human",              "taxid": 9606,   "name": "Homo sapiens",                 "phenotype": ["baseline"],                               "genome_assembly": "GCF_000001405.40"}
]
```

---

## 4. Database Schema

```python
# db/models.py — SQLAlchemy ORM

class Species(Base):
    id              str PK
    taxid           int
    scientific_name str
    phenotypes      ARRAY[str]
    genome_assembly str
    proteome_path   str          # local or S3 path

class Gene(Base):
    id              UUID PK
    human_gene_id   str          # NCBI Gene ID
    gene_symbol     str
    human_protein   str          # UniProt accession

class Ortholog(Base):
    id              UUID PK
    gene_id         UUID FK(Gene)
    species_id      str  FK(Species)
    protein_seq     text
    sequence_identity_pct  float  # vs human
    orthofinder_og  str          # OrthoGroup ID

class DivergentMotif(Base):
    id              UUID PK
    ortholog_id     UUID FK(Ortholog)
    start_pos       int
    end_pos         int
    animal_seq      str          # 10-20 aa window
    human_seq       str
    divergence_score float
    esm_distance    float        # embedding cosine distance
    # Peptide tractability (Layer 4)
    half_life_min   float
    logp            float
    immunogenic     bool
    synthesisable   bool

class EvolutionScore(Base):
    gene_id         UUID FK(Gene)
    dnds_ratio      float
    dnds_pvalue     float
    selection_model str          # aBSREL, BUSTED etc
    convergence_count int        # number of independent lineages
    phylop_score    float

class DiseaseAnnotation(Base):
    gene_id         UUID FK(Gene)
    disease_id      str          # EFO ID
    disease_name    str
    opentargets_score float
    gwas_pvalue     float
    gnomad_pli      float        # constraint score
    mouse_ko_phenotype str
    tissue_expression JSON       # {tissue: tpm_value}

class DrugTarget(Base):
    gene_id         UUID FK(Gene)
    pocket_count    int
    top_pocket_score float
    chembl_target_id str
    existing_drugs  ARRAY[str]
    cansar_score    float
    druggability_tier str        # A/B/C/undruggable

class GeneTherapyScore(Base):
    gene_id         UUID FK(Gene)
    gene_size_bp    int
    aav_compatible  bool
    tissue_tropism  ARRAY[str]   # matching AAV serotypes
    crispr_sites    int          # number of valid guide sites
    offtarget_risk  str          # low/medium/high

class SafetyFlag(Base):
    gene_id         UUID FK(Gene)
    is_essential    bool         # pLI > 0.9
    phewas_hits     JSON         # {trait: pvalue}
    network_degree  int          # STRING interaction count
    hub_risk        bool         # degree > 50
    family_size     int          # protein family size

class CandidateScore(Base):
    gene_id         UUID FK(Gene)   PK
    convergence_score  float
    selection_score    float
    disease_score      float
    druggability_score float
    expression_score   float
    safety_score       float
    composite_score    float
    tier               str       # Tier1 / Tier2 / Tier3
    updated_at         datetime
```

---

## 5. Phase 1 Implementation Order

Build in this exact sequence. Each step depends on the previous.

### Step 1 — Environment + Database (Day 1–2)
- Run `setup_local.sh` or `setup_cloud.sh`
- Validate all tools: `orthofinder --version`, `iqtree --version`, `hyphy --version`, `mafft --version`, `fpocket --version`
- Init PostgreSQL, run Alembic migrations, load species seed
- Verify GPU (local): `python -c "import torch; print(torch.cuda.is_available())"`

### Step 2 — Data Download (Day 3–5)
- `pipeline/layer1_sequence/download.py`
- For each species in registry: fetch reference proteome FASTA from NCBI using `taxid`
- Store to `./data/proteomes/{species_id}.faa` or `s3://bioresilient-data/proteomes/`
- Also download human reference proteome (UniProt canonical, ~20k proteins)
- Validate: protein count per species, no corrupt FASTA files

### Step 3 — Ortholog Finding (Day 6–9)
- `pipeline/layer1_sequence/orthofinder.py`
- Run OrthoFinder with `--diamond` flag on all proteomes simultaneously
- **Critical flag:** `-S diamond` (uses DIAMOND, not BLAST — 500x faster, less RAM)
- Output: `OrthoGroups.tsv` — one row per ortholog group, one column per species
- Parse output, load into `Ortholog` table
- Filter: keep only orthogroups where human protein is present

```bash
# OrthoFinder command (generated by wrapper)
orthofinder -f ./data/proteomes/ -S diamond -t 8 -a 4 -o ./data/orthofinder_out/
```

### Step 4 — Sequence Alignment + Divergence (Day 10–12)
- `pipeline/layer1_sequence/alignment.py` + `divergence.py`
- For each ortholog group: run MAFFT on the protein sequences
- Calculate whole-protein sequence identity % (human vs each animal)
- Run sliding window (size=15, step=5) to find divergent motifs
- Score each motif: divergence % + ESM-2 embedding cosine distance
- Load `DivergentMotif` table

```bash
# MAFFT command (per ortholog group)
mafft --auto --thread 4 input.faa > aligned.faa
```

### Step 5 — Phylogenetic Tree (Day 13–14)
- `pipeline/layer2_evolution/phylo_tree.py`
- Build species tree using IQ-TREE on concatenated alignment of single-copy orthologs
- **Required before HyPhy** — HyPhy needs a tree
- Model selection: `-m TEST` (auto-selects best substitution model)
- Output: `species.treefile` — Newick format

```bash
iqtree2 -s concat_alignment.faa -m TEST -bb 1000 -T AUTO -o human
```

### Step 6 — Evolutionary Selection (Day 15–19)
- `pipeline/layer2_evolution/selection.py`
- Run HyPhy aBSREL on each candidate ortholog group (not full proteome — only those passing divergence threshold from Step 4)
- Threshold: sequence identity < 85% in at least 2 protective species
- aBSREL tests for episodic positive selection on each branch
- Extract: dN/dS ratio, p-value, branches under selection
- Load into `EvolutionScore`

```bash
# HyPhy command (per gene, generated by wrapper)
hyphy aBSREL --alignment aligned.faa --tree species.treefile --output results.json
```

### Step 7 — Convergence Detection (Day 20–22)
- `pipeline/layer2_evolution/convergence.py`
- Query UCSC API for PhyloP scores at each divergent motif position
- Count how many independent evolutionary lineages show the same directional change
- Lineage groups (treat each as independent): Rodents, Cetaceans, Bats, Sharks, Primates, Salamanders
- A motif scoring in 3+ independent lineages = convergence flag
- Load `convergence_count` into `EvolutionScore`

### Step 8 — Expression Annotation (Day 23–26)
- `pipeline/layer1_sequence/expression.py`
- For each species with GEO data available (ground squirrel, naked mole rat, axolotl priority)
- Search GEO for datasets matching species + protective condition (torpor, tumour challenge, injury)
- Download raw counts, run DESeq2 via rpy2
- Flag candidate genes as upregulated (log2FC > 1, padj < 0.05) or not
- Load `expression_score` into `CandidateScore`

### Step 9 — Composite Score (Day 27–28)
- `pipeline/scoring.py`
- Assemble `CandidateScore` from all Layer 1 + Layer 2 outputs
- Apply weights from `config/scoring_weights.json`
- Assign tiers: Tier 1 (composite > 0.7), Tier 2 (0.4–0.7), Tier 3 (< 0.4)
- **Phase 1 complete at this point**

### Step 10 — Basic API (Day 29–30)
- `api/main.py` + routes
- Three endpoints sufficient for Phase 1:
  - `GET /candidates?tier=1` — ranked list with all scores
  - `GET /candidates/{gene_id}` — full detail including motifs
  - `GET /species` — species registry
- No auth for internal use. Add API key middleware before any external access.

---

## 6. Key Tool Flags Reference

| Tool | Critical Flag | Reason |
|---|---|---|
| OrthoFinder | `-S diamond` | DIAMOND instead of BLAST — 500x faster, less RAM |
| OrthoFinder | `-t 8 -a 4` | 8 search threads, 4 alignment threads — tune to available cores |
| IQ-TREE2 | `-m TEST` | Auto-selects best substitution model — don't hardcode |
| IQ-TREE2 | `-bb 1000` | 1000 ultrafast bootstrap replicates for branch support values |
| HyPhy | `aBSREL` | Branch-site model — tests episodic selection per branch, not whole tree |
| MAFFT | `--auto` | Auto-selects algorithm based on input size — correct for protein sets |
| ESMFold | `chunk_size=64` | Reduces VRAM usage for long proteins — set if hitting OOM errors |
| fpocket | `-m 3.0` | Minimum pocket volume (Å³) — reduces noise from tiny surface pockets |

---

## 7. Data Formats Between Steps

| From → To | Format | Key Fields |
|---|---|---|
| Download → OrthoFinder | FASTA `.faa` | `>{species_id}|{protein_id}` header |
| OrthoFinder → DB | TSV `OrthoGroups.tsv` | OG ID, one accession column per species |
| Alignment → Divergence | FASTA `.afa` (aligned) | Same headers, gap characters inserted |
| Divergence → ESMFold | FASTA `.faa` | Animal motif sequences only, 15aa windows |
| HyPhy → DB | JSON `.json` | `test results` → branch dN/dS, p-values |
| IQ-TREE → HyPhy | Newick `.treefile` | Species tree with branch labels |
| GEO → DESeq2 | Counts matrix `.csv` | Gene × sample, raw integer counts |
| All layers → Score | PostgreSQL | Join on `gene_id` across all score tables |

---

## 8. Testing Strategy

Each layer module has a corresponding test. Tests use a **fixture dataset of 3 species** (human, naked mole rat, ground squirrel) and 100 genes — runs in under 5 minutes on any hardware.

```bash
pytest tests/ -v                    # all tests
pytest tests/test_layer1.py -v      # layer 1 only
pytest tests/ -k "not slow"         # skip full-proteome integration tests
```

Slow integration tests (full 20-species run) are marked `@pytest.mark.slow` and run manually before phase completion, not in CI.

---

## 9. What Phase 1 Does NOT Include

These are intentionally excluded and belong to Phase 2. Do not build them now:
- OpenTargets / GWAS / gnomAD annotation (Layer 3)
- AlphaFold / ESMFold structure prediction (Layer 4)
- fpocket / ChEMBL / CanSAR (Layer 4)
- Gene therapy assessment (Layer 5)
- Safety pre-screening (Layer 6)
- Frontend UI
- Any external API authentication

Phase 1 ends when: the pipeline runs end-to-end on all 12 species without intervention, `CandidateScore` is populated for all genes passing Layer 1 threshold, and dN/dS values with p-values exist for all Tier 1 + Tier 2 candidates.

---

## 10. Reference Databases (All Free)

| Database | Access | Used In |
|---|---|---|
| NCBI RefSeq | E-utilities API + FTP | Step 2 |
| UniProt | REST API | Step 2 |
| NCBI GEO | FTP + GEOparse | Step 8 |
| UCSC PhyloP | REST API | Step 7 |
| OpenTargets | REST API (no key) | Phase 2 |
| GWAS Catalog | FTP download | Phase 2 |
| gnomAD | GraphQL API | Phase 2 |
| IMPC | REST API | Phase 2 |
| Human Protein Atlas | REST API | Phase 2 |
| AlphaFold DB | FTP download | Phase 2 |
| ESMFold | API or local model | Phase 2 |
| ChEMBL | REST API | Phase 2 |
| CanSAR | REST API (registration) | Phase 2 |
| STRING | REST API | Phase 2 |
| PDB | REST API + FTP | Phase 2 |
