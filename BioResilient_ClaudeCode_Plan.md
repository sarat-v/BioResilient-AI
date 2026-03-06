# BioResilient AI — Engineering Build Plan

> **Purpose:** Complete, Claude Code-ready build plan for the BioResilient AI platform. Every tool has verified programmatic access. Every inter-step data format is specified. Feed this to Claude Code step by step.

---

## Repository Structure

```
bioresillient/
├── README.md
├── main.nf                          # Master Nextflow pipeline
├── nextflow.config                  # Cloud + resource profiles
├── conda.yml                        # All tool dependencies
├── docker/
│   └── Dockerfile                   # Single image with all CLI tools
├── config/
│   ├── species_registry.json        # Step 1 input seed data
│   ├── scoring_weights.json         # Step 8 scoring config
│   └── api_keys.env                 # NCBI API key etc (gitignored)
├── pipeline/
│   ├── step0_species_discovery.py   # Optional: automated candidate species finder
│   ├── step1_registry.py
│   ├── step2_download.py
│   ├── step3_qc.py
│   ├── step4_orthologs.py
│   ├── step5a_alignment.py
│   ├── step5b_codon_alignment.py
│   ├── step5c_trim.py
│   ├── step6_selection.py
│   ├── step7a_normalise_ids.py
│   ├── step7b_annotate.py
│   ├── step8_scoring.py
│   └── orchestrator.py
├── db/
│   ├── models.py                    # SQLAlchemy ORM models
│   ├── migrations/                  # Alembic migration files
│   └── seed.py
├── api/
│   ├── main.py
│   ├── routes/
│   │   ├── targets.py
│   │   ├── species.py
│   │   └── validation.py
│   └── templates/
├── nextflow/
│   └── modules/
│       ├── busco.nf
│       ├── repeatmasker.nf
│       ├── orthofinder.nf
│       ├── mafft.nf
│       ├── pal2nal.nf
│       ├── trimal.nf
│       ├── iqtree2.nf
│       └── hyphy.nf
└── tests/
    ├── test_download.py
    ├── test_orthologs.py
    ├── test_annotation.py
    └── fixtures/                    # Small test genomes (2–3 species)
```

---

## Environment Setup

```bash
conda create -n bioresillient python=3.11
conda activate bioresillient
conda install -c bioconda -c conda-forge \
    orthofinder mafft trimal iqtree hyphy \
    busco repeatmasker pal2nal macse \
    biopython boto3 sqlalchemy alembic \
    fastapi uvicorn statsmodels pandas \
    requests openpyxl reportlab nextflow

pip install python-dotenv pydantic psycopg2-binary tenacity

# Verify all CLI tools
orthofinder --version
mafft --version
hyphy --version
busco --version

# Docker (for Nextflow cloud execution)
docker build -t bioresillient:latest ./docker/
```

---

## Step 0: Automated Species Discovery *(Optional — run before Step 1)*

### What it does
Automatically identifies candidate species for a given trait (e.g. cancer resistance, longevity, hypoxia tolerance) by mining three data sources: the AnAge longevity database for statistical outliers, PubMed for published evidence of the trait, and NCBI Assembly for genome availability. Output is a ranked candidate list ready to review and seed into Step 1.

This step replaces hours of manual literature review. It does not replace human judgement — the output is a shortlist for review, not an automatic registry population.

### Input
A trait keyword: `cancer_resistance`, `longevity`, `hypoxia_tolerance`, or any free-text biological trait

### Output
`data/candidate_species.csv` — ranked table of candidates with evidence scores and genome accession IDs, ready to copy into `config/species_registry.json` after human review

### Data Sources

**AnAge Database** — lifespan + body mass for ~4,000 species
```
URL: https://genomics.senescence.info/species/dataset.zip
Auth: None — free flat file download
Key fields: species, max_longevity, body_mass_g, kingdom, order
Use: Compute longevity quotient (LQ) = actual_lifespan / expected_lifespan_from_body_mass
     Species with LQ > 2.0 are statistical longevity outliers — strong cancer resistance candidates
```

**NCBI PubMed API** — literature evidence per candidate species
```
Base: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
Endpoint: esearch.fcgi?db=pubmed&term="{species_name}"AND("{trait_keyword}"OR"tumor resistance"OR"neoplasia")
Auth: API key (optional, raises limit to 10 req/sec)
Use: Count papers + retrieve titles as evidence summary
```

**NCBI Assembly API** — genome availability check
```
Base: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
Endpoint: esearch → elink → esummary on db=assembly
Use: Check whether a sequenced genome exists; retrieve assembly accession ID
     Species without a genome assembly cannot be used in the pipeline
```

### Longevity Quotient Calculation
```python
import numpy as np

def longevity_quotient(max_longevity_years, body_mass_g):
    # Expected lifespan from Speakman's body mass regression (mammals)
    # log(lifespan) = 0.209 * log(body_mass) + 0.959
    expected = 10 ** (0.209 * np.log10(body_mass_g) + 0.959)
    return max_longevity_years / expected

# LQ > 2.0: strong candidate | LQ > 4.0: exceptional (naked mole rat = ~9.0)
```

### Candidate Scoring Formula
```python
# Composite score combining statistical outlier status + evidence strength
candidate_score = (
    min(longevity_quotient / 10, 1.0) * 0.40 +   # LQ (capped at 10)
    min(pubmed_count / 20, 1.0)       * 0.35 +   # Literature evidence (capped at 20 papers)
    genome_available                   * 0.25     # Binary: genome exists in NCBI
)
```

### Output Format
```csv
rank,species_name,common_name,max_longevity_yrs,body_mass_g,longevity_quotient,pubmed_count,genome_available,assembly_accession,candidate_score,evidence_summary
1,Heterocephalus glaber,Naked mole rat,32,35,9.1,847,True,GCF_000247695.1,0.97,"Cancer resistance: HMM-HA contact inhibition..."
2,Spalax ehrenbergi,Blind mole rat,21,180,4.2,312,True,GCF_000622305.1,0.81,"Interferon-mediated concerted cell death..."
3,Loxodonta africana,African elephant,70,3600000,2.8,1203,True,GCF_000001905.1,0.79,"TP53 gene duplications..."
```

### Database Schema
```sql
CREATE TABLE species_candidates (
    id SERIAL PRIMARY KEY,
    species_name VARCHAR(255),
    common_name VARCHAR(255),
    ncbi_taxon_id INTEGER,
    trait VARCHAR(100),
    max_longevity_years FLOAT,
    body_mass_g FLOAT,
    longevity_quotient FLOAT,
    pubmed_count INTEGER,
    pubmed_titles JSONB,               -- top 5 paper titles as evidence
    genome_available BOOLEAN,
    assembly_accession VARCHAR(100),
    candidate_score FLOAT,
    reviewed BOOLEAN DEFAULT FALSE,    -- human sets this to True before Step 1 import
    approved BOOLEAN DEFAULT FALSE,    -- human approves for registry
    created_at TIMESTAMP DEFAULT NOW()
);
```

### Claude Code Task
```
Build pipeline/step0_species_discovery.py:
1. Download AnAge flat file from genomics.senescence.info/species/dataset.zip
   Parse CSV → load all mammal records into species_candidates table
2. Compute longevity_quotient for each species using Speakman regression
3. Filter: LQ > 2.0 OR body_mass_g > 5000 (large body mass = Peto's paradox candidate)
   This typically yields 80–150 candidates from ~4,000 species
4. For each candidate, query NCBI PubMed:
   term = f'"{species_name}" AND ("{trait}" OR "cancer resistance" OR "tumor resistance")'
   Store paper count + top 5 titles in pubmed_titles JSONB
5. For each candidate with pubmed_count > 0, query NCBI Assembly:
   Check genome availability, retrieve assembly accession if exists
6. Compute candidate_score using formula above
7. Export top 30 candidates to data/candidate_species.csv
8. CLI: python step0_species_discovery.py --trait cancer_resistance --top 30
         python step0_species_discovery.py --trait longevity --min-lq 3.0

Phase 1 note: This step is entirely local — no cloud, no heavy compute.
After running, human reviews data/candidate_species.csv and copies approved
species into config/species_registry.json before running Step 1.
```

---

## Step 1: Species-Trait Registry

### What it does
Creates and manages the database of species → extreme trait → evidence mappings. Starting point for every pipeline run.

### Input
Manual curation via `config/species_registry.json` + optional PubMed literature queries

### Output
PostgreSQL table `species_registry` queryable by trait

### Tools
- Python + SQLAlchemy
- Biopython Entrez API (PubMed lookup)
- PostgreSQL

### API
```
NCBI Entrez: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
Auth: API key optional — get at https://www.ncbi.nlm.nih.gov/account/
Rate: 3 req/sec without key, 10 req/sec with key
Biopython: Bio.Entrez.esearch(), Bio.Entrez.efetch()
```

### Database Schema
```sql
CREATE TABLE species_registry (
    id SERIAL PRIMARY KEY,
    ncbi_taxon_id INTEGER UNIQUE NOT NULL,
    species_name VARCHAR(255) NOT NULL,
    common_name VARCHAR(255),
    trait VARCHAR(100) NOT NULL,          -- 'longevity', 'cancer_resistance', 'hypoxia_tolerance'
    trait_evidence_score FLOAT,           -- 0.0 to 1.0
    genome_available BOOLEAN DEFAULT FALSE,
    genome_assembly_id VARCHAR(100),      -- e.g. GCF_000001405.40
    genome_quality_score FLOAT,           -- BUSCO completeness % (filled in Step 3)
    notes TEXT,
    created_at TIMESTAMP DEFAULT NOW()
);
```

### Claude Code Task
```
Build pipeline/step1_registry.py:
1. Reads config/species_registry.json
2. Connects to PostgreSQL via SQLAlchemy, upserts records into species_registry
3. Function lookup_pubmed(species_name, trait) → queries NCBI PubMed via
   Bio.Entrez, returns top 3 paper titles + PMIDs
4. CLI: python step1_registry.py --seed
         python step1_registry.py --lookup "naked mole rat" longevity
```

### Output to Step 2
```json
[
  {
    "ncbi_taxon_id": 10181,
    "species_name": "Heterocephalus glaber",
    "trait": "cancer_resistance",
    "genome_assembly_id": "GCF_000247695.1"
  }
]
```

---

## Step 2: Genomic Data Acquisition

### What it does
Downloads three file types per species from NCBI: genome FASTA, proteome FASTA (for OrthoFinder), and CDS FASTA (required for PAL2NAL codon alignment in Step 5b).

### Input
Species list from Step 1 (NCBI assembly accession IDs)

### Output
```
s3://bioresillient-data/species/{species_name}/
    genome.fa.gz        # Full genome assembly
    proteome.faa.gz     # Protein sequences (for OrthoFinder)
    cds.ffn.gz          # Coding sequences (for PAL2NAL)
    metadata.json
```

### Tools
- Biopython `Bio.Entrez` — metadata and FTP path lookup
- NCBI FTP — file download (preferred over Entrez API for large files)
- `boto3` — S3 upload

### API
```
NCBI Entrez:
  Base: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
  esearch → find assembly from species name
  efetch  → get assembly metadata including FTP path

NCBI FTP:
  Base: ftp.ncbi.nlm.nih.gov/genomes/all/GCF/
  Files per assembly:
    *_genomic.fna.gz              → genome
    *_protein.faa.gz              → proteome
    *_cds_from_genomic.fna.gz     → CDS sequences
  Auth: Anonymous FTP, no credentials required
```

### Claude Code Task
```
Build pipeline/step2_download.py:
1. Reads species with genome_available=False from PostgreSQL
2. For each species:
   a. Bio.Entrez.esearch(db='assembly', term=assembly_id) → get FTP path
   b. Downloads genome.fa.gz, proteome.faa.gz, cds.ffn.gz via FTP
   c. Validates checksums against NCBI md5checksums.txt
   d. Uploads to S3 at s3://{BUCKET}/species/{species_name}/
   e. Sets species_registry.genome_available = True in DB
3. Rate limit: max 10 req/sec (with NCBI API key)
4. Retry failed downloads up to 3 times with exponential backoff
5. CLI: python step2_download.py --species "all"
         python step2_download.py --species "Heterocephalus glaber"
```

---

## Step 3: Preprocessing & Quality Control

### What it does
Masks repetitive regions in genomes (RepeatMasker) and scores proteome completeness against a reference gene set (BUSCO). Marks low-quality genomes for exclusion before ortholog analysis.

### Input
Raw genome + proteome files from Step 2

### Output
```
s3://bioresillient-data/processed/{species_name}/
    genome.masked.fa.gz     # RepeatMasked genome
    proteome.faa.gz         # Unchanged (BUSCO runs on this)
    busco_report.json       # Completeness score
    qc_summary.json         # {busco_complete: 94.2, status: "passed"}
```

### Tools
- **RepeatMasker 4.1** — masks transposable elements and repeats
  - `conda install -c bioconda repeatmasker`
  - CLI: `RepeatMasker -pa 8 -species mammalia genome.fa`
  - nf-core module available
- **BUSCO 5.7** — scores proteome completeness
  - `conda install -c bioconda busco`
  - CLI: `busco -i proteome.faa -m proteins -l mammalia_odb10 -o busco_out`
  - nf-core module available

### QC Thresholds (configurable in `config/`)
```python
QC_THRESHOLDS = {
    "busco_completeness_min": 85.0,   # % complete BUSCOs required
    "busco_duplicated_max": 10.0,     # % duplicated BUSCOs allowed
    "genome_size_min_mb": 500,
}
```

### Claude Code Task
```
Build pipeline/step3_qc.py:
1. Pulls unprocessed species from S3
2. Per species:
   a. Run RepeatMasker via subprocess: RepeatMasker -pa {threads} -species mammalia {genome}
      Parse .tbl output for masking stats
   b. Run BUSCO via subprocess: busco -i {proteome} -m proteins -l mammalia_odb10 -o {out}
      Parse short_summary.json for completeness %
   c. Apply QC thresholds → set pass/fail
   d. Write qc_summary.json to S3
   e. Update PostgreSQL: genome_quality_score, qc_status='passed'/'failed'
3. Generate aggregate QC markdown table for all species
4. Build nextflow/modules/busco.nf and nextflow/modules/repeatmasker.nf
```

---

## Step 4: Ortholog Identification

### What it does
Runs OrthoFinder across all species proteomes to identify which genes are the same across species. Produces orthogroups (gene families) and 1:1 ortholog mappings to human genes. The species tree output feeds into Step 6.

### Input
All QC-passed proteome FASTA files from Step 3 (one `.faa` per species)

### Output
```
s3://bioresillient-data/orthogroups/
    OG0000001.faa       # Multi-species protein sequences per gene family
    OG0000002.faa
    ...
    species_tree.nwk    # Newick species tree (used in Step 6)
```

### Tools
- **OrthoFinder 3.0** — CLI only, no REST API
  - `conda install -c bioconda orthofinder`
  - CLI: `orthofinder -f proteomes_dir/ -t 32 -a 8 -o output_dir/`
  - Single job requiring 32 cores + 64 GB RAM; 2–8 hours for 30 species

### Supplementary API (optional: cross-validate results)
```
OrthoDB REST API — pre-built ortholog reference database
  Base: https://www.orthodb.org/
  Endpoints: /search?query={gene_name}&level={taxon_id}
             /group?id={OG_id}
  Auth: None
  Python wrapper: pip install orthodb_py
```

### Claude Code Task
```
Build pipeline/step4_orthologs.py:
1. Pull all QC-passed proteomes from S3 → stage to /tmp/proteomes/ (one .faa per species)
2. Run OrthoFinder:
   subprocess.run(["orthofinder", "-f", "/tmp/proteomes/",
                   "-t", "32", "-a", "8", "-o", "/tmp/orthofinder_out/"])
3. Parse output:
   a. Orthogroups.tsv → load into PostgreSQL table 'orthogroups'
   b. All Homo_sapiens__v__{species}.tsv → load into 'human_orthologs'
4. Copy per-orthogroup FASTAs to s3://bioresillient-data/orthogroups/{OG_id}.faa
5. Copy SpeciesTree_rooted.txt to S3 (needed for Step 6)
6. CLI: python step4_orthologs.py --run / --parse-only

Build nextflow/modules/orthofinder.nf
  label: 'high_cpu', cpus: 32, memory: '64 GB'
```

### Database Schemas
```sql
CREATE TABLE orthogroups (
    og_id VARCHAR(20) PRIMARY KEY,     -- e.g. OG0000042
    species_count INTEGER,
    gene_count INTEGER,
    has_human_ortholog BOOLEAN
);

CREATE TABLE human_orthologs (
    id SERIAL PRIMARY KEY,
    og_id VARCHAR(20) REFERENCES orthogroups(og_id),
    human_ensembl_id VARCHAR(20),      -- Populated in Step 7a
    human_gene_id VARCHAR(100),        -- OrthoFinder internal ID
    species VARCHAR(100),
    species_gene_id VARCHAR(100),
    ortholog_type VARCHAR(20)          -- '1:1', '1:many', 'many:1'
);
```

---

## Step 5a: Multiple Sequence Alignment

### What it does
Aligns protein sequences for each orthogroup across all species using MAFFT. Required prerequisite for codon alignment in Step 5b.

### Input
Per-orthogroup FASTA files from Step 4

### Output
`s3://bioresillient-data/aligned/OG0000001.aligned.faa` — one per orthogroup (~20,000 files)

### Tools
- **MAFFT 7.526**
  - `conda install -c bioconda mafft`
  - CLI: `mafft --auto --thread 4 OG0000001.faa > OG0000001.aligned.faa`
  - nf-core module: `nf-core/modules/mafft`

### Claude Code Task
```
Build pipeline/step5a_alignment.py:
1. List all OG*.faa in S3 orthogroups bucket
2. Submit each as a separate AWS Batch job (Nextflow scatter):
   mafft --auto --thread 4 {input}.faa > {output}.aligned.faa
3. Monitor completion, retry failed jobs up to 3 times
4. Upload to s3://bioresillient-data/aligned/
5. Track job status in PostgreSQL

Build nextflow/modules/mafft.nf:
  process MAFFT {
    label 'low_cpu'
    cpus 4; memory '8 GB'
    input: path fasta
    output: path "*.aligned.faa"
    script: "mafft --auto --thread ${task.cpus} ${fasta} > ${fasta.baseName}.aligned.faa"
  }
```

---

## Step 5b: Codon Alignment

### What it does
Converts protein alignments (Step 5a) into codon-level DNA alignments. **Non-negotiable.** HyPhy's dN/dS models (BUSTED, aBSREL) operate on codon-level data — they require knowing which DNA changes are synonymous vs non-synonymous, which cannot be computed from protein alignments alone. PAL2NAL maps protein alignment positions back to underlying codon triplets using the CDS sequences from Step 2.

### Input
- Protein alignment from Step 5a: `OG0000001.aligned.faa`
- CDS sequences from Step 2: `{species}_cds.ffn.gz`

### Output
`s3://bioresillient-data/codon_aligned/OG0000001.codon.fna`

### Tools
- **PAL2NAL v14** — primary tool
  - `conda install -c bioconda pal2nal`
  - CLI: `pal2nal.pl aligned.faa cds.fna -output fasta > codon_alignment.fna`
- **MACSE v2** — fallback for frameshifts/pseudogenes (PAL2NAL fails on these)
  - `conda install -c bioconda macse`
  - CLI: `java -jar macse.jar -prog alignSequences -seq cds.fna`

### Claude Code Task
```
Build pipeline/step5b_codon_alignment.py:
1. For each orthogroup collect:
   a. Protein alignment OG*.aligned.faa from S3
   b. CDS sequences for each species in the OG, extracted from per-species cds.ffn.gz
2. Run PAL2NAL:
   subprocess.run(["pal2nal.pl", aligned_faa, cds_fna, "-output", "fasta"],
                  capture_output=True)
3. If PAL2NAL exits non-zero (frameshift/stop codon), fall back to MACSE:
   subprocess.run(["java", "-jar", "macse.jar", "-prog", "alignSequences", "-seq", cds_fna])
4. Validate output: check for in-frame stop codons, minimum alignment length
5. Upload to s3://bioresillient-data/codon_aligned/
6. Log PAL2NAL vs MACSE usage and failure rates to PostgreSQL

Build nextflow/modules/pal2nal.nf with MACSE fallback logic
```

---

## Step 5c: Alignment Trimming

### What it does
Removes poorly aligned, gappy regions from codon alignments. Filters out orthogroups too short for reliable dN/dS calculation.

### Input
Codon alignments from Step 5b

### Output
`s3://bioresillient-data/trimmed/OG0000001.trimmed.fna`

### Tools
- **trimAl 1.4**
  - `conda install -c bioconda trimal`
  - CLI: `trimal -in codon.fna -out trimmed.fna -automated1`
  - nf-core module: `nf-core/modules/trimal`

### Claude Code Task
```
Build pipeline/step5c_trim.py:
- Run trimal -automated1 on all codon alignments from S3
- Filter: skip orthogroups where post-trim alignment < 100 codons (log and mark skipped)
- Upload trimmed files to s3://bioresillient-data/trimmed/
```

---

## Step 6: Evolutionary Selection Analysis

### What it does
The scientific core of the platform. For each of ~20,000 genes, tests whether it shows statistically significant positive selection in the trait-bearing species lineages (foreground) vs background species. Produces the evolutionary evidence score per gene.

### Input
- Trimmed codon alignments from Step 5c
- Species tree from Step 4: `species_tree.nwk`
- Trait-bearing species labels from Step 1 (foreground vs background designations)

### Output
Per-gene HyPhy JSON: `OG0000001.hyphy.json`
```json
{
  "BUSTED": {
    "test results": { "p-value": 0.0023, "LRT": 12.4 }
  },
  "branch attributes": {
    "Heterocephalus_glaber": {"dN/dS": 3.2},
    "Homo_sapiens": {"dN/dS": 0.3}
  }
}
```

### Tools

**IQ-TREE2 2.3** — builds high-quality species tree (runs once only)
- `conda install -c bioconda iqtree`
- CLI: `iqtree2 -s concatenated.fna -m TEST -bb 1000 -T AUTO`
- Note: Use OrthoFinder's SpeciesTree_rooted.txt directly if sufficient quality

**HyPhy 2.5** — per-gene selection analysis (~20,000 independent CPU jobs)
- `conda install -c bioconda hyphy`
- Tests to run:
  - **BUSTED** — gene-wide: "does this gene show any positive selection?"
  - **aBSREL** — branch-specific: "which branches show selection?"
  - **RELAX** — "is selection relaxed/intensified in trait-bearing lineages?"
- CLI: `hyphy busted --alignment trimmed.fna --tree labelled.nwk`
- Parallelism is at the job level — 20,000 independent CPU jobs on cloud batch

**FDR Correction** — runs once after ALL gene jobs complete, not per-batch
```python
from statsmodels.stats.multitest import multipletests

pvalues = [result['p_value'] for result in all_gene_results]
reject, pvals_corrected, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')
# fdr_bh = Benjamini-Hochberg
```

### Claude Code Task
```
Build pipeline/step6_selection.py:

Part A — Species tree (run once):
1. Concatenate subset of high-quality single-copy OG alignments
2. Run IQ-TREE2 (or use OrthoFinder's tree if quality sufficient)

Part B — Per-gene selection (parallelised across ~20,000 jobs):
1. For each trimmed codon alignment:
   a. Label foreground branches in species tree from trait tags in registry
      (HyPhy label format: {species_name}{FOREGROUND})
   b. Run BUSTED: hyphy busted --alignment {trimmed.fna} --tree {labelled.nwk}
   c. Run aBSREL: hyphy absrel --alignment {trimmed.fna} --tree {labelled.nwk}
   d. Parse JSON → extract p-value, LRT, per-branch dN/dS
   e. Store in PostgreSQL 'selection_results'

Part C — FDR correction (after ALL genes complete):
1. Load all p-values from PostgreSQL
2. Apply Benjamini-Hochberg via statsmodels
3. Mark significant genes (FDR < 0.05) in database

Build nextflow/modules/hyphy.nf:
  Separate process definitions for BUSTED and aBSREL
  label: 'medium_cpu', cpus: 8, memory: '16 GB'
  Run on spot instances via AWS Batch
```

### Database Schema
```sql
CREATE TABLE selection_results (
    id SERIAL PRIMARY KEY,
    og_id VARCHAR(20) REFERENCES orthogroups(og_id),
    test_type VARCHAR(20),             -- 'BUSTED', 'aBSREL', 'RELAX'
    trait VARCHAR(100),
    p_value FLOAT,
    p_value_corrected FLOAT,           -- BH-FDR corrected
    lrt_statistic FLOAT,
    is_significant BOOLEAN,
    foreground_dnds FLOAT,
    background_dnds FLOAT,
    convergent_species_count INTEGER,
    raw_json JSONB
);
```

---

## Step 7a: Gene ID Normalisation

### What it does
OrthoFinder outputs internal sequence IDs. All downstream annotation APIs require **Ensembl Gene IDs** (format: `ENSG00000141510`). This step maps OrthoFinder human gene IDs → Ensembl IDs via the Ensembl REST API. Without this, the annotation pipeline cannot function.

### Input
`human_orthologs.human_gene_id` from PostgreSQL (Step 4 output)

### Output
`human_orthologs.human_ensembl_id` populated in PostgreSQL

### API
```
Ensembl REST API — free, no auth
Base: https://rest.ensembl.org
Endpoints:
  POST /lookup/symbol/homo_sapiens
    body: {"symbols": ["TP53", "BRCA1", ...]}
    → batch lookup, up to 1000 symbols per request
  GET  /lookup/symbol/homo_sapiens/{symbol}
    → single gene lookup
Headers: Content-Type: application/json
Rate limit: 15 req/sec
```

### Claude Code Task
```
Build pipeline/step7a_normalise_ids.py:
1. Query PostgreSQL for all distinct human_gene_id values in human_orthologs
2. Batch into groups of 1000 gene symbols
3. POST https://rest.ensembl.org/lookup/symbol/homo_sapiens per batch
4. Parse response → map gene_symbol → ensembl_id
5. Handle genes with no Ensembl mapping (log, flag as unmappable)
6. Bulk update human_orthologs.human_ensembl_id in PostgreSQL
7. Use tenacity for retry logic (max 3 retries, exponential backoff)
```

---

## Step 7b: Signal Integration & Annotation

### What it does
For each significantly selected human gene (FDR < 0.05), queries six external databases to build the full annotation: disease relevance, druggability, existing compounds, tissue expression, protein function, and pathway membership.

### Input
Significant genes from `selection_results` with Ensembl IDs from Step 7a

### Output
PostgreSQL table `gene_annotations`

### APIs — All Verified, All Free

#### 1. OpenTargets Platform API
```
Endpoint: https://api.platform.opentargets.org/api/v4/graphql
Auth: None | Format: GraphQL
Returns: disease associations, tractability scores, approved drugs
Query:
  target(ensemblId: "ENSG00000141510") {
    approvedSymbol
    tractability { modality id value }
    associatedDiseases {
      rows { disease { id name } score }
    }
  }
```

#### 2. OpenTargets Genetics Portal
```
Endpoint: https://api.genetics.opentargets.org/graphql
Auth: None | Format: GraphQL
Returns: GWAS evidence, L2G scores — separate from Platform API, query independently
```

#### 3. ChEMBL REST API
```
Base: https://www.ebi.ac.uk/chembl/api/data/
Auth: None | Rate: ~10 req/sec
Endpoints:
  /target/search.json?q={gene_symbol}
  /mechanism.json?target_chembl_id={id}
  /drug_indication.json?molecule_chembl_id={id}
Returns: existing drugs, mechanisms, clinical trial phases
```

#### 4. Human Protein Atlas REST
```
Base: https://www.proteinatlas.org/api/
Auth: None
Endpoint: GET https://www.proteinatlas.org/{gene_symbol}.json
Returns: tissue-specific expression levels, subcellular location
Rate: add 0.2s delay between requests
```

#### 5. UniProt REST
```
Base: https://rest.uniprot.org/
Auth: None
Endpoint: GET /uniprotkb/search?query=gene:{symbol}+organism_id:9606&format=json
Returns: protein function, GO terms, pathway involvement
```

#### 6. Reactome REST
```
Base: https://reactome.org/ContentService/
Auth: None
Endpoint: GET /data/pathways/low/entity/{uniprot_id}?species=9606
Returns: pathway membership
```

### Database Schema
```sql
CREATE TABLE gene_annotations (
    ensembl_id VARCHAR(20) PRIMARY KEY,
    gene_symbol VARCHAR(50),
    gene_name TEXT,
    biotype VARCHAR(50),

    -- OpenTargets Platform
    ot_disease_associations JSONB,    -- [{disease_id, disease_name, score}]
    ot_tractability JSONB,            -- {small_molecule: bool, antibody: bool}
    ot_top_disease VARCHAR(200),
    ot_max_disease_score FLOAT,

    -- GWAS
    gwas_evidence JSONB,

    -- ChEMBL
    chembl_target_id VARCHAR(20),
    existing_drugs JSONB,             -- [{drug_name, phase, indication}]
    drug_count INTEGER,
    is_drugged BOOLEAN,

    -- Protein Atlas
    tissue_expression JSONB,
    subcellular_location TEXT[],

    -- UniProt
    uniprot_id VARCHAR(15),
    protein_function TEXT,
    go_terms JSONB,

    -- Reactome
    pathways JSONB,                   -- [{pathway_id, pathway_name}]

    updated_at TIMESTAMP DEFAULT NOW()
);
```

### Claude Code Task
```
Build pipeline/step7b_annotate.py:

OPTIMISATION 1 — Async parallel API calls:
Use asyncio + aiohttp so all 6 APIs are queried concurrently per gene,
and multiple genes are processed in parallel. Do NOT use sequential requests.

import asyncio, aiohttp

# Per-API semaphores enforce individual rate limits
SEMAPHORES = {
    "opentargets":    asyncio.Semaphore(10),
    "chembl":         asyncio.Semaphore(8),
    "protein_atlas":  asyncio.Semaphore(5),
    "uniprot":        asyncio.Semaphore(10),
    "reactome":       asyncio.Semaphore(10),
}

async def annotate_gene(session, ensembl_id):
    tasks = [
        fetch_opentargets_platform(session, ensembl_id),
        fetch_opentargets_genetics(session, ensembl_id),
        fetch_chembl(session, ensembl_id),
        fetch_protein_atlas(session, ensembl_id),
        fetch_uniprot(session, ensembl_id),
        fetch_reactome(session, ensembl_id),
    ]
    results = await asyncio.gather(*tasks, return_exceptions=True)
    return merge_annotations(ensembl_id, results)

async def main(gene_list, concurrency=50):
    semaphore = asyncio.Semaphore(concurrency)  # max 50 genes in parallel
    async with aiohttp.ClientSession() as session:
        tasks = [annotate_gene(session, g) for g in gene_list]
        return await asyncio.gather(*tasks)

OPTIMISATION 2 — Cache responses, skip already-annotated genes:
Before querying any API, check if gene_annotations record exists AND
updated_at > NOW() - INTERVAL '30 days'. If fresh cache exists, skip all
API calls for that gene. Only re-query if record is missing or stale.

if not annotation_is_fresh(ensembl_id, max_age_days=30):
    result = await annotate_gene(session, ensembl_id)
    upsert_annotation(result)
else:
    log.info(f"Cache hit for {ensembl_id}, skipping API calls")

Full task:
1. Query PostgreSQL for all significant genes with ensembl_ids
2. Split into: needs_annotation (missing/stale) vs cache_hit (fresh)
3. Run async annotation for needs_annotation list
4. Merge all API responses into gene_annotations record (upsert)
5. Retry logic via tenacity on each individual API coroutine
6. CLI: python step7b_annotate.py --gene ENSG00000141510
         python step7b_annotate.py --all
         python step7b_annotate.py --all --force-refresh   # ignore cache
```

---

## Step 8: Target Scoring & Ranking

### What it does
Combines evolutionary evidence (Step 6) with annotation richness (Step 7) into a weighted composite score. Produces the final ranked list of drug target candidates.

### Input
- `selection_results` (evolutionary scores)
- `gene_annotations` (druggability, disease relevance)
- `config/scoring_weights.json`

### Output
- `ranked_targets` table in PostgreSQL
- Per-target PDF report (ReportLab)
- Excel export of full ranked list (openpyxl)

### Scoring Formula
```python
# config/scoring_weights.json
{
  "evolutionary_signal": 0.35,    # FDR-corrected p-value significance
  "convergence_bonus": 0.15,      # signal in 3+ independent species
  "disease_relevance": 0.20,      # OpenTargets max disease association score
  "druggability": 0.15,           # OpenTargets tractability score
  "novelty": 0.10,                # penalty if existing approved drugs
  "gwas_support": 0.05            # bonus for human genetic evidence
}

def compute_score(gene, weights):
    evol    = 1 - gene.p_value_corrected
    conv    = min(gene.convergent_species / 5, 1)
    disease = gene.ot_max_disease_score
    drug    = gene.tractability_score
    novelty = 0 if gene.is_drugged else 1
    gwas    = 1 if gene.gwas_evidence else 0

    return (evol    * weights['evolutionary_signal'] +
            conv    * weights['convergence_bonus'] +
            disease * weights['disease_relevance'] +
            drug    * weights['druggability'] +
            novelty * weights['novelty'] +
            gwas    * weights['gwas_support']) * 100
```

### Claude Code Task
```
Build pipeline/step8_scoring.py:
1. Join selection_results + gene_annotations on ensembl_id
2. Compute composite score using formula above
3. Support multiple traits (longevity, cancer_resistance) — separate ranked lists per trait
4. Write to ranked_targets table
5. Excel export (openpyxl):
   - Sheet 1: Full ranked list with all score components
   - Sheet 2: Top 20 targets with evidence summary
6. PDF per top-10 target (ReportLab):
   - Gene name, rank, score breakdown, disease associations,
     evolutionary evidence summary, existing compounds
7. CLI: python step8_scoring.py --trait longevity --top 20 --export pdf,excel
```

### Database Schema
```sql
CREATE TABLE ranked_targets (
    id SERIAL PRIMARY KEY,
    ensembl_id VARCHAR(20),
    gene_symbol VARCHAR(50),
    trait VARCHAR(100),
    rank INTEGER,
    composite_score FLOAT,
    evolutionary_signal_score FLOAT,
    convergence_score FLOAT,
    disease_relevance_score FLOAT,
    druggability_score FLOAT,
    novelty_score FLOAT,
    gwas_score FLOAT,
    top_disease VARCHAR(200),
    existing_drug_count INTEGER,
    is_novel BOOLEAN,
    analysis_version VARCHAR(20),
    created_at TIMESTAMP DEFAULT NOW()
);
```

---

## Step 9: Nextflow Pipeline Orchestration

### What it does
Wraps all steps into a single reproducible Nextflow pipeline. One command triggers the full analysis, handles retries automatically, and runs on cloud (AWS Batch) with monitoring via Seqera Platform.

### Tools
- **Nextflow DSL2** — `conda install -c bioconda nextflow`
- **Seqera Platform** — monitoring: `nextflow run main.nf -with-tower`
- **AWS Batch** — cloud execution, configured via `nextflow.config`
- **Docker** — required for cloud execution

### Master Pipeline (main.nf)
```groovy
nextflow.enable.dsl=2

include { DOWNLOAD }      from './nextflow/modules/download'
include { QC }            from './nextflow/modules/qc'
include { ORTHOFINDER }   from './nextflow/modules/orthofinder'
include { MAFFT }         from './nextflow/modules/mafft'
include { PAL2NAL }       from './nextflow/modules/pal2nal'
include { TRIMAL }        from './nextflow/modules/trimal'
include { HYPHY_BUSTED }  from './nextflow/modules/hyphy'
include { ANNOTATE }      from './nextflow/modules/annotate'
include { SCORE }         from './nextflow/modules/score'

workflow {
    species_ch = Channel.fromPath('config/species_registry.json').splitJson()

    genomes      = DOWNLOAD(species_ch)
    qc_passed    = QC(genomes).filter { it.qc_status == 'passed' }
    orthogroups  = ORTHOFINDER(qc_passed.collect())
    aligned      = MAFFT(orthogroups.flatten())
    codon        = PAL2NAL(aligned, genomes)
    trimmed      = TRIMAL(codon)
    selection    = HYPHY_BUSTED(trimmed, orthogroups.species_tree)
    annotated    = ANNOTATE(selection.significant_genes)
    SCORE(annotated)
}
```

### nextflow.config
```groovy
profiles {
    local {
        process.executor = 'local'
        docker.enabled = true
    }
    aws {
        process.executor = 'awsbatch'
        process.container = '123456789.dkr.ecr.us-east-1.amazonaws.com/bioresillient:latest'
        aws.region = 'us-east-1'
        workDir = 's3://bioresillient-work/nextflow-tmp'

        process {
            withLabel: 'low_cpu'    { cpus = 4;  memory = '8 GB';  queue = 'spot' }
            withLabel: 'medium_cpu' { cpus = 8;  memory = '16 GB'; queue = 'spot' }
            withLabel: 'high_cpu'   { cpus = 32; memory = '64 GB'; queue = 'on-demand' }
        }
    }
    gcp {
        process.executor = 'google-lifesciences'
        google.region = 'us-central1'
        google.project = 'bioresillient'
        workDir = 'gs://bioresillient-work/nextflow-tmp'
    }
}
```

### Claude Code Task
```
Build the complete Nextflow pipeline:
1. main.nf — master workflow wiring all modules
2. nextflow.config — local, aws, gcp profiles
3. nextflow/modules/ — one .nf file per tool
4. docker/Dockerfile — all conda tools in one image
5. Test locally: nextflow run main.nf -profile local --input test_data/
6. Test on AWS: nextflow run main.nf -profile aws -with-tower
```

---

## Step 10: Validation Interface (FastAPI)

### What it does
Internal web app making ranked targets browsable, generating experimental briefs for wet-lab partners, and capturing validation results back into the system.

### Tools
- FastAPI + Uvicorn + PostgreSQL + Jinja2 + ReportLab

### API Routes
```python
GET  /targets                          # List ranked targets (filter by trait, score, novelty)
GET  /targets/{ensembl_id}             # Full target detail
GET  /targets/{ensembl_id}/brief       # Download PDF experimental brief
POST /targets/{ensembl_id}/validation  # Submit wet-lab result
GET  /species                          # Species in registry
GET  /runs                             # Pipeline run history
GET  /runs/{run_id}/status             # Run status
```

### Experimental Brief PDF Structure
```
Target: FOXO3 (ENSG00000118689)
Rank: #1 of 847 | Score: 94.2/100 | Trait: Longevity

Evolutionary Evidence:
  - Significant positive selection (BUSTED p=0.0003, FDR=0.004)
  - Convergent signal in 4 species: Bowhead whale, Naked mole rat,
    Greenland shark, Blanding's turtle
  - Foreground dN/dS: 2.8 vs Background: 0.12

Human Disease Relevance:
  - Top disease: Type 2 Diabetes (OT score: 0.82)
  - Also associated with: Alzheimer's, Cardiovascular disease

Druggability:
  - Tractable by small molecule (OpenTargets tractability: confirmed)
  - 3 existing compounds in ChEMBL (none approved for aging indication)
  - Novel opportunity: no approved drug for longevity indication

Pathways: Insulin signalling, DNA damage response, Autophagy regulation

Suggested Validation Assay:
  - Cell type: Human fibroblasts (longevity) or HCT116 (cancer resistance)
  - Method: CRISPR-KO + RNAseq to confirm transcriptomic phenotype
  - Expected: Accelerated senescence or reduced stress resistance
```

### Claude Code Task
```
Build the FastAPI app in api/:
1. api/main.py — FastAPI app with PostgreSQL via SQLAlchemy
2. api/routes/targets.py — listing, filtering, detail view
3. api/routes/validation.py — wet-lab result submission + history
4. api/templates/brief.html — Jinja2 template for experimental brief
5. api/pdf.py — ReportLab function to convert brief → PDF
6. Run: uvicorn api.main:app --host 0.0.0.0 --port 8000
7. Dockerfile.api for containerised deployment
```

---

## Complete Data Flow

```
[PostgreSQL: species_registry]
        │ species list + assembly IDs
        ▼
[step2_download.py]
  NCBI Entrez API + FTP → genome.fa.gz + proteome.faa.gz + cds.ffn.gz
        │
        ▼
[step3_qc.py]
  RepeatMasker + BUSCO → masked genome, QC report, pass/fail in DB
        │ QC-passed proteomes
        ▼
[step4_orthologs.py]
  OrthoFinder (32 CPU) → Orthogroups in DB + per-OG FASTAs + species_tree.nwk in S3
        │ ~20,000 OG FASTA files
        ▼
[step5a_alignment.py — AWS Batch, 20,000 parallel jobs]
  MAFFT → aligned protein FASTA per OG
        ▼
[step5b_codon_alignment.py — AWS Batch]
  PAL2NAL (+ MACSE fallback) → codon-level FASTA per OG
        ▼
[step5c_trim.py — AWS Batch]
  trimAl → trimmed codon FASTA per OG
        │ ~20,000 trimmed codon alignments + species tree
        ▼
[step6_selection.py — AWS Batch, 20,000 parallel CPU jobs]
  HyPhy BUSTED + aBSREL → per-gene JSON
  FDR correction (statsmodels BH) → selection_results in DB
        │ significant genes
        ▼
[step7a_normalise_ids.py]
  Ensembl REST API → OrthoFinder IDs mapped to Ensembl IDs in DB
        ▼
[step7b_annotate.py]
  OpenTargets Platform + Genetics, ChEMBL, Human Protein Atlas,
  UniProt, Reactome → gene_annotations in DB
        │ annotated genes
        ▼
[step8_scoring.py]
  Weighted formula → ranked_targets in DB
  → PDF reports + Excel export
        ▼
[FastAPI web app]
  → browsable ranked targets
  → downloadable experimental briefs
  → wet-lab validation intake
  → feedback loop into scoring model
```

---

## All External APIs — Quick Reference

| API | Endpoint | Auth | Format | Used In |
|-----|----------|------|--------|---------|
| NCBI Entrez | `eutils.ncbi.nlm.nih.gov/entrez/` | API key (optional) | XML/JSON | Steps 1, 2 |
| NCBI FTP | `ftp.ncbi.nlm.nih.gov/genomes/all/` | Anonymous | FASTA/gz | Step 2 |
| Ensembl REST | `rest.ensembl.org` | None | JSON | Step 7a |
| OrthoDB REST | `orthodb.org` | None | JSON/FASTA | Step 4 (validation) |
| OpenTargets Platform | `api.platform.opentargets.org/api/v4/graphql` | None | GraphQL | Step 7b |
| OpenTargets Genetics | `api.genetics.opentargets.org/graphql` | None | GraphQL | Step 7b |
| ChEMBL REST | `ebi.ac.uk/chembl/api/data/` | None | JSON | Step 7b |
| Human Protein Atlas | `proteinatlas.org/api/` | None | JSON | Step 7b |
| UniProt REST | `rest.uniprot.org` | None | JSON | Step 7b |
| Reactome REST | `reactome.org/ContentService/` | None | JSON | Step 7b |

---

## All CLI Tools — Installation Reference

| Tool | Version | Conda channel | Used In | nf-core module |
|------|---------|---------------|---------|----------------|
| OrthoFinder | 3.0 | bioconda | Step 4 | Custom |
| MAFFT | 7.526 | bioconda | Step 5a | Yes |
| PAL2NAL | v14 | bioconda | Step 5b | Custom |
| MACSE | v2 | bioconda | Step 5b fallback | Custom |
| trimAl | 1.4 | bioconda | Step 5c | Yes |
| IQ-TREE2 | 2.3 | bioconda | Step 6 | Yes |
| HyPhy | 2.5 | bioconda | Step 6 | Custom |
| RepeatMasker | 4.1 | bioconda | Step 3 | Yes |
| BUSCO | 5.7 | bioconda | Step 3 | Yes |

---

## Testing-First Build Strategy

This plan is written for production, but build it in two phases. Tell Claude Code which phase you're in at the start of each session.

### Phase 1: Test (Current) — 1–3 species, validate science, zero cloud

**Infrastructure changes:**
- Skip AWS Batch, Nextflow, and Docker entirely — use local `subprocess` calls only
- Replace S3 paths with local filesystem paths (`data/` folder)
- Use SQLite instead of PostgreSQL (`sqlite:///bioresillient_test.db`)
- Skip RepeatMasker (Step 3) — run BUSCO only, sufficient for test QC

**Step 6 — use Datamonkey instead of local HyPhy:**
```python
# In step6_selection.py, add --backend flag
# backend='datamonkey' → POST to api.datamonkey.org/msa/busted
# backend='local'      → subprocess hyphy (production)

def run_hyphy_datamonkey(alignment_path, tree_path):
    import requests, time
    payload = {
        "msa": open(alignment_path).read(),
        "tree": open(tree_path).read()
    }
    r = requests.post("http://api.datamonkey.org/msa/busted", json=payload)
    job_id = r.json()["id"]
    while True:
        status = requests.get(f"http://api.datamonkey.org/job/{job_id}").json()
        if status["status"] == "completed":
            return status["result"]
        time.sleep(10)
```

**Scoring weights:** Use equal weights (0.17 each) for now — the production weights in `scoring_weights.json` are assumptions that need wet-lab feedback to calibrate properly.

**Validation checklist for Phase 1:**
- OrthoFinder finds sensible orthogroups across your 3 species
- PAL2NAL success rate > 80% (failure rate tells you CDS download quality)
- HyPhy returns valid JSON with p-values (confirms foreground labelling is correct)
- At least 1–2 genes show significant selection (sanity check)
- PDF brief generates cleanly for top-ranked target

**Prompt to give Claude Code at the start of every Phase 1 session:**
```
We are in Phase 1 (testing). Use SQLite not PostgreSQL, local filesystem
not S3, direct subprocess calls not Nextflow or AWS Batch. Step 6 should
call the Datamonkey REST API (api.datamonkey.org) not local HyPhy.
Build each step as a standalone script with a CLI entry point.
```

### Phase 2: Production (Post-Funding)

Swap in: PostgreSQL → AWS RDS, local paths → S3, subprocess → Nextflow + AWS Batch, Datamonkey → local HyPhy. The core logic in each `pipeline/step*.py` stays identical — only the infrastructure layer changes. No science rewrite needed.

---

## Claude Code — Suggested Build Order

1. `db/models.py` — all SQLAlchemy models first (every other module imports from here)
2. `pipeline/step0_species_discovery.py` *(optional — run before seeding registry)*
3. `pipeline/step1_registry.py` + seed data
4. `pipeline/step2_download.py` — test with 2 species
5. `pipeline/step3_qc.py`
6. `pipeline/step4_orthologs.py` — test with 3 species
7. `pipeline/step5a_alignment.py`
8. `pipeline/step5b_codon_alignment.py` ← do not skip
9. `pipeline/step5c_trim.py`
10. `pipeline/step6_selection.py`
11. `pipeline/step7a_normalise_ids.py` ← do not skip
12. `pipeline/step7b_annotate.py` ← implement with async + cache from the start
13. `pipeline/step8_scoring.py`
14. `nextflow/` modules
15. `api/` — FastAPI last

---

## Key Rules for Claude Code Sessions

- **Test with small dataset first.** Use 3 species (human, naked mole rat, bowhead whale) and 100 orthogroups before scaling.
- **Never hardcode paths.** All paths via `config/` or environment variables.
- **Every subprocess call must capture stdout, stderr, and exit code.** Log failures with enough context to diagnose without re-running.
- **Each step must be re-entrant.** Re-running a failed step should skip completed genes (check PostgreSQL before submitting any job).
- **All S3 uploads must verify.** Read back file size after upload and compare to local.
- **FDR correction runs only after ALL genes are complete.** Never per-batch.
- **Steps 5b and 7a are non-negotiable.** Skipping codon alignment invalidates HyPhy; skipping ID normalisation breaks all annotation APIs.

---

*Version 5.0 — February 2026*
*All 10 external APIs verified — no auth barriers*
*All CLI tools confirmed subprocess-invocable*
*Step 0 added: automated species discovery via AnAge + PubMed + NCBI Assembly*
*Step 7b optimised: async parallel API calls + 30-day response caching*
