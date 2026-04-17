# Step 8 — Functional Evidence Gathering

**Pipeline phase:** Functional annotation — Layer 4  
**Substeps:** 8a (Open Targets), 8b (GTEx + DepMap)  
**Run timestamps:** 8a: 2026-04-06 17:24:00 | 8b: 2026-04-06 17:24:05 UTC  
**Status:** ✅ PASS

---

## What This Step Does

Steps 1–7 are purely evolutionary — they identify genes under adaptive selection and convergent evolution across cancer-resistant species. Step 8 adds the orthogonal layer of **human functional evidence**: do these genes actually connect to cancer biology in human cells, tissues, and patients?

Three independent databases are queried:
1. **Open Targets Platform** — human genetic and clinical evidence for gene–cancer associations
2. **GTEx v10** — human tissue expression profiling
3. **DepMap 24Q4** — CRISPR-based cancer cell line essentiality screens

---

## Step 8a — Open Targets Platform (Disease Associations)

### What is Open Targets?

The **Open Targets Platform** (www.opentargets.org) is a public–private partnership that integrates evidence from multiple sources to score the association between every human protein-coding gene and every disease:

| Evidence type | Examples |
|---|---|
| Genetic associations | GWAS hits, somatic mutations in cancer |
| Transcriptomic evidence | Differentially expressed in disease tissue |
| Literature mining | NLP extraction from PubMed abstracts |
| Animal models | Knockout phenotypes matching disease |
| Clinical trials | Gene product is a drug target in trial |
| Pathway evidence | Gene in pathway associated with disease |

Association scores range 0–1, with 1 = maximum evidence for gene–disease link.

### Query

For each of the 12,795 human genes (identified by Ensembl ID):
```graphql
query {
  target(ensemblId: "ENSG...") {
    associatedDiseases(page: {index: 0, size: 200}) {
      rows {
        disease { id name }
        score
      }
    }
  }
}
```

Filtered to oncology-relevant terms (cancer, neoplasm, carcinoma, tumour, sarcoma, lymphoma).

### Results

| Metric | Value |
|---|---|
| Genes with Open Targets evidence | **5,376** |
| Expression result rows (OT) | 5,376 |
| Score range | 0–1 (association confidence) |

The 5,376 genes with OT evidence represent those with documented connections to cancer in the literature, genetic studies, or clinical trials. The remaining ~7,400 genes have no current OT cancer association — but this does not mean they are unimportant; they may represent novel biology that evolutionary analysis can reveal.

---

## Step 8b — GTEx v10 (Tissue Expression)

### What is GTEx?

The **Genotype-Tissue Expression (GTEx) project** measured RNA-seq expression across **54 human tissues and cell types** from ~1,000 post-mortem donors. It provides the definitive reference for tissue-specific and broadly expressed genes in healthy human adults.

### Why tissue expression matters for cancer

A gene involved in cancer resistance needs to be expressed in the relevant tissues. Broadly expressed genes are more likely to reflect fundamental cellular processes. GTEx allows us to score:
- **Expression breadth** — is the gene expressed in many tissue types? (housekeeping vs. tissue-specific)
- **Expression level** — what is the median TPM in each tissue?

### API

**GTEx v10 REST API** — queries versioned `gencodeId` identifiers (GTEx uses GENCODE gene version IDs, not plain Ensembl IDs). A two-step lookup is required:
1. Map Ensembl ID → GENCODE versioned ID (e.g. `ENSG00000000001.5`)
2. Query GTEx median expression per tissue using the versioned ID

### Score derivation

```
expression_score = median_TPM / (median_TPM + 10)  (normalised Michaelis-Menten)
```

A gene expressed at 10 TPM → score 0.5; at 50 TPM → score 0.83; at 100 TPM → score 0.91.

### Results

| Metric | Value |
|---|---|
| Genes with GTEx data | **5,529** |
| Expression result rows (GTEx) | 5,529 |

---

## Step 8b — DepMap 24Q4 (Cancer Essentiality — Chronos Scores)

### What is DepMap?

The **Cancer Dependency Map (DepMap)** project performed genome-scale CRISPR-Cas9 knockout screens across **~1,100 human cancer cell lines** from ~30 cancer types. For every gene, the Chronos score reports how much each cell line *depends* on that gene for survival and proliferation.

### Chronos score

The **Chronos** dependency score (Dempster et al. 2021, *Nature Genetics*) is derived from guide RNA depletion kinetics in pooled CRISPR screens. Unlike earlier CERES scores, Chronos explicitly models cell division rate and guide RNA efficacy biases.

| Chronos score | Interpretation |
|---|---|
| 0 | Gene knockout has no effect on cancer cell growth |
| -0.5 | Moderate growth disadvantage — likely essential in some contexts |
| -1.0 | Strong growth defect — broadly essential across many cell lines |
| < -1.0 | Pan-essential gene (ribosomal proteins, splicing factors, etc.) |

**Dataset:** `CRISPRGeneEffect.csv` from DepMap 24Q4 release (Figshare, ~1,100 cell lines × ~18,000 genes).

### Score derivation

```python
expression_score = 1 / (1 + exp(5 × (chronos_score + 0.5)))
```

This sigmoid maps:
- Chronos 0 (non-essential) → score ≈ 0.08
- Chronos −0.5 → score ≈ 0.5
- Chronos −1.0 → score ≈ 0.92
- Chronos −1.5 (pan-essential) → score ≈ 0.99

**Why this sigmoid?** It captures the threshold-like nature of essentiality: most genes with Chronos > −0.3 are not meaningfully essential; genes with Chronos < −0.8 are almost certainly essential. The sigmoid provides smooth, continuous scores for composite scoring.

### Results

| Metric | Value |
|---|---|
| Genes with DepMap Chronos data | **5,574** |
| Expression result rows (DepMap) | 5,574 |

### Top Essential Genes (highest Chronos dependency, log2fc = 1.0 = max)

These are the most broadly essential genes across all ~1,100 cancer cell lines:

| Gene | Source | Biological role |
|---|---|---|
| DHX15_HUMAN | DepMap | RNA helicase; splicing |
| U2AF1_HUMAN | DepMap | Splicing factor; mutated in MDS |
| UBA1_HUMAN | DepMap | Ubiquitin-activating enzyme E1 |
| RCC1_HUMAN | DepMap | Regulator of chromosome condensation |
| SRSF3_HUMAN | DepMap | Serine/arginine-rich splicing factor |
| THOC2_HUMAN | DepMap | mRNA export complex |
| HCFC1_HUMAN | DepMap | Transcriptional co-activator |
| PUF60_HUMAN | DepMap | Splicing factor; DNA repair |

> These are expected pan-essential genes (splicing machinery, ubiquitin pathway, chromatin regulators) — their appearance as top-scored validates the DepMap data quality.

---

## Combined Expression Score

The three data sources (Open Targets, GTEx, DepMap) are combined into a single `expression_score` stored in `candidate_score`:

```python
expression_score = max(ot_score, gtex_score, depmap_score)
```

Taking the maximum rather than average means: a gene is scored by whichever evidence source provides the strongest signal for that specific gene. This is a **ceiling-union** approach — it rewards genes with any strong functional evidence without requiring all three databases to agree. The trade-off is that a broadly pan-essential gene (high DepMap score) will receive a high expression_score even if it has no cancer-specific OT association. Pan-essential genes are handled downstream by the Step 14 safety screen (DepMap Chronos floor), which zeroes the composite score for genes with `depmap_score > 0.7` across cell lines. The two mechanisms are therefore complementary rather than redundant.

### Final expression_score distribution

| Score range | Gene count | Interpretation |
|---|---|---|
| 0 (no evidence) | 7,020 | Not found in OT, GTEx, or DepMap |
| > 0 | **5,775** | Has at least one source of functional evidence |
| > 0.3 | **817** | Moderate-to-strong functional evidence |
| > 0.5 | **58** | Strong functional evidence (highly essential or well-linked to cancer) |
| Max (0.766) | RTF2_HUMAN | Strongest functional evidence in dataset |

**Summary statistics (non-zero genes only):**
- Average expression score: 0.166
- Maximum expression score: 0.766 (RTF2_HUMAN — replication termination factor)

### Expression Evidence by Source

| Source | Genes with data | Description |
|---|---|---|
| DepMap (essentiality) | 5,574 | CRISPR dependency in cancer cell lines |
| GTEx (tissue expression) | 5,529 | Human tissue expression breadth |
| Open Targets (disease) | 5,376 | Human gene–cancer evidence |

---

## Database State After Step 8

```
expression_result table:    16,479 rows
  DEPMAP:essentiality:       5,574 rows
  GTEX:tissue_expression:    5,529 rows  
  OT:disease_association:    5,376 rows

candidate_score.expression_score:
  5,731 genes with expression_score > 0  (DB confirmed April 2026)
  7,064 genes with expression_score = 0 (no functional evidence found)
```

---

## Interpretation Notes

**Why do only 5,731/12,795 genes have expression evidence?**

The 7,020 genes without expression scores either:
1. Have no Open Targets cancer association (novel candidates not yet linked to cancer in the literature)
2. Are not expressed in the 54 GTEx tissues (e.g., testis-specific or immune-cell-specific genes not in GTEx's tissue panel)
3. Were not included in the DepMap 24Q4 screen (some very small genes, pseudogenes, or poorly-annotated loci are excluded)

These genes are not penalised — they receive `expression_score = 0.0`, which means they contribute zero signal to the expression layer but can still rank based on selection + convergence alone.

---

## Next Step

→ [Step 9: Composite Scoring & Tier Assignment](step9_composite_scoring.md)
