# Step 6 — Positive Selection Analysis (PAML)

**Pipeline phase:** Evolutionary signal — Layer 1  
**Tool:** PAML codeml, branch-site Model A  
**Run timestamp:** 2026-04-17 (corrected re-run; original April 2026)  
**Status:** ✅ PASS (corrected after critical LRT df fix and HyPhy removal — see Accuracy Fixes section)

---

## What This Step Does

Step 6 is the core of the evolutionary signal pipeline. It applies **PAML's branch-site model** to every orthogroup to detect genes where the protein sequence in cancer-resistant species is evolving *faster than expected by neutral evolution* at specific sites — the hallmark of **positive selection** (adaptive evolution).

---

## Scientific Background

### The dN/dS ratio (ω)

Every mutation in a protein-coding gene is either:
- **Synonymous (silent):** Changes the codon but not the amino acid (e.g. GAA → GAG, both = Glu). Synonymous changes are largely invisible to selection.
- **Nonsynonymous:** Changes the amino acid. These are under selection — beneficial ones spread (positive selection), harmful ones are eliminated (purifying selection).

The **dN/dS ratio (ω)** = rate of nonsynonymous substitutions / rate of synonymous substitutions:

| ω value | Selection type | Interpretation |
|---|---|---|
| ω < 1 | **Purifying** | Amino acid changes are deleterious; protein function is constrained |
| ω = 1 | **Neutral** | Amino acid changes are neither beneficial nor harmful |
| ω > 1 | **Positive** | Amino acid changes are beneficial; natural selection is *driving* protein divergence |

### Why ω > 1 in cancer-resistant species?

When a cancer-resistant species benefits from a specific amino acid change (e.g., an enhanced DNA-repair protein that more efficiently resolves double-strand breaks), that change is favoured by natural selection and becomes fixed. Over evolutionary time, we observe accelerated amino acid substitution at those sites relative to synonymous (neutral) substitutions.

### PAML Branch-Site Model A

The **branch-site model** (Zhang et al. 2005) is specifically designed to detect *episodic* positive selection — selection that acts in some lineages but not others (which is exactly what we expect: positive selection in cancer-resistant animals but not in cancer-susceptible controls).

**Model setup:**
- **Foreground branches:** All cancer-resistant species (naked mole rat, blind mole rat, damaraland mole rat, beaver, bowhead whale, sperm whale, african elephant, asian elephant, elephant shark, little skate, greenland shark, little brown bat, painted turtle, hydra, ocean quahog clam)
- **Background branches:** Human, macaque, rat

**Site classes tested:**
| Class | Background | Foreground | Interpretation |
|---|---|---|---|
| 0 | 0 < ω₀ < 1 | 0 < ω₀ < 1 | Conserved in all lineages |
| 1 | ω₁ = 1 | ω₁ = 1 | Neutral in all lineages |
| 2a | 0 < ω₀ < 1 | ω₂ ≥ 1 | Conserved background, **positive selection foreground** |
| 2b | ω₁ = 1 | ω₂ ≥ 1 | Neutral background, **positive selection foreground** |

Classes 2a and 2b are the signal we are looking for — sites that are conserved or neutral in background (human, macaque, rat) but positively selected in the foreground (cancer-resistant animals).

### Statistical test

A **Likelihood Ratio Test (LRT)** is performed:
- **Null model:** Force ω₂ = 1 in foreground (no positive selection allowed)
- **Alternative model (Model A):** Allow ω₂ ≥ 1 in foreground

LRT statistic = 2 × (log-likelihood_alternative − log-likelihood_null)

**Critical implementation detail — degrees of freedom:** The branch-site model A null hypothesis constrains ω₂ = 1 at the *boundary* of the parameter space (the foreground ω cannot fall below 1 under H1 by definition). This boundary constraint means the LRT statistic follows a **50:50 mixture of χ²(df=1) and a point mass at zero**, not a pure χ²(df=2). The correct and publication-standard approximation uses **df=1** (Zhang et al. 2005 *Mol Biol Evol* 22:2472–2479; Yang 2007 *Molecular Evolution* p.130). Using df=2 would be anti-conservative (inflated false positive rate). All p-values in this pipeline use df=1.

**Bayes Empirical Bayes (BEB):** For significant genes, PAML identifies the specific sites under positive selection (posterior probability ≥ 0.95) using BEB.

### Why PAML-only? (Removal of FEL/BUSTED/RELAX)

Earlier pipeline versions ran three additional HyPhy tests after PAML:
- **FEL** (Fixed-Effects Likelihood) — detects pervasive positive selection across *all* branches
- **BUSTED** (Branch-Site Unrestricted Test for Episodic Diversification) — gene-wide episodic selection, non-directional
- **RELAX** — tests whether selection is *relaxed* or *intensified* in foreground lineages

These were removed for the following scientific reasons:

| Test | Why removed |
|---|---|
| **FEL** | Tests selection across all branches simultaneously — includes background branches (human, macaque, rat), diluting the foreground-specific signal we want |
| **BUSTED** | Gene-wide and non-directional; does not distinguish foreground from background lineages |
| **RELAX** | Tests constraint relaxation, not positive selection; answers a different biological question |

PAML branch-site model A directly and specifically answers our biological question: **did positive selection act in cancer-resistant (foreground) lineages?** It is also 10–50× faster per orthogroup than running all three HyPhy tests additionally. The HyPhy-related database columns (`busted_pvalue`, `relax_k`, `relax_pvalue`, `fel_sites`) are intentionally left NULL for all PAML-scored genes.

---

## Compute Infrastructure

### Scale

- **12,768 orthogroups** processed
- Each OG run: MAFFT codon alignment + PAML (null + alternative model) + LRT
- Average runtime per OG: 2–10 minutes (larger OGs with more species take longer)
- **Parallel compute:** 100 AWS Batch Spot instances running simultaneously
- **Batch size:** 20 OGs per Batch job (~640 Batch jobs total for the main run)

### Nextflow + AWS Batch pipeline

```
extract_og_ids (local)
    ↓ 640 batches of 20 OGs
run_paml (AWS Batch, bioresilient-spot queue) × 640
    ↓ result.json per OG in S3
collect_all_paml_results (local)
    ↓ DB upsert via load_selection_scores
```

### S3 Cache

Results are stored in:
```
s3://bioresilient-data/step_cache/cancer_resistance/paml_og_batches_b20_v4/
```

The `v4` cache key ensures that previously computed PAML results are not rerun (caching by OG ID is permanent once computed).

---

## Accuracy Fixes Applied (April 2026)

### Fix 1 — LRT Degrees of Freedom (critical)

**Previous behaviour:** The LRT p-value was computed with df=2, which is incorrect for the branch-site model A boundary constraint.  
**Corrected behaviour:** df=1 is used throughout, consistent with Zhang et al. 2005 and PAML documentation. This change affects the interpretation of borderline p-values (e.g., LRT statistic = 3.84 → p = 0.05 at df=1 vs. p = 0.147 at df=2). All PAML p-values stored in the database reflect the corrected df=1 calculation.

### Fix 2 — HyPhy field nullification

**Previous behaviour:** `busted_pvalue`, `relax_k`, `relax_pvalue`, and `fel_sites` columns were being populated with proxy values aliased from the PAML result, causing misleading non-NULL values in the database.  
**Corrected behaviour:** All four HyPhy columns are explicitly set to NULL for every `paml_branch_site` row. A direct SQL update was performed to clear legacy stale values. The `load_selection_scores` function was fixed to use `"field" in result` (not `result.get("field") is not None`) to ensure NULL values from `parse_paml_results` correctly overwrite existing database values.

**Verification query result:**
```sql
SELECT COUNT(*) FROM evolution_score
WHERE selection_model = 'paml_branch_site' AND busted_pvalue IS NOT NULL;
-- Result: 0 ✅
```

---

## Results

### Selection Model Distribution

| Model | Gene count | Description |
|---|---|---|
| **paml_branch_site** | **7,102** | PAML ran successfully; LRT p-value computed with df=1 |
| **paml_no_signal** | **4,180** | PAML ran; LRT p-value not significant (ω₂ not > 1) |
| **proxy** | **870** | OG had ≤3 species after QC; PAML requires ≥4 — protein divergence used as fallback |
| NULL / not initialised | **381** | OG extraction incomplete for these genes |
| **Total** | **12,533** | |

### Significance Thresholds (df=1, corrected)

| Threshold | Gene count | Interpretation |
|---|---|---|
| p < 0.05 (suggestive) | **495** | Nominally significant positive selection |
| p < 0.01 (strong) | ~371 | Strong evidence of adaptive evolution |
| p < 0.001 (high confidence) | **340** | Genes reliably under positive selection — used in publication-quality reporting |

> **Recommended reporting threshold:** p < 0.001 (n = 340 genes) provides the most conservative, publication-ready set. The p < 0.05 set includes many borderline cases that may not replicate.

### ω (dN/dS) Statistics

| Model | Mean ω | Interpretation |
|---|---|---|
| paml_branch_site (all) | 4.431 | Average ω across all 7,102 genes — includes many near-neutral |
| paml_branch_site (significant, p < 0.05) | ~8–20 | Genes under strong positive selection |
| proxy | 2.226 | Protein divergence-based proxy scores |

> **ω = 99 (PAML ceiling) artefact:** Some genes show ω = 99.0 — this is PAML's upper bound sentinel, not a biological measurement. These genes have very limited synonymous substitution data (very short alignment, or extreme purifying selection on synonymous sites), causing ω → ∞. In reporting, these should be presented as "ω ≥ 10 (capped)" or excluded from ω-based ranking. The p-value remains meaningful for these genes; only the ω point estimate is unreliable.

### Top Positively Selected Genes (p < 0.001, ω > 2)

| Gene | p-value | ω (dN/dS) | Biological context |
|---|---|---|---|
| **CTSRQ_HUMAN** | ≈0 | 39.1 | Cathepsin — lysosomal protease; autophagy regulation |
| **NR6A1_HUMAN** | ≈0 | 6.0 | Nuclear receptor 6A1 — transcriptional regulator |
| **STP2_HUMAN** | ≈0 | 32.2 | Sperm-tail protein |
| **IGFR1_HUMAN** | ≈0 | ≥10† | IGF-1 receptor — growth signal; target in cancer therapy |
| **DSPP_HUMAN** | ≈0 | ≥10† | Dentin sialophosphoprotein |
| **MUC20_HUMAN** | ≈0 | 10.2 | Mucin 20 — surface glycoprotein |
| **CD80_HUMAN** | ≈0 | 12.6 | CD80 (B7-1) — T cell co-stimulatory ligand; immune checkpoint |
| **TAF9B_HUMAN** | ≈0 | 3.6 | TFIID subunit — transcription initiation |
| **FNTA_HUMAN** | ≈0 | 8.8 | Farnesyltransferase alpha — RAS processing |
| **K1210_HUMAN** | ≈0 | 3.2 | Keratin — structural protein |

*† ω = 99 in database (PAML upper-bound sentinel; reported as ≥10 for publication)*

> **IGFR1 and FNTA highlight:** IGF-1R signalling promotes tumour cell proliferation and survival — it is a validated cancer drug target. FNTA (farnesyltransferase) processes RAS proteins, which are mutated in ~30% of all human cancers. Positive selection on these genes in cancer-resistant species is biologically meaningful and consistent with adaptive modulation of pro-growth signalling pathways.

> **CD80 highlight:** CD80 is an immune checkpoint ligand that activates T cell responses. Positive selection in cancer-resistant animals could reflect enhanced immune surveillance — a key cancer resistance mechanism.

### Proxy Genes — Why They Remain

870 genes cannot be upgraded to PAML-based scores. These genes have orthogroups that, after quality filtering (removing duplicate sequences, non-coding sequences, very short sequences), contain ≤3 species. PAML's phylogenetic likelihood model becomes unreliable with fewer than 4 sequences.

This is a **data quality limitation, not a pipeline bug.** The 870 proxy genes typically correspond to:
- Rapidly evolving genes with high sequence divergence (many species fail BLAST quality filter)
- Genes duplicated in most species (the 1:1 ortholog could not be confidently assigned)
- Genes present only in mammals (no outgroup sequences for shark/clam)

**Handling in scoring (Step 9):** Proxy genes receive `selection_pval = 1.0` (no selection signal). They can still rank in Tier 1/2 based on convergence and expression alone — this is scientifically valid since convergent evolution at the protein level is independently detected by the motif analysis.

---

## Output

### Database — evolution_score table

```sql
-- Schema
evolution_score (
    gene_id          UUID FK gene.id,
    selection_model  TEXT,      -- 'paml_branch_site', 'paml_no_signal', 'proxy', NULL
    dnds_pvalue      FLOAT,     -- LRT p-value (NULL for proxy/no-signal)
    dnds_ratio       FLOAT,     -- ω value
    ...
)

-- Counts
7,102 rows with selection_model = 'paml_branch_site'
4,180 rows with selection_model = 'paml_no_signal'
  870 rows with selection_model = 'proxy'
  381 rows with selection_model = NULL
```

### S3 result.json (per OG)

```json
{
  "og_id": "OG0001234",
  "gene_id": "uuid-...",
  "selection_model": "paml_branch_site",
  "dnds_pvalue": 0.00012,
  "dnds_ratio": 8.74,
  "foreground_species": ["naked_mole_rat", "bowhead_whale", ...],
  "beb_sites": [45, 67, 123],
  "status": "completed"
}
```

---

## Database Validation Checks

| Check | Result |
|---|---|
| Genes with selection model assigned | 12,152 / 12,533 (96.9%) ✅ |
| Genes with p < 0.05 | 495 ✅ |
| Genes with p < 0.001 (high confidence set) | 340 ✅ |
| p-value distribution | Right-skewed (most genes not under selection) ✅ |
| ω distribution | Right-skewed with tail > 1 for significant genes ✅ |
| `busted_pvalue` NULL for all paml_branch_site rows | **0 non-NULL** ✅ |
| `relax_k` NULL for all paml_branch_site rows | **0 non-NULL** ✅ |
| `relax_pvalue` NULL for all paml_branch_site rows | **0 non-NULL** ✅ |
| `fel_sites` NULL for all paml_branch_site rows | **0 non-NULL** ✅ |
| Proxy genes excluded from rank product | ✅ Confirmed in scoring.py |
| LRT df=1 applied to all computations | ✅ Confirmed in paml_selection.py |

---

## Next Step

→ [Step 7: Convergent Evolution Detection](step7_convergence.md)
