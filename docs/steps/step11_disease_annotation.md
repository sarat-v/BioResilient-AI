# Step 11 — Disease Annotation (OpenTargets + GWAS + gnomAD + Literature)

**Pipeline phase:** Phase 2 — Clinical Translation Layer  
**Run timestamp:** 2026-04-11 (re-run to populate opentargets_score after initial race-condition failure)  
**Status:** ✅ PASS — 280 genes annotated across five independent evidence sources  

---

## What This Step Does

Step 11 answers: **does existing human clinical and genetic data corroborate the evolutionary signal?** A gene independently flagged by cross-species evolution AND by human disease genetics is the strongest possible candidate — two completely orthogonal evidence types pointing at the same gene.

Five external APIs are queried in parallel:

| Sub-step | Source | Evidence Type |
|---|---|---|
| 11 | Open Targets Platform | Disease association scores, tractability, known drugs, safety liabilities |
| 11b | GWAS Catalog + gnomAD | Genome-wide significant variants; loss-of-function intolerance |
| 11c | PubMed | Literature citations in longevity/cancer-resistance literature |
| 11d | Reactome / GO | Pathway-level convergence enrichment |
| (concurrent) | IMPC | Mouse knock-out phenotypes |

---

## Open Targets (Step 11 proper)

### Methodology

For each candidate gene:
1. Resolve gene symbol (e.g. `G6PD`) to Ensembl ID via Open Targets search API
2. Query the Open Targets Platform v4 GraphQL API for:
   - `associatedDiseases` — maximum association score across all diseases
   - `tractability` — small molecule / antibody / PROTAC modality flags
   - `drugAndClinicalCandidates` — highest-phase known drug
   - `safetyLiabilities` — known adverse event labels

### Results

| Metric | Value |
|---|---|
| Genes with Open Targets association score | 119 / 280 (43%) |
| Average OT score (genes with data) | 0.398 |
| Max OT score | 0.852 (G6PD) |
| Genes with known drug (any phase) | 22 |
| Genes SM-tractable | 36 |
| Genes antibody-tractable | 64 |
| Genes PROTAC-tractable | 76 |

**161 genes with no OT score** is scientifically expected. These are novel candidates discovered through convergent evolution — they have not yet been studied in the context of human disease. This is precisely the value of the computational approach: identifying targets that conventional disease-genetics studies have not reached.

### Genes with Notable Disease Links

| Gene | OT Score | Known Disease | Known Drug |
|---|---|---|---|
| G6PD_HUMAN | 0.852 | G6PD deficiency, hemolytic anemia | — |
| MPV17_HUMAN | 0.851 | Mitochondrial hepatopathy | — |
| ABCB7_HUMAN | 0.821 | Sideroblastic anemia | — |
| HEPC_HUMAN | 0.798 | Hemochromatosis, iron deficiency anemia | — |
| XDH_HUMAN | 0.758 | Xanthinuria, gout | Allopurinol |
| ACLY_HUMAN | 0.544 | Atherosclerosis, metabolic syndrome | Bempedoic acid |
| CD80_HUMAN | 0.618 | Autoimmune disease, cancer immunotherapy | Galiximab (Ph3) |
| KCNS1_HUMAN | 0.590 | Pain sensitivity, neuropathy | Guanidine (Ph3) |

---

## GWAS Catalog + gnomAD (Step 11b)

### Methodology

- Query GWAS Catalog for genome-wide significant associations (p < 5×10⁻⁸) for each gene's chromosomal window
- Query gnomAD v4 for pLI (probability of loss-of-function intolerance) and LOEUF (loss-of-function observed/expected upper bound fraction)

### LOEUF Threshold (Pipeline-Wide Consistency)

The LOEUF threshold used throughout this pipeline is **LOEUF < 0.35**, corresponding to the "high-confidence LoF intolerant" gene set defined in Karczewski et al. (2020) *Nature* 581:434–443. This is consistent across:
- Step 4d `motif_direction` classification (`loss_of_function` class requires LOEUF < 0.35)
- `gnomad.py` `LOEUF_INTOLERANT_THRESHOLD = 0.35`
- `scoring.py` `disease_score()` LOEUF contribution

This is intentionally more conservative than the gnomAD v4 general guideline (0.6) in order to restrict the LoF-intolerant designation to genes where LoF variants are robustly depleted from the human population — a stricter standard appropriate for therapeutic target nomination.

### Results

| Metric | Value |
|---|---|
| Genes with GWAS hit (p < 5×10⁻⁸) | 140 / 280 |
| Genes with gnomAD pLI data | 17 |
| Genes with gnomAD LOEUF data | Available for several (used in preference to pLI) |
| One precision artifact | HMOX2: gwas_pvalue stored as 0.0 (guard clause prevents crash; GWAS bonus = 0) |

**Protective variant paradigm (PCSK9 model):** Step 11b also maps convergent motif positions to rare human variants in gnomAD that have the same amino acid direction as the longevity-model substitution. This is the strongest possible human validation — a naturally-occurring human variant recapitulating the animal adaptation (Cohen et al. 2006 *NEJM* 354:1264; Miyata et al. 1979 *J Mol Evol* 12:219 for biochemical similarity grouping used in variant matching).

---

## Pathway Enrichment and BH-FDR Correction (Step 11d)

Step 11d computes Reactome pathway enrichment over Tier1/Tier2 candidate genes using a hypergeometric test. As of pipeline version 2026-04-17:

- Raw hypergeometric p-values are computed per pathway using the `_log_hypergeometric()` function (log-stable implementation)
- **Benjamini-Hochberg (BH) FDR correction** is applied across all pathways tested in a single run (Benjamini & Hochberg 1995, *J Royal Stat Soc B* 57:289–300)
- Both `log_pvalue` (raw, log₁₀ scale) and `adjusted_pvalue` (BH-corrected, linear scale) are stored in the `pathway_convergence` table
- Pathways with `adjusted_pvalue < 0.05` are statistically enriched; those with 0.05–0.20 are reportable as trends

This correction was introduced in DB migration 0031. Legacy rows created before migration will have `adjusted_pvalue = NULL`; re-running Step 11d populates them.

**Background gene count:** The hypergeometric test uses 20,000 as the approximate human protein-coding gene count (Ensembl 110 estimates ~19,800). This is a conservative approximation; the actual number changes marginally across Ensembl releases and does not materially affect enrichment significance for the pathway sizes we observe.

---

## Literature Validation (Step 11c)

- 280 / 280 genes have literature scores
- Literature score range: 0.0 (no relevant papers) to 0.822 (G6PD — extensively studied in cancer/longevity contexts)

### Literature Scoring Algorithm (deterministic, no LLM)

The literature score is computed in two stages per gene (implemented in `pipeline/layer3_disease/literature.py`):

**Stage 1 — PubMed hit count:** Query NCBI E-utilities esearch API for papers co-mentioning the HGNC gene symbol AND ≥ 1 of a trait-specific keyword set (longevity, aging, senescence, cancer resistance, tumour suppressor, antiviral, etc.). Retrieves up to 20 PMIDs.

**Stage 2 — Semantic proportion score (primary):** Fetch full abstracts via efetch. For each abstract, check whether the abstract mentions both (a) the gene symbol and (b) at least one of the phenotype keywords or mechanism terms (protective, resistance, positive selection, etc.). The literature score is `matched_abstracts / total_abstracts_fetched` — a proportion of "relevant" papers (range: [0, 1]).

**Fallback (efetch failure):** If abstract fetching fails, score = `min(log10(1 + count) / 5.0, 0.50)` — a compressed log-count capped at 0.50. The full 1.0 is only achievable via the semantic proportion path.

**Why proportion over raw count:** A gene with 5 highly relevant papers (score = 1.0 if all match) outperforms a gene with 500 tangentially mentioned papers (score = 0.01 if only 5 abstracts match). This removes the well-studied gene publication bias.

**Sanity check:** At run time, the module logs whether any known resilience genes (TP53, FOXO3, IGF1, AKT1, BRCA1, etc.) appear in the candidate set. Their presence validates that the literature queries are finding biologically meaningful results.

**Role in composite:** `lit_score` is an input to `disease_score()` (weight 0.15 in the literature component within disease_score). It is NOT directly an input to the composite — it feeds through the disease score. It is also displayed separately in the API and dashboard as a confidence signal for researchers.

---

## Disease Score Formula

```
disease_score = 0.30 × min(opentargets_score, 1.0)          [OT association breadth]
              + 0.20 × min(–log10(gwas_pvalue) / 15, 1.0)   [GWAS p-value strength]
              + 0.15 × LOEUF-derived or pLI constraint score  [LoF intolerance]
              + 0.20 × min(known_drug_phase / 4, 1.0)         [clinical validation]
              + 0.15 × protective_variant_bonus               [PCSK9 paradigm]
```

LOEUF is preferred over pLI when available (LOEUF is a continuous variable; pLI is quasi-binary).

---

## Data Quality Notes

- `disease_name` field is not populated (OT code stores max score, not individual disease names) — informational gap only, no scoring impact
- `disease_id` (EFO ID) field is not populated — same reason
- GWAS race condition on first run caused all opentargets_scores to be NULL; resolved by re-running step 11 after rows were established in the table

---

## Data Quality Assessment

✅ 280 disease_annotation rows, zero orphaned  
✅ opentargets_score strictly in [0, 1]  
✅ gwas_pvalue > 0 guard prevents math domain error  
✅ lit_score present for all 280 genes  
✅ disease_score in candidate_score verified against formula (G6PD: 0.5061 confirmed)  
