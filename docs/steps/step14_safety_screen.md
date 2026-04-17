# Step 14 — Safety Screen (PheWAS + STRING + DepMap + GTEx)

**Pipeline phase:** Phase 2 — Clinical Translation Layer  
**Run timestamp:** 2026-04-11  
**Status:** ✅ PASS — 280 genes screened; safety floor applied in Step 15  

---

## What This Step Does

Step 14 applies a clinical safety lens to candidate genes before final scoring. A gene that is highly conserved across longevity-model species and druggable is still not a viable target if modulating it would cause serious side effects in humans. Safety screening identifies genes that are:

- **Broadly essential** for normal cell viability (DepMap)
- **Widely expressed** across many tissues (GTEx) — high off-target risk
- **Network hubs** — highly connected in protein interaction networks (STRING)
- **PheWAS flagged** — associated with adverse phenotypes in UK Biobank / BioBank Japan

Safety data informs a **multiplicative floor** in the composite score, not a simple additive penalty. Any gene with `safety_score < 0.40` has its composite score zeroed regardless of evolutionary signal. This prevents high-scoring unsafe targets from appearing in the final candidate list.

---

## Evidence Sources

### Step 14 — PheWAS + STRING (safety.py)

- **PheWAS (phenome-wide association study):** Query UK Biobank / BioBank Japan for any association between common variants near each gene and adverse phenotypes (e.g., renal failure, cardiac events)
- **STRING:** Count high-confidence protein-protein interactions (score > 700) as a proxy for hub status — highly connected proteins are more likely to have widespread effects if modulated
- **IMPC:** Mouse knock-out lethality / severe phenotype as an essentiality proxy

### Step 14b — DepMap + GTEx (depmap_gtex.py)

- **DepMap CRISPR Chronos scores:** Cancer cell line essentiality. Chronos score → 1.0 means knocking out the gene kills all cell lines (pan-essential). Score < 0.5 = selectively essential or non-essential.
- **GTEx:** Expression breadth (how many tissues) and maximum TPM — very broadly expressed genes have higher off-target risk

---

## Results

### Safety Flag Summary (safety_flag table)

| Metric | Value |
|---|---|
| Total genes screened | 280 |
| Genes flagged essential (is_essential = True) | 2 (AKT1, ARPC4) |
| Genes flagged hub-risk | 41 / 280 (15%) |
| Genes with DepMap data | 154 / 280 (55%) |
| Average DepMap Chronos score | 0.308 |
| High DepMap risk (score > 0.5) | 43 genes |
| Genes with GTEx expression data | 154 / 280 |
| Average tissue expression count | 41.1 tissues |
| Widely expressed genes (> 20 tissues) | 120 / 154 with data |

### Essential Genes

| Gene | is_essential | DepMap Score | Composite | Outcome |
|---|---|---|---|---|
| AKT1_HUMAN | ✅ True | 0.099 | 0.406 | At safety floor (0.4) — remains in list with caution flag |
| ARPC4_HUMAN | ✅ True | 0.706 | 0.000 | **Zeroed out** — pan-essential cytoskeletal component |

**ARPC4** (actin-related protein 2/3 complex subunit 4) is correctly excluded — it is a core component of branched actin networks essential for cell division in virtually all cell types. Its high convergence score in longevity models likely reflects pleiotropic evolutionary effects, not a viable therapeutic target.

**AKT1** is at the safety floor boundary. It is a validated cancer drug target (multiple PI3K pathway inhibitors in clinical use), but essentiality requires tissue-specific delivery approaches. The composite score of 0.406 (floor value) correctly signals: high interest, high caution.

### 126 Genes Missing DepMap Data

These genes are not represented in the DepMap cancer cell line dataset (Broad Institute DepMap 23Q4). This is expected for:
- Very recently characterized proteins with no gene dependency data
- Genes expressed primarily in post-mitotic cells not well-represented in cancer lines
- Pseudogenes or non-coding gene annotations

**Default `safety_score = 0.5` for missing DepMap data** is a conservative neutral prior. This means that genes with no essentiality evidence are treated as "unknown risk" — they can still reach the final candidate list but are flagged for experimental essentiality screening. A higher default (e.g., 0.7) would be too optimistic; a lower one (e.g., 0.3) would incorrectly penalise novel proteins that happen to not be in the cancer line library.

**Rationale for the 0.5 default:** DepMap CRISPR covers ~18,000 genes but is enriched for cancer-relevant proteins. Proteins absent from DepMap are more likely to be tissue-restricted, developmentally regulated, or low-confidence annotations — all of which argue for caution rather than dismissal. The floor (0.40) is still enforced; no DepMap-absent gene can exceed a safety score of 1.0 × (all other terms passing), so the default itself is not a backdoor to a falsely safe score.

---

## Safety Score Formula

```
safety_score = 1.0                           [starting value: safe]
             − 0.30 × is_essential           [broad lethality]
             − 0.20 × (depmap_score > 0.7)   [pan-essential in cancer lines]
             − 0.10 × hub_risk               [high STRING connectivity]
             − 0.10 × (gtex_tissue_count > 30 AND gtex_max_tpm > 100)  [broad expression]
             
if safety_score < 0.40: composite_score = 0  [hard floor — excluded from ranking]
```

---

## Phase 2 Composite: Safety as a Multiplicative Floor

Unlike all other Phase 2 evidence layers (which additively contribute to the composite), safety acts as a **multiplicative floor**:

- If `safety_score ≥ 0.40`: composite proceeds normally
- If `safety_score < 0.40`: `composite_score = 0` regardless of evolutionary or clinical evidence

This design choice prevents a scenario where a highly essential gene (pan-lethal in all cell types) scores highly simply because it evolves rapidly in longevity models.

**Result:** 0 Tier1/Tier2/Validated genes have safety_score < 0.40. The floor is enforced without eliminating any high-confidence candidates — the two essential genes (AKT1 floor-limited, ARPC4 zeroed) are correctly handled.

---

## Data Quality Assessment

✅ 280 safety_flag rows, 0 orphaned  
✅ 0 Tier1/Tier2/Validated genes with safety_score < 0.40 and non-zero composite  
✅ DepMap Chronos scores in [0, 1]  
✅ GTEx tissue counts plausible (0–54 GTEx tissues)  
✅ Safety floor correctly applied in Step 15 rescore  
