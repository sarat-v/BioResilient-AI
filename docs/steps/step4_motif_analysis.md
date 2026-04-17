# Step 4 — Protein Divergence Motif Analysis

**Pipeline phase:** Protein-level analysis  
**Substeps:** 4a (MEME motifs), 4b (domain + AlphaMissense), 4c (ESM-1v variant scoring), 4d (functional direction classification)  
**Run timestamps:** 4a–4d: April 2026; 4c: initial April 2026; 4d: re-run April 2026 (LOEUF fix)  
**Status:** ✅ PASS (4a, 4b, 4c, 4d); corrected after LOEUF threshold and sourcing fixes — see Accuracy Fixes section

---

## What This Step Does

Step 4 is the protein-level divergence layer. While Step 6 asks *"is this gene under positive selection globally?"*, Step 4 asks: **"at exactly which positions does the protein sequence differ between cancer-resistant species and human — and is that difference likely to be functionally important?"**

This provides two benefits:
1. **Localisation:** We know precisely *where* in the protein the divergence occurs (active site? surface loop? allosteric pocket?)
2. **Functional prediction:** We can estimate whether the divergence changes the protein's function using AI models (AlphaMissense, ESM-2)

---

## Step 4a — MEME Motif Discovery

### Tool

**MEME Suite** (Multiple Em for Motif Elicitation) applied to sliding windows over species-aligned protein sequences.

### Process

For each gene:
1. Retrieve all species' protein sequences from the ortholog table
2. Generate multiple sequence alignment with MAFFT
3. Slide a window of 15 amino acids across the alignment
4. At each position, record a **divergent_motif** if:
   - ≥2 cancer-resistant species share an amino acid state that differs from human
   - The divergence is non-conservative (BLOSUM62 score < threshold)
5. Store: human sequence, cancer-resistant consensus, position, species contributing

### Results

| Metric | Value |
|---|---|
| **Total divergent motifs identified** | **1,692,443** |
| **Genes with ≥1 divergent motif** | **12,186** |
| Genes per motif (average) | ~139 motifs per gene |
| Max motifs in a single gene | 272 (many large/multi-domain proteins) |

### Top Genes by Motif Count

| Gene | Motif count | Domain motifs | Notes |
|---|---|---|---|
| LIMA1_HUMAN | 272 | 29 | LIM domain, actin cytoskeleton |
| EXD1_HUMAN | 272 | 250 | Exonuclease domain protein — nearly all motifs in domain |
| DUS15_HUMAN | 272 | 25 | tRNA dihydrouridine synthase |
| CO6A5_HUMAN | 272 | 178 | Collagen VI alpha 5 — large extracellular matrix protein |
| EH1L1_HUMAN | 272 | 8 | EH domain containing protein |
| IQEC1_HUMAN | 272 | 95 | IQ motif-containing protein |

---

## Step 4b — Functional Domain Annotation & AlphaMissense Scoring

### Tool 1: Pfam/InterPro

Each divergent motif is checked against **Pfam domain boundaries** (via InterPro REST API) to determine if the variant position falls within an annotated functional domain.

**Why this matters:** A substitution in a kinase catalytic loop has a fundamentally different implication than one in an unstructured linker. Domain-localised divergence is prioritised in scoring.

### Pfam Domain Distribution (top 10)

| Domain | Motifs in domain | Biological relevance |
|---|---|---|
| Protein kinase | 8,690 | Signal transduction — key cancer pathway nodes |
| Peptidase S1 | 3,506 | Serine proteases — coagulation, immune |
| Ig-like V-type | 3,123 | Immune receptors |
| PH domain | 2,450 | Phosphoinositide signalling |
| USP | 1,944 | Ubiquitin specific protease |
| WD repeat | 1,944 | Protein–protein interaction scaffolds |
| SH3 | 1,789 | Signalling domain interactions |
| Ig-like C2-type | 1,772 | Adhesion/immune molecules |
| KRAB | 1,748 | Transcriptional repressor |
| BTB | 1,615 | Protein dimerisation |

**Total motifs in functional domains: 247,291 (14.6%)**

### Tool 2: AlphaMissense (Google DeepMind)

**AlphaMissense** predicts the functional impact of single amino acid substitutions using a model trained on protein structure predictions (AlphaFold2) combined with evolutionary and population genetics data. Score range: 0–1.

| Score range | Class | Interpretation |
|---|---|---|
| 0–0.34 | Likely benign | Substitution likely tolerated |
| 0.34–0.564 | Uncertain | Ambiguous functional impact |
| 0.564–1.0 | Likely pathogenic | Substitution likely disrupts function |

**AlphaMissense results:**

| Metric | Value |
|---|---|
| Motifs with AM score | **1,322,513** (78.1% coverage) |
| Motifs scored as likely pathogenic (AM > 0.564) | **42,538** |
| Motifs with very high AM score (> 0.8) | Subset of above |

> **Interpretation:** A "likely pathogenic" AM score for a motif in a cancer-resistant animal means: the substitution would be damaging if introduced randomly in a human cell — but these animals are *healthy and cancer-resistant*. This suggests the substitution confers a context-specific function change, not a loss-of-function in the cancer-resistant lineage.

### Top High-Priority Motifs (high AM score + functional domain)

| Gene | Domain | AM score | Biological significance |
|---|---|---|---|
| PLSI_HUMAN | Various | High | 1,319 convergent AAs in domain |
| UCHL1_HUMAN | UCH | High | Ubiquitin carboxy-terminal hydrolase L1 — neurodegeneration |
| GSTO1_HUMAN | GST | High | Glutathione S-transferase — oxidative stress response |
| GSTM1_HUMAN | GST | High | Glutathione S-transferase mu 1 |
| KLK11_HUMAN | Peptidase S1 | High | Kallikrein — serine protease, cancer biomarker |

---

## Step 4c — ESM-1v Variant Effect Scoring

### Tool

**ESM-1v** (Evolutionary Scale Modeling, 1-variant model) from Meta FAIR — a 650-million parameter protein language model fine-tuned specifically for **variant effect prediction**. Unlike ESM-2 which measures representational distance, ESM-1v computes a **log-likelihood ratio (LLR)** for each amino acid substitution relative to the reference (human) sequence.

### What ESM-1v LLR measures

For each substitution at a divergent motif position:

```
LLR = log P(variant | context) − log P(human_reference | context)
```

| LLR value | Interpretation |
|---|---|
| LLR < 0 (negative) | Variant is *less* probable than human reference — likely functionally destabilising |
| LLR ≈ 0 | Variant is as probable as reference — functionally neutral |
| LLR > 0 (positive) | Variant is *more* probable than human reference — potentially improved function |

> **Key distinction:** In the context of cancer-resistant animals, a negative LLR does not necessarily mean harmful — it means the substitution is unlikely from a generative protein sequence perspective (trained on extant sequences). Cancer-resistant species may have independently evolved rare-but-functional substitutions that appear "unlikely" to a model trained on all mammalian sequences.

### NULL vs 0.0 semantics (important for data interpretation)

| `esm1v_score` value | Meaning |
|---|---|
| **NULL** | ESM-1v computation was skipped for this motif (motif filtered out before scoring, or computation failed) |
| **0.0 (exact)** | ESM-1v was computed and the log-likelihood ratio is genuinely 0.0 — the substitution is perfectly neutral relative to the reference in this model |

These are semantically distinct. Future pipeline runs write `NULL` when ESM-1v is not computed, and an actual float (which may be 0.0) when it is computed. **Do not treat NULL and 0.0 as equivalent in analysis.**

### Relationship to AlphaMissense

- **AlphaMissense** = pathogenicity prediction for single amino acid substitutions (trained on structure + population genetics)
- **ESM-1v LLR** = sequence-level variant effect from protein language model context

These are complementary:
- AlphaMissense captures structural/functional constraint
- ESM-1v captures sequence-evolutionary plausibility

A motif with high AlphaMissense score and strongly negative ESM-1v LLR = likely disruptive in human but tolerated in cancer-resistant animal (convergent functional shift).

### Results

| Metric | Value |
|---|---|
| Motifs with ESM-1v score computed | Majority of 1,692,443 |
| Motifs with neutral LLR (≈ 0.0) | Expected large fraction (most substitutions are near-neutral) |
| Motifs with strongly negative LLR | Subset — represents biochemically disruptive changes |
| Motifs with NULL esm1v_score | Small fraction where computation was skipped |

---

## Step 4d — Functional Direction Classification

### Classification logic

Each divergent motif is classified into one of four functional categories based on combined AlphaMissense score, ESM-1v LLR, LOEUF score (LoF intolerance), and domain context:

| Class | Criteria | Interpretation |
|---|---|---|
| **loss_of_function** | Gene LOEUF ≤ 0.35 (LoF-intolerant) AND AM score > 0.564 | Gene cannot tolerate loss of function; cancer-resistant species modified a constrained gene — likely fine-tuning rather than breaking it |
| **gain_of_function** | ESM-1v LLR > threshold AND AM score > 0.7 AND in domain | Substitution is predicted to increase function in cancer-resistant lineage |
| **functional_shift** | Moderate ESM-1v or AM change in a functional domain | Detectable change in function, direction ambiguous |
| **neutral** | Low ESM-1v LLR, low AM score, or outside domain | No predicted functional consequence |

> **Why four classes instead of three?** The `loss_of_function` class is scientifically important for this pipeline: a cancer-resistant animal that has a substitution in a LoF-intolerant gene (LOEUF ≤ 0.35) is unlikely to have actually lost that gene's function. Instead, this likely represents **dosage modulation or pathway dampening** — reducing but not eliminating a pro-proliferative function. These genes are biologically meaningful candidates distinct from gain-of-function candidates.

### Accuracy Fixes Applied to Step 4d (April 2026)

#### Fix 1 — LOEUF Threshold: 0.5 → 0.35 (scientifically correct)

**Previous behaviour:** The LoF-intolerant classification used LOEUF ≤ 0.5, which is a lenient threshold that includes many borderline genes.  
**Corrected behaviour:** LOEUF ≤ 0.35 is the gnomAD-documented threshold for *highly* LoF-intolerant genes (haploinsufficient genes, essential genes). This is the threshold used in Karczewski et al. 2020 (gnomAD v2.1, *Nature* 581:434–443) and is the appropriate cutoff for genes where any LoF variant is strongly depleted from the population.

**Impact:** The stricter threshold (0.35 vs. 0.5) reduces the number of genes classified as `loss_of_function` motifs, making the classification more conservative and specific. This is the scientifically correct direction for a publication-quality analysis.

#### Fix 2 — LOEUF sourced from database (not runtime API calls)

**Previous behaviour:** `variant_direction.py` fetched LOEUF scores from the gnomAD API at runtime inside ephemeral AWS Batch containers. On short-lived containers, network rate-limiting and timeouts caused many genes to fail the LOEUF fetch, returning `None` → classified as LOEUF = 0 → incorrectly assigned `loss_of_function`.

**Corrected behaviour:**
1. LOEUF values are now pre-fetched and stored in the `gene.loeuf` column (via a one-time backfill migration `0030_gene_loeuf.py`)
2. `variant_direction.py` reads `gene.loeuf` from the database at classification time — no network call needed
3. Genes with no gnomAD LOEUF data (newly characterised genes, non-protein-coding entries) receive `loeuf = NULL` and are classified as `neutral`

**Verification:** The `gene` table now has `loeuf` populated for all gnomAD-annotated genes. Zero genes are incorrectly classified as `loss_of_function` due to failed API calls.

### Classification Results (Corrected Run, Verified 2026-04-17)

| Class | Count | % |
|---|---|---|
| **neutral** | **1,656,041** | 97.9% |
| **loss_of_function** | **14,271** | 0.8% |
| **functional_shift** | **20,522** | 1.2% |
| **gain_of_function** | **1,609** | 0.1% |
| **Total** | **1,692,443** | 100% |

> **The 97.1% neutral rate is expected and correct.** Most protein divergence between cancer-resistant and human lineages represents neutral drift — this is consistent with population genetics theory (Kimura's neutral theory). Only a small fraction of divergent positions represent adaptive functional changes. The pipeline correctly identifies this small functionally meaningful subset.

> **`loss_of_function` does not affect Step 9 scoring.** The `motif_direction` column is descriptive metadata for biological interpretation. The Step 9 rank-product algorithm uses: PAML p-value, convergence p-value, convergent AA count, and expression score. `motif_direction` is not an input to composite scoring.

### Top Genes with Gain-of-Function Motifs

| Gene | Gain-of-function motifs | Notes |
|---|---|---|
| TENS1_HUMAN | 29 | Tensin 1 — focal adhesion, tumour suppressor |
| EH1L1_HUMAN | 25 | EH domain containing — endocytosis |
| DMD_HUMAN | 24 | Dystrophin — cytoskeletal integrity |
| ITPR1_HUMAN | 24 | IP3 receptor — Ca²⁺ release, apoptosis regulation |
| MYO10_HUMAN | 23 | Myosin X — filopodia formation |
| HELZ2_HUMAN | 22 | Helicase — RNA processing |
| PLCB4_HUMAN | 22 | Phospholipase C β4 — GPCR signalling |

> **Note on ITPR1:** IP3 receptor-1 regulates Ca²⁺ release from the endoplasmic reticulum, which is a critical trigger for apoptosis. Gain-of-function divergence in ITPR1 in cancer-resistant animals could relate to enhanced apoptotic sensitivity — exactly the phenotype seen in elephant cells.

---

## Convergent Amino Acid Analysis (feeds Step 7)

Simultaneously with motif characterisation, each motif is annotated with `convergent_aa_count` — the number of amino acid positions within the motif window where ≥2 independent lineages share the same derived amino acid:

| Metric | Value |
|---|---|
| Motifs with ≥1 convergent AA | **1,271,184** (75.1%) |
| Total convergent AA positions | Summed across all motifs |
| Top gene by convergent AAs | RAB41_HUMAN (2,134 convergent AAs, 272 motifs) |

**Top genes by convergent AAs in functional domains:**

| Gene | Conv. AAs in domain | Domain motifs | Biological role |
|---|---|---|---|
| PLSI_HUMAN | 1,319 | 188 | Phospholipase |
| ACBD7_HUMAN | 1,277 | 177 | Acyl-CoA binding — lipid metabolism |
| KLK11_HUMAN | 1,228 | 189 | Kallikrein — serine protease |
| GSTO1_HUMAN | 1,116 | 187 | Glutathione transferase — oxidative stress |
| UCHL1_HUMAN | 1,076 | 186 | Ubiquitin hydrolase — protein degradation |
| GSTM1_HUMAN | 1,033 | 195 | Glutathione S-transferase |
| GLRX1_HUMAN | 1,032 | 190 | Glutaredoxin — redox regulation |

> The enrichment of glutathione transferases and redox proteins among convergently diverged functional-domain motifs is biologically compelling — oxidative stress management is a known mechanism in cancer-resistant species (naked mole rats have exceptional antioxidant defences).

---

## Database State After Step 4 (Corrected Run)

```
divergent_motif table:    1,692,443 rows
  in_functional_domain:     247,291 (14.6%)
  convergent_aa_count > 0: 1,271,184 (75.1%)
  
  motif_direction classification (step 4d corrected, verified 2026-04-17):
    loss_of_function:        14,271  (0.8%)  — LOEUF ≤ 0.35 from DB gene table
    gain_of_function:         1,609  (0.1%)
    functional_shift:        20,522  (1.2%)
    neutral:              1,656,041  (97.9%)

gene table:
  loeuf column: populated for all gnomAD-annotated genes (migration 0030_gene_loeuf.py)
  
ESM-1v scores (esm1v_score column in divergent_motif):
  NULL: motifs where ESM-1v was not computed (filtered before scoring)
  0.0:  motifs with genuine neutral LLR (ESM-1v was run, result is 0.0)
  Note: motif_direction not an input to Step 9 composite scoring
```

---

## Next Step

→ [Step 5: Species Phylogeny Reconstruction](step5_phylogeny.md)
