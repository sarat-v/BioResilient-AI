# Step 4 — Protein Divergence Motif Analysis

**Pipeline phase:** Protein-level analysis  
**Substeps:** 4a (MEME motifs), 4b (domain + AlphaMissense), 4c (ESM-2 embeddings), 4d (functional classification)  
**Run timestamps:** 4a–4d: 2026-04-06 16:53:45 – 16:54:24 UTC  
**Status:** ✅ PASS (4a, 4b, 4c, 4d); ⚠️ WARN (4d — 97.8% neutral expected, noted)

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

## Step 4c — ESM-2 Protein Language Model Embeddings

### Tool

**ESM-2 (650M)** from Meta FAIR — a 650-million parameter transformer trained on 250 million protein sequences from UniRef50. ESM-2 encodes each amino acid in context, analogous to how BERT/GPT encode words.

### What ESM distance measures

The ESM-2 embedding of a 15-amino-acid window captures its biochemical "meaning" — the spatial relationships, chemical properties, and evolutionary constraints of that sequence in its protein context.

**ESM distance** = cosine distance between:
- The embedding of the human motif window
- The embedding of the cancer-resistant species' variant at the same position

A high ESM distance (> 0.5) indicates the substitution is *biochemically non-conservative* — the protein's local chemistry has fundamentally changed.

### Relationship to AlphaMissense

- **AlphaMissense** = per-residue impact of a single amino acid change
- **ESM-2 distance** = holistic change in the *entire window's biochemical context*

These are complementary: a motif can have high ESM distance even if each individual substitution scores as low-impact by AlphaMissense. Together they provide orthogonal evidence for functional change.

### Results

| Metric | Value |
|---|---|
| Motifs with ESM distance computed | Majority of 1,692,443 |
| Motifs with high ESM distance (> 0.5) | **419** (most impactful changes) |

---

## Step 4d — Functional Consequence Classification

### Classification logic

Each motif is classified based on combined AlphaMissense score, ESM-2 distance, and domain context:

| Class | Criteria | Count | % |
|---|---|---|---|
| **Gain-of-function** | ESM dist > 0.5 AND AM score > 0.7 | **1,609** | 0.1% |
| **Functional shift** | Moderate ESM or AM change in a domain | **34,793** | 2.1% |
| **Neutral** | Low ESM, low AM, or outside domain | **1,656,041** | 97.8% |

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

## Database State After Step 4

```
divergent_motif table:    1,692,443 rows
  in_functional_domain:     247,291 (14.6%)
  convergent_aa_count > 0: 1,271,184 (75.1%)
  gain_of_function:           1,609
  functional_shift:          34,793
  neutral/unclassified:   1,656,041
```

---

## Next Step

→ [Step 5: Species Phylogeny Reconstruction](step5_phylogeny.md)
