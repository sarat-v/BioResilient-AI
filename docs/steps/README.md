# BioResilient Pipeline — Step-by-Step Reports

**Phenotype:** Cancer Resistance | **Pipeline version:** v4 | **Date:** April 2026

---

## Step Reports

| Step | Title | Status | Key Output |
|---|---|---|---|
| [Step 1](step1_environment_setup.md) | Environment Setup & Validation | ✅ PASS | Infrastructure verified |
| [Step 2](step2_species_selection.md) | Species Selection & Proteome Assembly | ✅ PASS | 18 species, 16 cancer-resistant + 2 controls |
| [Step 3](step3_ortholog_detection.md) | Ortholog Detection & Genomic Regions | ✅ PASS | 12,795 genes, 200,979 orthologs |
| [Step 4](step4_motif_analysis.md) | Protein Divergence Motif Analysis | ✅ PASS | 1,692,443 motifs, 1,271,184 convergent |
| [Step 5](step5_phylogeny.md) | Species Phylogeny Reconstruction | ⚠️ WARN | Tree built; bootstrap mean 43.9 (mammalian backbone reliable) |
| [Step 6](step6_positive_selection.md) | Positive Selection Analysis (PAML) | ✅ PASS | 495 genes p<0.05; 7,102 with selection signal |
| [Step 7](step7_convergent_evolution.md) | Convergent Evolution Detection | ✅ PASS | 6,601 genes in ≥3 lineages; 782 in all 5 |
| [Step 8](step8_functional_evidence.md) | Functional Evidence Gathering | ✅ PASS | 16,479 evidence records, 5,775 genes |
| [Step 9](step9_composite_scoring.md) | Composite Scoring & Tier Assignment | ✅ PASS | **36 Tier 1, 96 Tier 2 candidates** |

---

## Top-Line Results

- **36 Tier 1 candidates** (FDR < 5%): highest confidence genes for experimental follow-up
- **96 Tier 2 candidates** (FDR 5–20%): secondary priority list
- **Best-in-class benchmark validated:** CDK4 (known cancer driver) appears Tier 1 rank 12 ✅

---

## Comprehensive Report

For the full integrated analysis with cross-step interpretation:  
→ [BioResilient_Pipeline_Report_Steps1_9.md](../BioResilient_Pipeline_Report_Steps1_9.md)
