# BioResilient Pipeline — Step-by-Step Reports

**Phenotype:** Cancer Resistance | **Pipeline version:** v4 | **Date:** April 2026

---

## Phase 1 Step Reports (Evolutionary Genomics)

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

## Phase 2 Step Reports (Clinical Translation)

| Step | Title | Status | Key Output |
|---|---|---|---|
| [Step 9b](step9b_structural_annotation.md) | Structural Context Annotation | ✅ PASS | 2,092 motifs; 69 pocket-adjacent; AlphaMissense + fpocket |
| [Step 10b](step10b_regulatory_divergence.md) | Regulatory Divergence (AlphaGenome) | ✅ PASS | 202 records; 45 genes with regulatory signal |
| [Step 11](step11_disease_annotation.md) | Disease Annotation (OpenTargets + GWAS + gnomAD) | ✅ PASS | 119 genes with OT scores; 140 with GWAS; 22 with known drugs |
| [Step 12](step12_druggability.md) | Druggability Assessment (fpocket + P2Rank) | ✅ PASS | 279/280 genes with pockets; 76 PROTAC-tractable |
| [Step 13](step13_gene_therapy.md) | Gene Therapy Feasibility | ✅ PASS (informational) | 15 genes evaluated; all AAV-compatible |
| [Step 14](step14_safety_screen.md) | Safety Screen (DepMap + GTEx + PheWAS) | ✅ PASS | 2 essential genes; safety floor enforced |
| [Step 15](step15_final_scoring.md) | Phase 2 Final Scoring & Tier Assignment | ✅ PASS | **74 Validated, 81 Tier 2 — 155 total candidates** |

---

## Phase 3 Step Reports (Translational Readiness)

| Step | Title | Status | Key Output |
|---|---|---|---|
| Step 16 | Translational Status (ClinicalTrials.gov + ChEMBL) | ⏳ PENDING run | Known drug phases; clinical trial mapping |
| Step 17 | Preclinical Readiness Scoring | ⏳ PENDING run | Synthesis score, model score, assay score, readiness tier |
| Step 18 | Target Dossier Generation (Markdown) | ⏳ PENDING run | Per-gene dossier: all evidence + recommended first experiments |

> **Pipeline boundary:** Step 18 (target dossier) is the **final computational output**. The pipeline delivers a ranked, annotated list of therapeutic targets with supporting evidence from 18 integrated computational steps. Wet-lab validation (binding assays, cellular phenotype, in-vivo studies) is the responsibility of the downstream experimental programme.

> **Out of scope:** Steps 19–25 (virtual screening, ADMET, lead optimisation) require wet-lab binding data as input and constitute a separate product. Contact the BioResilient team to initiate the Lead Discovery programme.

---

## Phase 2 Top-Line Results

- **1 Tier 1 gene:** G6PD (composite 0.709) — only gene clearing the 0.70 threshold
- **74 Validated genes:** Tier1/Tier2 with corroborating human genetic evidence (GWAS / OT / protective variants)
- **81 Tier 2 genes:** Strong evolutionary signal, awaiting human genetic corroboration
- **G6PD independently rediscovered** — validates the cross-species evolutionary approach
- **3 genes with Phase 3 drugs at convergent sites:** CD80 (Galiximab), KCNS1 (Guanidine), XDH (Allopurinol)

---

## Comprehensive Reports

**Phase 1:**  
→ [BioResilient_Pipeline_Report_Steps1_9.md](../BioResilient_Pipeline_Report_Steps1_9.md)  
→ [BioResilient_Pipeline_Report_Steps1_9.docx](../BioResilient_Pipeline_Report_Steps1_9.docx)

**Phase 2:**  
→ [Phase2_Pipeline_Report.md](../Phase2_Pipeline_Report.md)  
→ [Phase2_Pipeline_Report.docx](../Phase2_Pipeline_Report.docx)
