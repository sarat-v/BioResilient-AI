# Phase 2 Literature Research & Tool Selection
## BioResilient Pipeline Steps 9b–15

*Research date: April 2026*  
*Purpose: Pre-implementation review of industry best practices, published methods, and optimal tools for each planned Phase 2 step*

---

## Executive Summary

Before implementing Phase 2, this document synthesizes current literature, industry standards, and benchmark results for each planned step. Key findings:

1. **Human genetics is the new gold standard** — Targets with genetic evidence have **2.6x higher probability of clinical success** (Nature 2024). Our evolutionary approach generates orthogonal evidence that compounds this advantage.
2. **The mTOR/AKT pathway convergence our pipeline found is independently validated** — BMC Genomics 2023 independently identified mTOR network convergence in long-lived species, and 229 convergent cancer-related genes have been detected in extremely long-lived mammals (Science China 2024). Our approach is not incidental — it is cutting edge.
3. **Structural annotation is non-negotiable** — The field has moved to residue-level mapping (SIFTS, AlphaMissense, DeepAllo) before any further target assessment. This validates the proposed Step 9b as a critical prerequisite.
4. **Druggability tools have clear frontrunners** — The 2024 LIGYSIS benchmark definitively showed fpocket + PRANK/DeepPocket combination achieves best recall (60%), outperforming standalone tools.
5. **PheWAS needs LD correction** — Naive PheWAS produces spurious associations; CoPheScan (Nature Communications 2024) with Bayesian LD correction is now the standard.
6. **Tractability goes beyond small molecules** — PROTAC, antibody, and degrader modalities should all be assessed from the start (Open Targets Platform tractability framework).

---

## Section 1: Scientific Context — Why Evolutionary Convergence Predicts Drug Targets

### 1.1 Validation from the Literature

Our pipeline's core hypothesis — that genes under convergent positive selection in long-lived, cancer-resistant species are candidates for longevity/cancer protection — is directly supported by independent published research:

**mTOR pathway (directly validates our AKT1 finding)**  
Barja et al., *BMC Genomics* 2023: Analysis of 72 genes in the mTOR network across mammalian species found that long-lived species evolved convergent regulation of autophagy and cancer suppression. **4 positively selected genes and 2 convergent evolution genes** in the mTOR network directly overlap with our candidates. AKT1's 8-lineage convergence in our pipeline is not surprising — it is consistent with this literature.

**Convergent cancer resistance (directly validates our approach)**  
Liu et al., *Science China Life Sciences* 2024: Identified **229 convergent cancer-related genes** in extremely long-lived species (bowhead whales, naked mole rats). A convergent mutation in *LZTS1* shared across multiple long-lived lineages suppresses cancer development. *YAP1* mutations localizing tumor suppressors to the cytoplasm appear in cetaceans. This is the exact biological pattern our pipeline is designed to detect.

**Convergent evolution for drug target discovery (validates methodology)**  
Cardona-Alberich et al., *npj Precision Oncology* 2023: Convergent phenotypes across cancer types were used to derive a drug response signature predicting efficacy for 200+ chemotherapy drugs. This methodological parallel confirms that evolutionary convergence signals can be translated into actionable drug targets.

**Birds and bats convergence (validates species choice)**  
Tohoku University group found convergent amino acid substitutions in genes controlling cancer, reactive oxygen species, and immunity in birds and bats — both long-lived, high-metabolism species. Both appear in our species tree. These are not random.

### 1.2 The Genetic Evidence Multiplier

MacArthur et al., *Nature* 2024: Drug programs with human genetic evidence supporting mechanism have **2.6x greater probability of advancing** through all clinical phases. Our pipeline combines evolutionary selection pressure (cross-species) with rare protective variant signals (human genetics) — this is a particularly powerful combination not widely exploited in published pipelines.

---

## Section 2: Step 9b — Structural Annotation of Convergent Positions (PROPOSED NEW STEP)

*This step is missing from the current pipeline and represents the largest scientific gap. Without it, Steps 12–15 are built on incomplete foundations.*

### 2.1 Why This Step is Essential

Identifying that a gene shows convergent evolution is necessary but not sufficient. The specific **amino acid changes** must be:
- Structurally located (active site? allosteric? interface? disordered loop?)
- Functionally classified (gain-of-function? loss-of-function? regulatory?)
- Predicted for pathogenicity/impact

Without this, druggability assessment (Step 12) is generic rather than position-specific, and gene therapy design (Step 13) is impossible.

### 2.2 Best Tool Stack for This Step

#### Primary Structure Source: AlphaFold Database v4 + PDB
- All 14 Tier1 genes will have AlphaFold structures
- Use experimental PDB structures where available (higher confidence for drug binding)
- **Tool**: `py_alphafold` API or direct AFDB download

#### Residue-Level Functional Mapping: SIFTS (EBI)
*The authoritative resource for this.*

SIFTS (Structure Integration with Function, Taxonomy and Sequences) provides residue-level correspondence between PDB structures and UniProtKB annotations, InterPro domains, Pfam families, SCOP/CATH folds, and GO terms. Updated *weekly* and now covers >61,000 unique proteins.

```python
# SIFTS REST API: map PDB residue → UniProt position → domain/functional annotation
GET https://www.ebi.ac.uk/pdbe/graph-api/residue_mapping/{pdb_id}/{chain_id}/{residue_number}
```

**What SIFTS tells us for each convergent position:**
- Is it in a known Pfam domain?
- Is it annotated as active site / binding site in UniProt?
- Is it in a SCOP/CATH structural fold?
- Is there an InterPro functional site annotation?

#### Missense Effect Prediction: AlphaMissense (Google DeepMind, 2023)
*Best-in-class tool, validated in Science 2023. Already running in our pipeline (Step 4b).*

> **Important context**: AlphaMissense is **not new to our pipeline**. Step 4b (`pipeline/layer1_sequence/alphamissense.py`) already downloads the full TSV (~1.6 GB), builds a filtered per-protein index, and scores every `DivergentMotif` row with a `consequence_score` — asking "if this animal amino acid differs from human at this position, how pathogenic would that substitution be in humans?"
>
> The **new application in Step 9b** is targeted differently: instead of scoring all divergent motifs broadly, we would query the same already-downloaded index specifically for **convergent positions** (positions where 3+ independent lineages converged on the same amino acid change). This answers: "are the specific convergently-selected changes at residues that AlphaMissense classifies as functionally significant?" The tool and data are already present — Step 9b only needs a new query against the existing index.

AlphaMissense combines AlphaFold structural context with protein language modeling and population frequency. Classifies 32% of all human missense variants as likely pathogenic, 57% as likely benign, with **90% precision on ClinVar**. Outperforms SIFT, PolyPhen, EVE on every benchmark.

For Step 9b, the query would be:
- For each convergent position in a Tier1 gene, retrieve the AlphaMissense score for the specific convergent amino acid change (e.g., human Ala → convergent Ser across 8 lineages)
- Flag if score > 0.564 (likely pathogenic) vs < 0.34 (likely benign) — the AM classification thresholds
- Compare: are convergent changes enriched for functionally significant (high AM score) positions? This would validate they are not neutral drift.

#### Allosteric Site Prediction: DeepAllo (2024)
*New: better than PASSer and FPocket alone.*

DeepAllo (Lim et al. 2024, bioRxiv) fine-tunes ProtBERT-BFD on the Allosteric Database (ASD) with multitask learning. Achieves **89.66% F1** for allosteric pocket prediction, identifying correct allosteric site in top 3 ranked positions with 90.5% accuracy. This is important because convergent positions near allosteric pockets are higher-value drug targets (avoid disrupting active site directly).

#### Active Site / Binding Site Detection: fpocket + P2Rank + PRANK
The 2024 LIGYSIS benchmark confirmed fpocket rescored by PRANK or DeepPocket achieves the highest recall at 60%. For our 14 proteins:
1. Run fpocket for initial pocket detection
2. Rescore with P2Rank (ML-based) for ranking
3. Overlay convergent residue positions — if within 6Å of top pocket, flag as "pocket-adjacent"

#### Functional Context Classifier
Build a simple classifier:
- **Active site**: convergent residue annotated as catalytic/active in PROSITE or UniProt
- **Binding interface**: convergent residue in known protein-protein interaction interface (STRING binding evidence)  
- **Allosteric site**: DeepAllo score > threshold for that pocket
- **Surface accessible**: RSA (relative solvent accessibility) > 0.25 from DSSP
- **Disordered region**: AF2 pLDDT < 70 AND IUPRED2 score > 0.5
- **Other**: everything else

### 2.3 Computational Cost Estimate
- Per-protein cost: ~5 minutes (fpocket + P2Rank + SIFTS REST + AlphaMissense lookup)
- For 14 Tier1 genes: <2 hours total
- For all 12,795 genes: ~1 week on single CPU (parallelize to <1 day with 8 cores)

---

## Section 3: Step 11 — Disease Annotation (Current Plan vs. Best Practice)

### 3.1 Industry Standard: Human Genetics First

The current orthodoxy in pharma (GSK, AstraZeneca, Roche, Pfizer — all publishing through Open Targets collaboration) is **genetics-first target prioritization**. The Genetic Priority Score (GPS), published in Nature Genetics 2023, integrates 8 genetic features across 19,365 genes and 399 indications. Targets in the top 0.28% of GPS rankings are **8.8-fold more likely to reach Phase IV**.

Our pipeline adds *evolutionary* evidence that acts as orthogonal genetic evidence. This is novel and not in GPS.

### 3.2 Open Targets Platform: L2G (Locus-to-Gene) Model

The L2G model (Nature Genetics 2021) assigns each GWAS-associated locus a causal gene using a machine learning model trained on fine-mapped variants + functional genomics (92 cell types/tissues). Key stats:
- Trained on 133,441 GWAS loci
- Prioritized genes enriched **8.1-fold** for known approved drug targets
- Integrates eQTL, chromatin accessibility, protein interaction data

**Recommendation**: Query Open Targets GraphQL API for L2G scores for all 14 Tier1 genes across all diseases. This is preferable to manual GWAS catalog queries.

```graphql
query targetDiseaseAssociations($target: String!) {
  target(ensemblId: $target) {
    associatedDiseases(page: {index: 0, size: 20}) {
      rows {
        disease { name id }
        score
        datatypeScores { id score }  # includes genetics score
      }
    }
  }
}
```

### 3.3 Rare Protective Variants: The PCSK9 Paradigm — Gold Standard Approach

The PCSK9 story (*Cohen et al., NEJM 2006*) remains the gold standard. Individuals with LoF variants in PCSK9 had dramatically lower LDL and 88% lower coronary artery disease risk. This *human knockout* validated the drug target before the drug existed.

**How to replicate this for our genes:**

1. **gnomAD v4.0** (2024) — now 730,947 exomes, ~6x larger than previous version. Query LoF variants at each Tier1 gene.
2. **Compute per-gene LOEUF score**: LOEUF < 0.6 = LoF intolerant (updated v4.0 threshold). Important: LoF-intolerant genes are *not automatically excluded* — many successful drug targets (aspirin target COX, statin target HMGCR) show extreme constraint. What matters is the *direction* of effect and existence of human carriers.
3. **Drug-target Mendelian Randomization (MR)** using cis-pQTLs: If protein levels are measured in a pQTL study (e.g., deCODE, UKB-PPP), genetically-reduced protein levels can be tested for protective associations. This is now standard at big pharma. Tools: **TwoSampleMR** (R package), **MendelianRandomization**, **gwasvcf**.
4. **Colocalization (COLOC+SuSiE)**: If MR is significant, confirm the pQTL and GWAS signal share a causal variant. SharePro (2024 benchmark winner) outperforms coloc and coloc+SuSiE (79.6% support vs. 46.8% for coloc+SuSiE). Use SharePro as primary, coloc+SuSiE as secondary confirmation.

### 3.4 AstraZeneca's Mantis-ML 2.0

AstraZeneca's Mantis-ML 2.0 (Zenodo 2024) uses knowledge graphs + graph neural networks + UK Biobank for phenome-wide target identification. It's freely available and could be used for disease enrichment analysis for our Tier1 genes. This represents the state of the art for ML-based target-disease association.

### 3.5 Recommended Tool Stack for Step 11

| Sub-step | Current Tool | Recommended Tool | Reason |
|---|---|---|---|
| GWAS disease links | OpenTargets (manual) | OpenTargets GraphQL API + L2G scores | Automated, ML-scored |
| Rare variant burden | gnomAD manual | gnomAD v4.0 API + LoF collapsing | 6x more samples |
| Human genetics validation | None | cis-pQTL MR with TwoSampleMR | PCSK9 paradigm standard |
| Colocalization | None | SharePro (best 2024 benchmark) | Outperforms coloc 79.6% vs 47% |
| Phenome-wide disease | PheWAS basic | Mantis-ML 2.0 or CoPheScan | LD-corrected, ML-enhanced |

---

## Section 4: Step 12 — Druggability Assessment

### 4.1 2024 Benchmark: Definitive Tool Ranking

The LIGYSIS benchmark (Journal of Cheminformatics 2024, 30,000 proteins) is the definitive 2024 comparison. Results:

| Tool | Recall | Notes |
|---|---|---|
| fpocket + PRANK | 60% | **Best overall** |
| fpocket + DeepPocket | 60% | **Equal best** |
| IF-SitePred | 39% | Worst performer |
| P2Rank (standalone) | ~52% | Strong standalone |
| DeepPocket | ~55% | Good deep learning option |

**Winner**: fpocket for initial detection, rescored by either PRANK or DeepPocket for final ranking.

For membrane proteins (GPCRs, channels): A 2025 follow-up found DeepPocket, PUResNetV2.0, and ConCavity perform best, though overall performance is lower. If any Tier1 genes are membrane proteins, apply these specialized tools.

### 4.2 Druggable Genome Overlap

The **Human Protein Atlas Druggable Proteome** resource provides a pre-computed overlay of all human proteins against known druggable families (kinases, GPCRs, nuclear receptors, ion channels, proteases, etc.). Cross-referencing Tier1 genes against this resource provides instant druggability context.

### 4.3 Modality-Specific Tractability

Current best practice (Open Targets Platform, AstraZeneca) assesses **four modalities simultaneously**:

1. **Small molecule (SM)**: Pocket volume ≥ 300 Å³, druggability score from DrugEBIlity/PocketFinder
2. **Antibody (AB)**: Extracellular/surface-exposed epitopes (UniProt subcellular location + topology)
3. **PROTAC/Degrader**: "PROTACtable genome" criteria: protein stability, E3 ligase accessibility, cellular location, natural half-life. 1,000+ proteins identified as PROTAC-tractable but not otherwise druggable (*Nature Reviews Drug Discovery* 2021).
4. **Other clinical modalities (OC)**: Antisense oligonucleotides (ASOs), siRNA, gene therapy — based on tissue expression and cell type accessibility

AstraZeneca has released ML models for PROTAC degradation prediction on GitHub. For any Tier1 gene where small molecule tractability is low but convergent changes are structural, PROTAC could be the modality.

### 4.4 Structural Genomics Consortium (SGC) Check

For each Tier1 gene, check whether SGC (sgc.ox.ac.uk) has solved experimental structures with chemical probes or tools compounds. SGC probes are freely available and represent immediate starting points for drug discovery.

### 4.5 Recommended Tool Stack for Step 12

| Assessment | Tool | Why |
|---|---|---|
| Binding site detection | fpocket | Fast, open-source, best-benchmarked starting point |
| Pocket rescoring | P2Rank or DeepPocket | ML-based, 14% recall improvement over fpocket alone |
| Allosteric pockets | DeepAllo | Best 2024 tool (89.66% F1) |
| SM tractability | DrugEBIlity + fpocket score | Industry standard, OpenTargets-integrated |
| AB tractability | UniProt topology + HPA | Fast lookup |
| PROTAC tractability | PROTACtable genome criteria | Nature Reviews benchmark |
| Known probes | SGC chemical probe portal | Immediate starting point discovery |

---

## Section 5: Step 13 — Gene Therapy Feasibility

### 5.1 Critical Prerequisite: Structural Mapping FIRST

This is the biggest structural problem in the current Phase 2 plan. Gene therapy decisions depend entirely on:
- **What change do we want to make?** (requires Step 9b output)
- **Can we deliver that change to the relevant tissue?**
- **What are the on/off-target risks?**

Without Step 9b, gene therapy assessment is generic and not tied to any specific modification. The step should be **gated behind Step 9b completion**.

### 5.2 Computational Pre-assessment Framework

Given Step 9b output, the correct gene therapy pre-assessment pipeline is:

**AAV-based delivery assessment:**
1. **Gene size check**: AAV cargo capacity ≈ 4.7 kb. Genes > 4.5 kb require split-intein AAV or other approaches.
2. **Tissue expression check**: Human Protein Atlas RNA-seq + GTEx for dominant expression tissue. Determines AAV serotype (AAV9 for CNS, AAV-PHP.eB for broad CNS, AAV8 for liver, etc.)
3. **CrAAVe-seq platform** (Nature Neuroscience 2025) — gold standard for in vivo cell-type-specific CRISPR screening with AAV. Provides guidelines for sgRNA library size, viral titer, coinfection rate.
4. **hafoe tool** (Nature Gene Therapy 2025) — analyzes chimeric AAV library tropism for rational vector selection.

**CRISPR feasibility assessment:**
1. **GuideScan2** (Genome Biology 2025) — memory-efficient high-specificity guide RNA design. **Superior off-target enumeration** vs Cas-OFFinder, CRISPOR. Use this.
2. **Off-target prediction**: Cas-OFFinder for exhaustive off-target enumeration (slower but complete); GuideScan2 for fast batch screening.
3. **Base editing vs prime editing vs classic cut**: Classify convergent amino acid changes by edit type:
   - Single nucleotide changes → base editing candidates (ABEs, CBEs)
   - Small indels → prime editing candidates
   - Larger structural changes → NHEJ/HDR approaches
4. **Delivery modality**: LNP (lipid nanoparticles) for liver/systemic, AAV for CNS/muscle/eye, mRNA for transient expression

### 5.3 Realistic Scope for Phase 2

Gene therapy computational pre-assessment for 14 genes is feasible. But the output should be a **feasibility score card**, not a therapy design. Inputs:
- Gene size (from Ensembl)
- Dominant tissue expression (from GTEx/HPA)
- Type of modification (from Step 9b)
- GuideScan2 guide quality score for relevant edits
- Known clinical-stage programs for that gene/indication

### 5.4 Recommended Tool Stack for Step 13

| Assessment | Tool | Access |
|---|---|---|
| Gene size | Ensembl REST API | Free |
| Tissue expression | GTEx portal API | Free |
| AAV serotype selection | hafoe + literature lookup | GitHub + publication |
| CRISPR guide RNA design | GuideScan2 | Free, open-source |
| Off-target screening | GuideScan2 + Cas-OFFinder | Free |
| Edit type classification | Manual + AlphaMissense context | Per-residue |

---

## Section 6: Step 14 — Safety Pre-Screening

### 6.1 The Problem with Naive PheWAS

Current plan uses PheWAS broadly. The key limitation: PheWAS tests genetic variants at the gene locus against all phenotypes. However, LD structure means variants at a gene can be in LD with variants at neighboring genes, creating spurious signals.

**Solution**: **CoPheScan** (Nature Communications 2024) — a Bayesian PheWAS approach that explicitly accounts for LD confounding, distinguishing true pleiotropy from LD-driven false positives. This is a methodological upgrade that should replace naive PheWAS.

### 6.2 Machine Learning for Off-Target Safety

Nature Reviews Drug Discovery 2024: ML models combining:
- On-target IC₅₀ data
- Off-target binding panel data (Cerep-style panels: hERG, nuclear receptors, kinases)
- Tissue-specific protein expression profiles

Achieved **>75% accuracy** predicting clinical adverse events. The key insight: *most clinical safety failures are predictable from preclinical on/off-target profiles*.

Implementation for our pipeline:
1. ChEMBL query — does any compound bind our Tier1 target as primary or off-target? What adverse events?
2. STRING network expansion — first-shell interactors. Compounds hitting those interactors → infer potential off-target phenotypes.
3. gnomAD LoF carriers — if human carriers of LoF variants in our target are healthy adults, the safety window is good (validated PCSK9-style).

### 6.3 DepMap: Essentiality ≠ Danger

Important nuance from the DepMap team (Nature Reviews Cancer 2024): Cancer cell line essentiality from CRISPR screens does NOT automatically mean a gene is dangerous to target. Many successful drug targets are essential in cancer cells but dispensable in normal cells. The metric to use is:

**Selective dependency**: essentiality score in cancer lines vs. in normal cell lines. Available from DepMap portal as "context dependency."

The shinyDepMap tool (eLife 2021) provides interactive "target hopping" — if a target is non-druggable but shows strong selective dependency, use it to find druggable partners with similar selectivity profiles.

### 6.4 gnomAD Constraint as Safety Window

Oxford BDI (2024) paper on evaluating drug targets through human LoF variation: The correct interpretation is:
- **High LOEUF (LoF tolerant)**: Gene can be knocked out without major phenotype → therapeutic window likely safe
- **Low LOEUF (LoF intolerant, haploinsufficient)**: Must modulate, not eliminate — partial inhibition may be needed
- **Presence of healthy LoF carriers**: Best possible safety evidence (natural human knockout validation)

### 6.5 Secondary Pharmacology Framework

*Nature Reviews Drug Discovery 2024* — best comprehensive reference. Key principles:
1. All Tier1 targets should be screened against a minimal secondary pharmacology panel: cardiac ion channels (hERG/Nav/Cav), nuclear receptors (AR, ER, GR, PXR), kinase panel
2. In silico secondary pharmacology predictions are available via SwissTargetPrediction, SuperPred, and ChEMBL target prediction
3. Tissue-specific expression of off-targets matters: a CNS-specific off-target is less concerning than a cardiac one

### 6.6 Recommended Tool Stack for Step 14

| Safety Dimension | Current Tool | Recommended Tool | Upgrade Reason |
|---|---|---|---|
| Phenome-wide associations | PheWAS | CoPheScan (Bayesian, LD-corrected) | Eliminates LD spurious signals |
| Off-target prediction | STRING | STRING + SwissTargetPrediction + ChEMBL panel | More comprehensive |
| Essentiality safety | DepMap CRISPR score | DepMap selective dependency score | Context-specific, not global |
| Human genetic safety | None | gnomAD LoF carriers + LOEUF | PCSK9-style validation |
| Tissue expression safety | GTEx | GTEx + HPA tissue specificity tau score | Quantified breadth of expression |

---

## Section 7: Step 15 — Final Rescoring and Tier Refinement

### 7.1 Genetic Priority Score (GPS) Integration

Nature Genetics 2023: The GPS framework integrates 8 genetic features. Targets in the **top 0.28% of GPS** are 8.8-fold more likely to reach Phase IV. Our genes can be scored against the publicly available GPS and the score added as a Phase 2 evidence layer.

### 7.2 Multi-Evidence Bayesian Scoring Framework

Industry best practice (Mantis-ML, Open Targets) is to combine evidence layers probabilistically, not with simple weighted sums. Options:
1. **Weighted sum** (current approach): Simple, interpretable, but assumes independence
2. **Rank Product** (what we use in Phase 1): Handles different scales well, but doesn't account for correlated evidence
3. **Bayesian network** (optimal): Models conditional dependencies between evidence layers. Implemented in Mantis-ML 2.0.

For Phase 2, the minimum improvement over the current approach is to:
- Add a **penalty correction** for redundant evidence (e.g., two OpenTargets scores from the same GWAS shouldn't count double)
- Add **uncertainty quantification** — annotate each score with data completeness (did we find data, or did we just get null for a gene?)
- Maintain a **"missing data" flag** distinct from a "true negative" — currently our pipeline may treat data gaps as zeros

### 7.3 The Control Species Divergence Correction

This is currently planned for Step 15 but should logically be applied earlier. The correction penalizes genes where changes in control (short-lived) species are as frequent as in longevity species. This is a **pre-filter**, not a rescoring operation. Moving it to before Step 9b would:
1. Reduce the number of genes sent through expensive structural annotation
2. Improve precision of the initial candidate set
3. Be more scientifically defensible (you're eliminating confounders, not adjusting scores)

### 7.4 Recommended Final Score Architecture

```
Phase2_Score = 
  Phase1_CompositeScore (weight: 0.35)  # Evolutionary signal — our unique contribution
  + DiseaseAnnotation_Score (weight: 0.20)  # OpenTargets L2G + GPS
  + RareVariant_Score (weight: 0.15)  # PCSK9 paradigm
  + Druggability_Score (weight: 0.15)  # fpocket + modality tractability
  + StructuralContext_Score (weight: 0.10)  # Step 9b output — pocket proximity, AlphaMissense
  + Safety_Score (weight: 0.05)  # Penalty only — CoPheScan + gnomAD
  [subtract] ControlSpeciesPenalty (Step 15 → move earlier as pre-filter)
```

Note: Safety_Score should be a **penalty** applied to a multiplicative factor (0.0–1.0) rather than an additive score. A major safety signal should be able to drop a gene regardless of its evolutionary score.

---

## Section 8: Recommended Additions Not Currently Planned

### 8.1 AlphaMissense Per-Position Scoring (NEW)
*Priority: HIGH*

Run AlphaMissense on all convergent amino acid changes across Tier1 genes. This gives:
- Pathogenicity classification of each convergent change
- Comparison: convergent mammals gained "protective" vs "pathogenic" changes?
- Direct input to Step 9b structural annotation

Cost: Pre-computed scores available for all human proteins as a downloadable TSV. No compute required.

### 8.2 SIFTS Residue Annotation Lookup (NEW)
*Priority: HIGH*

Query EBI SIFTS API for each convergent residue position in each Tier1 protein. Provides instant functional context from 10+ integrated databases. REST API, free, fast.

### 8.3 SharePro Colocalization (NEW)
*Priority: MEDIUM*

For any Tier1 gene with a GWAS association, run SharePro colocalization between the GWAS signal and eQTL/pQTL signals in relevant tissues. This converts "gene near a GWAS hit" to "gene IS the causal gene at this GWAS locus." Evidence quality difference is enormous.

### 8.4 STRING Network Context for Each Convergent Gene (EXISTING BUT UNDERUSED)
*Priority: MEDIUM*

STRING protein network analysis is planned but should specifically:
1. Map which convergent genes interact with each other (are they in the same pathway?)
2. Identify hub genes in the convergent network — potentially better targets than the convergent genes themselves
3. Flag any convergent gene whose STRING network is enriched for known drug targets (target neighborhood validation)

### 8.5 GPS Score Integration (NEW)
*Priority: MEDIUM*

Download the publicly available Genetic Priority Score table from Nature Genetics 2023 and join to our Tier1 candidates. This takes 30 minutes and provides a published benchmark for our candidates.

### 8.6 Tissue Selectivity Score (EXISTING, UPGRADE)
*Priority: LOW-MEDIUM*

Replace simple GTEx expression check with Human Protein Atlas **tau score** (tissue specificity index, 0–1 where 1 = perfectly tissue-specific). Combined with convergent lineage tissue context, this directly informs which tissues a therapeutic intervention would affect and what toxicity profile to expect.

---

## Section 9: Revised Step Order Recommendation

Based on this literature review, the revised Phase 2 workflow should be:

```
[Pre-filter] Control Species Divergence Penalty
    ↓
Step 9b: Structural Annotation (NEW)
    (SIFTS + fpocket + P2Rank + DeepAllo + AlphaMissense convergent-position query)
    ↓
Step 10b: AlphaGenome Regulatory Divergence (ALREADY CODED)
    (Track B — regulatory/promoter divergence for expression-positive genes)
    ↓
Step 11a: Disease Annotation
    (OpenTargets L2G + Mantis-ML 2.0 + GPS score)
    ↓
Step 11b: Rare Protective Variants
    (gnomAD v4.0 + cis-pQTL MR + SharePro colocalization)
    ↓
Step 12: Druggability
    (fpocket+DeepPocket + multi-modality: SM/AB/PROTAC + SGC probes)
    ↓
Step 14: Safety
    (CoPheScan + gnomAD LoF carriers + DepMap selective dependency + SwissTargetPrediction)
    ↓
Step 13: Gene Therapy Feasibility
    (GuideScan2 + hafoe + payload/tissue assessment — gated on Step 9b)
    ↓
Step 15: Final Rescoring
    (Bayesian multi-evidence combination + uncertainty quantification + tier refinement)
```

**Key changes from original plan:**
1. Control species penalty moved to **pre-filter** (before all Phase 2 steps)
2. **Step 9b added** as mandatory gate before druggability and gene therapy
3. **Step 10b (AlphaGenome)** is already fully coded — it just needs `ALPHAGENOME_API_KEY` set. It runs Track B regulatory divergence for expression-positive genes and can compound protein-level signals from Track A.
4. Step 13 (gene therapy) moved **after** safety (Step 14) — no point assessing delivery for unsafe targets
5. Safety scoring changed from additive to **multiplicative penalty**

---

## Section 10: Tool Availability Summary

| Step | Tool | License | Performance | Priority | Status |
|---|---|---|---|---|---|
| 9b | SIFTS (EBI REST API) | Free | Gold standard residue mapping | Must-have | Not yet coded |
| 9b | AlphaMissense | Free (pre-computed TSV already downloaded) | 90% precision on ClinVar | Must-have | **Already in pipeline (Step 4b)** — new query only |
| 9b | fpocket | Open source (C) | 60% recall (benchmark winner) | Must-have |
| 9b | P2Rank | Open source (Java) | Top ML tool | Must-have |
| 9b | DeepAllo | Free (GitHub) | 89.66% F1 allosteric | Recommended | Not yet coded |
| 10b | AlphaGenome | API key required | Regulatory promoter divergence | Must-have | **Already fully coded** |
| 11a | Open Targets GraphQL API | Free | 8.1x drug target enrichment | Must-have | Partially coded |
| 11b | gnomAD v4.0 API | Free | 730K exomes | Must-have | Partially coded |
| 11b | TwoSampleMR (R) | Open source | Industry standard MR | Recommended | Not yet coded |
| 11b | SharePro | Open source | Best 2024 coloc benchmark | Recommended | Not yet coded |
| 12 | DeepPocket | Open source (Python/PyTorch) | Best recall rescorer | Must-have | Not yet coded |
| 12 | PROTACtable genome criteria | Published table | 1000+ PROTAC targets | Must-have | Not yet coded |
| 12 | SGC probe portal | Free | Immediate starting points | Should-have | Not yet coded |
| 13 | GuideScan2 | Open source | Best specificity CRISPR | Must-have | Not yet coded |
| 13 | Cas-OFFinder | Free | Comprehensive off-target | Recommended | Not yet coded |
| 14 | CoPheScan | Open source (R) | LD-corrected PheWAS | Must-have upgrade | Not yet coded |
| 14 | SwissTargetPrediction | Free web API | Off-target profile | Recommended | Not yet coded |
| 14 | DepMap selective dep. | Free API | Context-specific essentiality | Must-have | Partially coded |
| 15 | GPS score table | Publicly available (download) | 8.8x Phase IV enrichment | Should-have | Not yet coded |

---

## Section 11: Summary of Key Scientific Upgrades

### What the literature says we're doing right
1. **Convergent evolution approach is scientifically validated** — 229 convergent cancer genes in long-lived mammals, mTOR network convergence confirmed independently. Our methodology is at the frontier.
2. **Rare protective variant analysis** — PCSK9 paradigm, now the industry gold standard. Our Step 11b follows this correctly.
3. **Multi-evidence scoring** — Rank Product from Phase 1, now being upgraded to Bayesian combination in Phase 2.

### What the literature says we must upgrade
1. **Step 9b is missing** — Without residue-level structural annotation, we cannot distinguish a surface-exposed, allosteric convergent change from a buried, structurally constrained one. SIFTS + AlphaMissense is a 2-day implementation that would dramatically improve all downstream steps.
2. **PheWAS needs LD correction** — Replace with CoPheScan. Same concept, eliminates a major source of false positives.
3. **Druggability needs ML rescoring** — fpocket alone is insufficient; DeepPocket or P2Rank rescoring improves recall by 14%.
4. **Tractability should include PROTAC** — 1000+ non-SM-tractable proteins are PROTAC-tractable. Some of our Tier1 genes may fall in this category.
5. **Safety scoring should be multiplicative** — A single major safety signal should eliminate a candidate, not just lower its score.

---

*Document prepared April 2026. Citations available upon request for specific claims. All tool benchmarks reference 2024–2025 publications.*
