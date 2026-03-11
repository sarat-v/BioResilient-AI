# BioResilient AI — Platform Updates Since the Researcher Brief

> **Brief version**: Confidential Partner Brief · Phase 1 System · March 2026  
> **This document covers**: All additions, fixes, and enhancements made after the brief was written.

---

## Summary of Changes

The brief described the *intended* architecture. This document records what was built, what was corrected, and what was added beyond the spec. It is organised by category and includes the rationale for each change.

---

## 1. Critical Bug Fixes

### 1.1 Twelve Pipeline Steps Were Not Executing

**Problem**: The pipeline orchestrator defined all 30 steps and listed them in its step registry, but the main dispatch loop (`run_pipeline()`) was missing `elif` clauses for 12 sub-steps. Every time the pipeline ran, these steps were silently marked "complete" without executing any logic.

**Affected steps** (all now fixed):

| Step | What It Does |
|------|-------------|
| step4b | Pfam domain annotation + AlphaMissense pathogenicity scores |
| step4c | ESM-1v structural variant effect scoring |
| step4d | Variant direction inference (GoF / LoF / neutral) |
| step6b | FEL + BUSTED supplementary selection tests |
| step6c | RELAX branch-specific selection acceleration |
| step7b | True convergent amino acid substitution detection |
| step8b | Bgee cross-species expression supplement |
| step11b | Rare protective variant mapping (gnomAD / PCSK9 paradigm) |
| step11c | Literature validation (PubMed sanity check) |
| step11d | Pathway convergence enrichment (hypergeometric test) |
| step12b | P2Rank ML-based pocket prediction |
| step14b | DepMap CRISPR essentiality + GTEx expression breadth |

**Impact**: All motif-level annotations (domain, AlphaMissense, ESM-1v, variant direction), supplementary selection tests (FEL, BUSTED, RELAX), and Phase 2 safety/druggability data were completely absent from every prior run.

### 1.2 RELAX Step Reading Wrong Pickle File

**Problem**: `step6c_relax()` loaded the alignment cache with `pickle.load(f)` directly instead of unpacking the dict keys `"aligned"` and `"motifs_by_og"`. It also searched for a separate `motifs_by_og.pkl` file that never existed (step4 writes both into a single `aligned_orthogroups.pkl`). Additionally it looked for the species treefile at the wrong path (`data/species.treefile` instead of `data/phylo/species.treefile`).

**Fix**: Corrected all three — proper dict unpacking, removed the phantom second pickle file, corrected tree path.

### 1.3 gnomAD GraphQL Query Used Wrong Field Names

**Problem**: The query used `LOEUF: gnomad_constraint { loe_uf }` — a field alias syntax and field name (`loe_uf`) that do not exist in gnomAD v4's GraphQL schema. Every variant direction call silently returned `None` for LOEUF, causing most motifs to default to `lof_tolerant = True` (the fallback when LOEUF is unknown).

**Fix**: Updated to gnomAD v4 field names: `gnomad_constraint { pLI, oe_lof_upper }`.

---

## 2. Species Registry Corrections

### 2.1 Rougheye Rockfish ID and Lineage Corrected

**Problem**: The species ID was `"bowhead_rockfish"` (copy-paste from the bowhead whale entry) and the lineage group was `"Sharks"`. *Sebastes aleutianus* is a ray-finned fish (Actinopterygii), not a shark (Chondrichthyes). Misclassifying it as "Sharks" caused the rockfish and Greenland shark to count as the same lineage, reducing the available convergence power and under-weighting their combined signal.

**Fix**:
- ID: `"bowhead_rockfish"` → `"rougheye_rockfish"`
- Lineage group: `"Sharks"` → `"Fishes"`
- All hardcoded references updated in `convergence.py` and `alphagenome.py`

### 2.2 "Fishes" Lineage Added to Phylogenetic Distance Table

The `LINEAGE_DIVERGENCE_MY` divergence-time lookup table (28 pairs) did not include "Fishes" as a separate lineage from "Sharks". Seven new pairs were added with TimeTree.org-sourced estimates (e.g., Fishes–Rodents: 430 My; Fishes–Sharks: 420 My).

**Result**: The platform now correctly represents 8 independent resilient lineages as specified in the brief (Rodents, Cetaceans, Bats, Sharks, **Fishes**, Birds, Proboscideans, Salamanders).

### 2.3 AlphaGenome Lineage Map Was Incomplete

The regulatory divergence module (`layer_regulatory/alphagenome.py`) had its own hardcoded species-to-lineage map that was missing 6 species: blind mole rat, right whale, Brandt's bat, amazon parrot, budgerigar, human. All were added.

---

## 3. New Features Beyond the Brief

### 3.1 Gate 1 / Gate 2 Dual Filtering (Step 4)

The brief describes alignment and divergence detection as a single step. In practice, running MAFFT + HyPhy on every orthogroup across 20 species would take days and produce mostly noise. Two gates were added:

- **Gate 1** (before MAFFT): keeps only orthogroups where ≥2 resilient species diverge >15% from human, using a fast global alignment. This reduces the MAFFT workload by ~80%.
- **Gate 2** (after divergence detection): keeps only orthogroups with motifs in ≥N independent lineages (`convergence_min_lineages` from config). This gates what enters HyPhy.

**Config knobs**: `divergence_min_species`, `divergence_identity_max`, `convergence_min_lineages` in `config/environment.yml`.

### 3.2 trimAl Safety Guard

trimAl masking uses `--gt 0.1 --cons 60 -fasta`. A guard was added: if trimAl removes >80% of alignment columns (over-masking an already-short protein), the original unmasked alignment is used instead.

### 3.3 Composite Score Denominator Normalisation

When Phase 2 annotations are absent (zero scores for disease, druggability, safety, regulatory), the composite score denominator now excludes those weights. Without this, a gene with a perfect Phase 1 score would appear at ~0.40 composite because 60% of the denominator was occupied by Phase 2 components that hadn't been computed yet.

The three Phase 1 signals (convergence, selection, expression) always anchor the denominator regardless of their value.

### 3.4 Validated Tier — Precise Logic

The brief describes "Validated" as "Tier 1 score + human genetic evidence". The implementation uses:

```
human_genetics_score ≥ 0.3  →  Tier upgraded to Validated
```

Where `human_genetics_score` is computed separately from the disease score and accumulates:
- GWAS pvalue < 5×10⁻⁸: up to +0.50
- gnomAD pLI: up to +0.30
- OpenTargets genetic association: up to +0.30
- **PCSK9 paradigm**: rare protective variant (gnomAD, pvalue < 5×10⁻⁸) matching animal divergence direction → immediate +0.50, fast-tracks to Validated regardless of composite tier

This logic is separate from the `disease_score()` component used in the composite formula.

### 3.5 Graduated DepMap and GTEx Safety Penalties

The brief mentions DepMap and GTEx as binary flags. The implementation uses graduated thresholds:

**DepMap** (CRISPR chronos score):
- < −0.7: −0.30 (very broadly essential — serious toxicity concern)
- < −0.5: −0.15 (broadly essential)
- < −0.3: −0.05 (moderate essentiality)

**GTEx** (tissue count with TPM > 1):
- > 40 tissues: −0.20 (ubiquitous expression)
- > 25 tissues: −0.10
- > 15 tissues: −0.05

### 3.6 Control Species Divergence Penalty — Formula

Applied in step15 (Phase 2 rescore only). For each gene:

```
if control_divergence_fraction > 0.5:
    score = score × (1 − 0.6 × control_divergence_fraction)
```

At 100% control divergence (gene changes equally in all control species), the score is reduced by 60%. This suppresses non-specific evolutionary divergence signals without eliminating the gene entirely.

### 3.7 Convergence Min-Lineages Auto-Adjustment

If the number of registered non-human resilient lineages in the active species panel is smaller than `convergence_min_lineages` (e.g., running a 3-species test), the threshold is automatically reduced to `max(1, n_available_lineages − 1)`. This prevents every gene from being filtered out during Gate 2 in small validation runs.

### 3.8 RELAX Scoring Formula

The brief lists RELAX as contributing "0.15 × RELAX to the selection score" without defining the component. The formula is:

```
relax_score = min((k − 1) / 2.0, 1.0) × min(−log10(p) / 10, 1.0)
```

Applied only when k > 1 (selection intensification, not relaxation). The k−1 normalisation gives full weight at k = 3 (doubling of selection intensity).

### 3.9 P2Rank Contributes Up to 10% Bonus to Druggability

P2Rank pocket probability contributes `min(score, 1.0) × 0.1` on top of the existing druggability score from fpocket + ChEMBL + CanSAR.

---

## 4. API Enhancements

### 4.1 New Candidate Filters

| Parameter | Description |
|-----------|-------------|
| `?variant_direction=loss_of_function` | Filter by motif direction (GoF / LoF / likely_pathogenic / neutral) |
| `?trait_id=cancer_resistance` | Filter by trait preset (was already supported in scoring but not exposed) |
| `?in_functional_domain=true` | Only genes with ≥1 motif inside a Pfam/InterPro domain |
| `?min_score=0.5` | Minimum composite score threshold |
| `?species_id=bowhead_whale` | Genes with an ortholog in this species |

### 4.2 Expanded Gene Detail Response

The `/candidates/{gene_id}` endpoint now returns three additional data blocks previously missing from the response:

- **`gene_therapy`**: gene size, AAV compatibility, CRISPR sites, off-target risk, tissue tropism
- **`safety`**: hub risk, network degree, DepMap score, GTEx tissue count, PheWAS hits
- **`regulatory`**: per-species AlphaGenome promoter divergence table

---

## 5. Frontend Updates

### 5.1 Pipeline Tracking Page — All 30 Steps Visible

The Pipeline page previously showed 18 entries, hiding all 12 sub-steps (4b, 4c, 4d, 6b, 6c, 7b, 8b, 11b, 11c, 11d, 12b, 14b). These are now visible in the step tracker with correct Phase 1 / Phase 2 grouping and live status from the SSE feed.

### 5.2 Gene Detail Page — Four Missing Sections Added

Per Section 6.2 of the brief, the gene detail view now shows:

1. **Safety Profile** — essentiality, hub risk, DepMap score, GTEx tissue count, PheWAS associations
2. **Gene Therapy Feasibility** — AAV compatibility, gene size, CRISPR sites, tissue tropism
3. **Regulatory Divergence table** — per-species AlphaGenome scores with promoter divergence and expression log2FC
4. **RELAX k and p-value** — added to the Evolutionary Evidence section alongside FEL and BUSTED

### 5.3 Candidates Page — Three Missing Filters Added

Per Section 6.1 of the brief, the candidates table now supports:

1. **Trait preset selector** — runs Cancer Resistance / Longevity / Viral Immunity / Hypoxia / DNA Repair presets
2. **Variant direction filter** — filter by ↓ LoF / ↑ GoF / ⚠ Likely pathogenic / Neutral
3. **"Validated" tier button** — the tier filter now includes Validated (emerald green badge) alongside Tier1/2/3

### 5.4 Pathway Convergence Page (New)

A dedicated `/pathways` page was added, accessible from the sidebar, showing:
- Bar chart of top 20 enriched pathways by combined enrichment + evolutionary weight score
- Full sortable table with pathway ID (Reactome links), candidate count, background gene count, log10 hypergeometric p-value, evolutionary weight, and member gene symbols
- Minimum candidate count slider filter

This implements Section 6.3 of the brief.

### 5.5 Dashboard — Validated Tier Card Added

The dashboard summary cards now show **Validated** (emerald), Tier1, Tier2, Tier3 counts separately, and the bar chart colour-codes Validated genes distinctly.

---

## 6. Validation Test Configuration

A lightweight test configuration (`config/species_registry_test.json`) was created for pipeline validation without a full 19-species run. It was intentionally designed to produce *meaningful* convergence results rather than just confirming the pipeline doesn't crash:

| Species | Lineage | Key Phenotype |
|---------|---------|---------------|
| Naked mole rat | Rodents | Cancer resistance, longevity |
| Bowhead whale | Cetaceans | Longevity, cancer resistance |
| African elephant | Proboscideans | Cancer resistance (TP53 copies) |
| Little brown bat | Bats | Viral tolerance, longevity |
| Greenland shark | Sharks | Longevity (400+ yr) |
| Human | Primates | Reference |
| Rat | Rodents | Control |
| Macaque | Primates | Control |

**Rationale**: 5 independent lineages spanning 450 million years of divergence (mammals to sharks) is the minimum needed for the convergence scoring to produce genuine Tier1 signals. Known benchmark genes (TP53, ATM, ERCC1) are well-documented in exactly this combination of species. The run takes 4–8 hours on a 10-core server.

Run with: `bash scripts/run_validation_test.sh`

---

## 7. What the Brief Says That Is Intentionally Not Yet Implemented

These items are acknowledged limitations, not bugs:

| Item | Status |
|------|--------|
| Gene-tree reconciliation (GeneRax / Notung) | Planned for Phase 2 |
| ML-based scoring (gradient boosting on DrugBank) | Planned for Phase 3 |
| CYP interactions and hERG liability | Not in current safety screen |
| Elephant LIF6 / TP53 retrogene copies as separate entries | Currently excluded (1:1 filter) |

---

*This document should be updated after each production run or significant code change.*  
*Last updated: March 2026*
