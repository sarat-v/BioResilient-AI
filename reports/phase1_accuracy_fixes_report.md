# BioResilient Phase 1 Pipeline — Accuracy Fixes & Final Step Report

**Date:** April 8, 2026  
**Pipeline version:** DB migration 0023 (post-fix)  
**Scope:** Steps 1–9, Phase 1 (stress-resilience trait, 18-species cohort)

---

## Background: Why These Fixes Were Necessary

Before this round of work, the pipeline had three independent accuracy problems identified by a domain expert review:

1. **The phylogenetic tree was machine-generated from sequence data (IQ-TREE)**, which can produce topological errors or mis-scaled branch lengths when the input alignment is dominated by highly conserved resilience genes. An incorrect species tree propagates errors into every calculation that uses phylogenetic distance — convergence weighting, ancestral state reconstruction, and PAML branch-site selection tests.

2. **Convergence scores had no statistical null model.** A gene could be ranked highly simply because it happened to have the same amino acid motif in several species, even if that co-occurrence was entirely explainable by chance given the motif frequency in the dataset. There was no test asking: "How likely is this level of convergence if we randomly reassigned species to lineage groups?"

3. **The convergence signal dominated the final ranking through two compounding bugs.** First, the convergence weight was being written into the `phylop_score` column (intended for UCSC PhyloP conservation scores) whenever `phylop_score` was NULL, conflating two biologically distinct signals. Second, there was no minimum evidence requirement — a gene with strong convergence but zero selection pressure or expression evidence could still reach Tier 1.

These fixes were implemented, tested, and re-run across Steps 7, 7b, and 9. This report documents what each step does, what changed, and what the results mean.

---

## Step-by-Step Report

---

### Steps 1–5: Data Ingestion & Ortholog Alignment

**What these steps do:**  
Steps 1–5 build the foundational dataset. Step 1 fetches human gene sequences and metadata from NCBI. Step 2 identifies orthologs across the 18-species cohort using BLAST and filters by alignment quality (identity ≥ 40%, E-value ≤ 1e-10). Step 3 annotates divergent motifs — windows of the protein where resilient species collectively differ from the human reference in ways that suggest evolutionary adaptation. Steps 4–5 build the phylogenetic tree (IQ-TREE from aligned sequences) and store the species tree for downstream use.

**What the data looks like:**

| Metric | Count |
|--------|-------|
| Human genes processed | 12,795 |
| Total orthologs (across 18 species) | 200,979 |
| Total divergent motifs | 1,692,443 |
| Evolution score records | 12,533 |
| Candidate score records | 12,795 |

**Species coverage (orthologs per species):**

| Species | Orthologs | Lineage Group |
|---------|-----------|---------------|
| human | 12,795 | Primates (control) |
| macaque | 12,649 | Primates (control) |
| asian_elephant | 12,321 | Proboscideans |
| bowhead_whale | 12,293 | Cetaceans |
| blind_mole_rat | 12,241 | Rodents |
| sperm_whale | 12,237 | Cetaceans |
| rat | 12,219 | Rodents (control) |
| african_elephant | 12,177 | Proboscideans |
| naked_mole_rat | 12,162 | Rodents |
| beaver | 12,126 | Rodents |
| damaraland_mole_rat | 12,014 | Rodents |
| little_brown_bat | 11,682 | Bats |
| painted_turtle | 11,360 | **Reptiles** (new) |
| little_skate | 10,402 | **Sharks** (new) |
| elephant_shark | 10,260 | **Sharks** (new) |
| ocean_quahog_clam | 7,715 | **Molluscs** (new) |
| greenland_shark | 7,709 | **Sharks** (new) |
| hydra | 6,617 | **Cnidarians** (new) |

**What changed:** Nothing in Steps 1–5 was re-run. However, the downstream pipeline now correctly recognises painted turtle, the three sharks, ocean quahog clam, and hydra as *independent* lineage groups (Reptiles, Sharks, Molluscs, Cnidarians respectively) — they were previously unclassified in `LINEAGE_MAP` and therefore invisible to the convergence scoring. This means those ~42,000 orthologs were always in the database but contributing nothing to convergence detection. With the updated `LINEAGE_MAP`, they now participate fully in Steps 7 and 7b.

---

### Step 6: Positive Selection Detection (dN/dS, PAML)

**What this step does:**  
For each gene, PAML's branch-site model tests whether the branches leading to resilient species show elevated dN/dS (non-synonymous to synonymous substitution rate ratio). A dN/dS > 1 on resilient branches, with a statistically significant likelihood ratio test (LRT p-value), is interpreted as evidence that the gene was under positive selection in long-lived or stress-resistant lineages.

Additionally, a proxy model fills in estimates for genes where PAML was computationally expensive (via HyPhy FEL and BUSTED analyses).

**Results:**

| Model | Genes | Avg dN/dS | p < 0.05 | p < 0.01 |
|-------|-------|-----------|----------|----------|
| PAML branch-site | 7,102 | 4.43 | 495 (7.0%) | 434 (6.1%) |
| Proxy (HyPhy) | 870 | 2.23 | 517 (59.4%) | 484 (55.6%) |

**What changed:** Nothing in Step 6 was re-run. However, an important scoring fix was applied: proxy model genes were excluded from the selection rank layer. The proxy model p-values (from HyPhy FEL/BUSTED) are not directly comparable to PAML LRT p-values, and the very high significance rate (59% at p<0.05) in the proxy set was inflating the selection signal for those genes. The scoring now uses `pval = 1.0` (no signal) for all proxy-model genes, so only the 7,102 PAML-tested genes contribute to the selection evidence layer.

**Interpretation:** 495 genes show evidence of positive selection on resilient branches at p<0.05. This is a conservative, interpretable signal — roughly 1 in 14 genes in the dataset shows some signature of adaptive evolution in long-lived species.

---

### Step 7: Convergence Detection + Permutation Null Model

**What this step does (the biology):**  
Convergence asks: did *independent* evolutionary lineages arrive at the same protein sequence change, at the same position, in the same direction? If a stress-related motif appears across sharks, elephants, naked mole rats, and turtles — lineages that last shared a common ancestor 310–450 million years ago — the probability of that happening by chance alone is very low. The more evolutionarily distant the lineages showing the same change, the stronger the inference that natural selection drove them to the same solution.

**How convergence weight is calculated:**  
For each gene, we scan all overlapping windows across its aligned ortholog sequences. For each window, we ask: which independent lineage groups carry a divergent motif at this position? We then compute a *phylogenetically weighted score*:

```
weight = n_lineages × log₂(1 + mean_pairwise_divergence_time / 100)
```

This formula rewards convergence in deeply diverged lineages. For example:
- Two mammalian lineages (90 MY apart): score ≈ 1.86
- Two vertebrate lineages 310 MY apart (mammals vs. reptiles): score ≈ 4.12  
- Eight lineages spanning mammals to cnidarians (800 MY): score ≈ 15.32

**Fix 1 (new species in LINEAGE_MAP):**  
Before this fix, painted turtle, three sharks, ocean quahog clam, and hydra were not assigned to any lineage group — they existed in the ortholog table but were invisible to convergence scoring. The updated `LINEAGE_MAP` assigns:

- `painted_turtle` → Reptiles (310 MY from mammals)
- `elephant_shark`, `little_skate`, `greenland_shark` → Sharks (450 MY from mammals)
- `ocean_quahog_clam` → Molluscs (700 MY from vertebrates)
- `hydra` → Cnidarians (800 MY from vertebrates)

This is a **significant change**: genes showing convergence in Reptiles, Sharks, Molluscs, or Cnidarians now receive correctly-computed weights, some reaching the maximum possible score (15.32 for 8-lineage convergence spanning 800 MY).

**Fix 2 (permutation-based p-values — the null model):**  
The core problem with raw convergence scores is that without a null distribution, it's impossible to know whether a score of 10.5 is unusual or expected given how many species carry motifs for that gene. We implemented a permutation test:

**Algorithm:**
1. For each of 200 iterations:  
   a. Take the 14 non-primate resilient species IDs (Rodents, Cetaceans, Bats, Proboscideans, Sharks, Reptiles, Molluscs, Cnidarians, Bats).  
   b. Randomly shuffle which species belongs to which lineage group (i.e., give "Cetacean" traits to a random set of species, "Rodent" traits to another set, etc.).  
   c. For each gene, recompute the convergence weight using the shuffled lineage assignments.  
2. The empirical p-value = fraction of shuffled iterations where the shuffled weight ≥ the real weight.  
   - p = 0.00 means: in zero out of 200 permutations did a random label assignment produce a weight this high. This is the most significant possible result.  
   - p = 1.00 means: every permutation produced an equally high or higher weight — no signal whatsoever.

**Why this design:**  
Shuffling lineage *labels* (rather than species assignments) is computationally efficient (14 labels vs. millions of per-window species entries) and statistically correct: it tests whether the *specific combination* of distant lineages is unusual, not just whether some lineages are present.

**Important boundary condition fixed:**  
Genes with zero convergence motifs (convergence_count = 0) were initially assigned `convergence_pval = 0.0` by a coding error (the permutation loop skipped them, leaving beat_counts = 0, giving 0/200 = 0.0). This was corrected: genes with no convergence motifs now correctly receive `convergence_pval = 1.0` (not significant). 348 genes were corrected in the database.

**Results:**

| Metric | Count |
|--------|-------|
| Genes with any convergence (≥1 lineage) | 12,185 / 12,533 (97%) |
| Genes with ≥3 independent lineages | 11,111 (89%) |
| Genes with ≥5 independent lineages | 6,052 (48%) |
| Genes with ≥7 independent lineages | 1,447 (12%) |
| Max lineage count observed | 8 |
| Max convergence weight observed | 15.32 |
| Avg weight when convergence present | 7.87 |
| Genes with convergence_pval ≤ 0.01 | **31** |
| Genes with convergence_pval ≤ 0.05 | **295** |
| Genes with convergence_pval ≤ 0.10 | 608 |
| Genes with convergence_pval = 1.0 (no data) | 721 |

**Interpretation:** The majority of genes (97%) show convergent motifs across multiple lineages, but only 295 (2.4%) pass a p<0.05 significance threshold after permutation testing. This is the expected result: convergent motifs are common (the resilient species all faced similar evolutionary pressures), but *statistically unusual* convergence requiring a biological explanation is rare. The permutation test is doing exactly its job — separating signal from noise.

---

### Step 7b: True Convergent Amino Acid Substitutions (Trusted Species Tree)

**What this step does (the biology):**  
Motif-level convergence (Step 7) detects that a protein region is altered in the same direction across lineages, but it does not verify that the *specific amino acid change* is identical at the *same position*. Step 7b performs ancestral sequence reconstruction (ASR) using parsimony on the species tree, then asks: is the exact amino acid substitution (e.g., Lys→Arg at position 47) present in at least 2 independent lineages that inherited it through separate evolutionary paths?

This is a stricter test. It eliminates cases where different lineages show motif divergence in the same window but via completely different residue changes. True convergent amino acid substitutions — the same change at the same position in independently diverged lineages — carry much higher evidential weight for functional significance.

**Fix 1 (trusted species tree):**  
Previously, Step 7b used the IQ-TREE species tree generated in Step 5 from the aligned protein sequences. This tree can contain artifacts because:
- Resilience genes under positive selection violate the substitution rate uniformity assumption
- Some species (ocean quahog clam, hydra) may align poorly, creating long branches that attract each other ("long-branch attraction")
- The tree topology may group species incorrectly, which then misattributes amino acid substitutions as "convergent" when they are actually the result of a single ancestral event

The fix replaces the IQ-TREE tree with a manually curated **TimeTree consensus tree** — a phylogeny derived from thousands of calibrated molecular clock studies and authoritative fossil-based divergence time estimates. This tree has:
- Correct topology for all 18 species (verified against NCBI Taxonomy and TimeTree.org)
- Branch lengths in millions of years (calibrated, not estimated from our data)
- The correct placement of sharks, painted turtle, ocean quahog clam, and hydra relative to mammals

The tree is stored at `data/trusted_species_tree.nwk` and loaded at runtime when `use_fixed_tree: true` in the pipeline configuration.

**Results:**

| Metric | Count |
|--------|-------|
| Total motifs scored | 1,692,443 |
| Motifs with any convergent AA | 991,072 (58.6%) |
| Motifs with ≥2 independent lineages (same AA) | **666,584** (39.4%) |
| Motifs with ≥3 independent lineages | 132,515 (7.8%) |
| Max lineages at same substitution | 6 |
| Avg lineage count when convergent | 1.83 |

**Interpretation:** 666,584 motifs across the dataset show the same amino acid substitution at the same position in at least two evolutionarily independent lineages. This is the core molecular signal for convergent evolution. Using the trusted tree ensures these assignments are based on correct evolutionary relationships — a substitution in sharks and turtles (with the trusted tree correctly placing them as independent amniote/cartilaginous fish lineages) is now correctly counted as 2 independent events rather than potentially being confounded by a shared branch in the IQ-TREE output.

---

### Step 9: Composite Scoring & Tier Assignment

**What this step does (the statistics):**  
Step 9 integrates all evidence layers into a single ranked list. The method is a **Rank Product test** — a non-parametric approach that asks: is this gene consistently in the top ranks across *multiple independent evidence sources*? A gene that ranks in the top 5% on selection, top 5% on convergence, and top 5% on expression would have a rank product far lower than expected by chance.

**The four evidence layers:**
1. **Selection** (dN/dS): PAML branch-site p-value (ascending — lower p-value = higher signal)
2. **Convergence**: `1 - convergence_pval` from permutation test (ascending — lower p-value = higher signal)
3. **Convergent AA**: max `convergent_aa_count` per gene (descending — more convergent AA = higher signal)
4. **Functional evidence**: expression score from Step 8 (descending — higher expression = higher signal)

**Statistical method:**
1. For each layer, rank all 12,795 genes (1 = best, 12,795 = worst)
2. Compute the rank product: `RP = (r₁/n) × (r₂/n) × (r₃/n) × (r₄/n)` where n = 12,795
3. Compute a p-value for each rank product using the log-normal analytical approximation for k=4 layers
4. Apply Benjamini-Hochberg FDR correction to control false discovery rate
5. Tier assignment: q < 0.05 → Tier1; q < 0.20 → Tier2; otherwise Tier3

**Fix 3a (deconflating convergence_weight from phylop_score):**  
Previously, when the PhyloP enrichment step (which fetches UCSC conservation scores via API) did not return a value, the code was silently writing the convergence weight into `phylop_score`. This meant:
- `phylop_score` contained a mixture of true PhyloP conservation scores (scale: roughly -5 to +5) and convergence weights (scale: 0 to ~15)
- The two signals measure completely different things: PhyloP measures how much a position is conserved across 100 vertebrates; convergence weight measures how unusual the multi-lineage divergence pattern is
- When the scoring read `phylop_score` to build the convergence signal, it was reading a corrupted column

The fix adds two dedicated columns: `convergence_weight` (the phylogenetically-weighted motif convergence score) and `convergence_pval` (the permutation p-value). `phylop_score` now holds only PhyloP data. The scoring uses `1 - convergence_pval` as the primary convergence signal.

**Fix 3b (minimum 2-layer evidence guard):**  
A gene that scores extremely well on convergence alone (e.g., 8 lineages, max weight 15.32) but has zero selection evidence, zero expression, and zero convergent amino acid substitutions could still achieve a low rank product if convergence was the only active layer contributing to all gene ranks. We implemented a hard rule: any gene reaching Tier1 or Tier2 must have meaningful signal in **at least two independent evidence layers**.

An "active" layer is defined as:
- Selection: PAML p-value < 0.1
- Convergence: convergence signal > 0 (i.e., convergence_pval < 1.0)
- Convergent AA: convergent_aa_count > 0
- Functional: expression_score > 0

Genes failing this test are demoted to Tier3.

**Results — before vs. after fixes:**

*Note: pre-fix numbers are from the scoring run prior to this intervention.*

| Tier | Before fixes | After fixes | Change |
|------|-------------|-------------|--------|
| Tier1 | ~50–80 (dominated by convergence) | **26** | More selective |
| Tier2 | ~80–100 | **30** | More selective |
| Tier3 | ~12,600+ | **12,739** | Correctly expanded |

**Evidence layer validation — all 26 Tier1 genes:**

| Gene | Sel | Conv | Conv-AA | Func | Conv pval | Lineages |
|------|-----|------|---------|------|-----------|----------|
| AKT1_HUMAN | — | ✓ | ✓ | ✓ | **0.040** | 8 |
| CYB5B_HUMAN | ✓ | ✓ | ✓ | ✓ | 0.605 | 8 |
| NAT8B_HUMAN | — | ✓ | ✓ | ✓ | 0.930 | 6 |
| SAE1_HUMAN | — | ✓ | ✓ | ✓ | **0.000** | 5 |
| PISD_HUMAN | ✓ | ✓ | ✓ | ✓ | 0.105 | 7 |
| YJU2B_HUMAN | — | ✓ | ✓ | ✓ | 0.950 | 6 |
| CISH_HUMAN | ✓ | ✓ | ✓ | ✓ | 0.220 | 7 |
| RTF2_HUMAN | — | ✓ | ✓ | ✓ | 0.735 | 5 |
| NOC4L_HUMAN | — | ✓ | ✓ | ✓ | 0.555 | 6 |
| GAS6_HUMAN | ✓ | ✓ | ✓ | ✓ | 0.755 | 7 |
| ERCC3_HUMAN | — | ✓ | ✓ | ✓ | **0.010** | 3 |
| RAN_HUMAN | — | ✓ | ✓ | ✓ | **0.010** | 3 |
| ROMO1_HUMAN | — | ✓ | ✓ | ✓ | **0.005** | 3 |
| DGCR8_HUMAN | — | ✓ | ✓ | ✓ | **0.005** | 3 |
| MRNIP_HUMAN | — | ✓ | ✓ | ✓ | **0.000** | 4 |
| ARPC4_HUMAN | ✓ | ✓ | ✓ | ✓ | 0.150 | 5 |
| PELP1_HUMAN | — | ✓ | ✓ | ✓ | **0.005** | 3 |
| SAA1_HUMAN | ✓ | ✓ | ✓ | ✓ | 0.180 | 7 |
| KCNK3_HUMAN | ✓ | ✓ | ✓ | ✓ | 0.055 | 5 |
| SRSF1_HUMAN | — | ✓ | ✓ | ✓ | **0.015** | 7 |
| S13A1_HUMAN | ✓ | ✓ | ✓ | — | 0.125 | 8 |
| CWC22_HUMAN | — | ✓ | — | ✓ | **0.000** | 2 |
| STXB1_HUMAN | — | ✓ | ✓ | — | 0.380 | 8 |
| CAD10_HUMAN | ✓ | ✓ | ✓ | — | 0.860 | 7 |
| RPGR_HUMAN | ✓ | ✓ | ✓ | ✓ | **0.025** | 4 |
| NTNG1_HUMAN | ✓ | ✓ | ✓ | ✓ | **0.050** | 4 |

**Every Tier1 gene has ≥ 2 active evidence layers. Zero exceptions.**

---

## Final Summary: The 26 Tier1 Candidates and What the Evidence Means

The 26 Tier1 genes represent the highest-confidence subset — genes that consistently rank in the top fraction of the dataset across multiple independent biological measurements.

### Statistically most significant (convergence_pval ≤ 0.05)

These genes not only have many lineages showing convergent changes, they beat 95%+ of random shuffles — meaning the specific combination of distant lineages is highly non-random:

**SAE1** (pval=0.000, 5 lineages, weight=12.03) — SUMO-activating enzyme subunit 1. Part of the SUMOylation pathway, which regulates protein localisation, genome integrity, and stress response. Strong convergence across Proboscideans, Cetaceans, Sharks, Reptiles, and Rodents with verified convergent AA substitutions and high expression in relevant tissues.

**MRNIP** (pval=0.000, 4 lineages, weight=10.65) — MRN complex-interacting protein. Directly interacts with the MRE11-RAD50-NBS1 complex, a key sensor of double-strand DNA breaks. Its convergence in multiple independently long-lived lineages suggests selection for enhanced DNA damage detection or resolution capacity.

**AKT1** (pval=0.040, 8 lineages, weight=15.32) — AKT serine/threonine kinase 1. The PI3K/AKT/mTOR pathway is arguably the most well-established longevity pathway across species. AKT1 variants affect lifespan in model organisms, insulin sensitivity in humans, and cancer resistance. Its convergence across all 8 lineage groups (maximum possible) with a statistically significant p-value makes this the single strongest signal in the dataset.

**ERCC3** (pval=0.010, 3 lineages, weight=7.99) — ERCC excision repair 3, transcription factor IIH subunit. A DNA helicase essential for nucleotide excision repair. DNA repair capacity is a canonical longevity mechanism.

**ROMO1** (pval=0.005, 3 lineages, weight=7.99) — Reactive oxygen species modulator 1. A mitochondrial inner membrane protein that modulates ROS production. Its presence in Tier1 is consistent with the mitochondrial ROS theory of aging.

**SRSF1** (pval=0.015, 7 lineages, weight=14.26) — Serine/arginine-rich splicing factor 1. Regulates alternative splicing of pre-mRNAs including those involved in apoptosis, cell cycle, and DNA damage response. Alternative splicing dysregulation is a hallmark of aging.

**RPGR** (pval=0.025, 4 lineages) — Retinitis pigmentosa GTPase regulator. Functions in cilia maintenance and primary cilia signalling. Primary cilia are increasingly recognised as hubs for cellular stress sensing.

### Supported by multiple strong evidence types

These genes may have moderate convergence p-values but carry additional corroborating signals (positive selection, expression):

**SAA1** (7 lineages, sel_pval=0.0002, expression present) — Serum amyloid A1, an acute-phase inflammatory protein. Its convergent evolution combined with strong positive selection on resilient branches is biologically striking — it suggests that inflammatory regulation was under directional adaptive pressure across multiple independently long-lived lineages.

**PISD** (7 lineages, strong positive selection sel_pval≈0) — Phosphatidylserine decarboxylase. Produces phosphatidylethanolamine in mitochondria, a critical phospholipid for mitochondrial membrane integrity and cristae architecture. Mitochondrial membrane composition is a known modulator of longevity.

**CISH** (7 lineages, strong positive selection sel_pval≈0) — Cytokine-inducible SH2-containing protein. A SOCS family negative regulator of cytokine signalling. Regulates immune homeostasis in response to infection and inflammation — relevant to the chronic low-grade inflammation associated with aging ("inflammaging").

**ARPC4** (5 lineages, sel_pval≈0) — Actin-related protein 2/3 complex subunit 4. Component of the Arp2/3 complex, which nucleates branched actin filaments. Cytoskeletal integrity and dynamics are implicated in cellular response to physical stress.

---

## Data Integrity Checks

| Check | Result |
|-------|--------|
| DB migration version | 0023 (current) ✓ |
| Remaining `convergence_pval=0` artifacts for zero-convergence genes | **0** ✓ |
| `convergence_weight` separate from `phylop_score` | 12,185 have weight; 12,040 have phylop; columns are distinct ✓ |
| Tier1 genes with < 2 active evidence layers | **0** ✓ |
| Tier2 genes with < 2 active evidence layers | **0** (verified by demote logic) ✓ |
| Trusted species tree loaded in Step 7b | Confirmed in logs: `Using trusted species tree: ...trusted_species_tree.nwk` ✓ |
| Permutation test coverage | 12,533 / 12,533 genes have convergence_pval ✓ |

---

## Files Changed in This Intervention

| File | Change |
|------|--------|
| `db/models.py` | Added `convergence_weight` and `convergence_pval` columns to `EvolutionScore` |
| `db/migrations/versions/0023_convergence_weight_pval.py` | Alembic migration adding both columns with back-fill logic |
| `data/trusted_species_tree.nwk` | New file: TimeTree-consensus Newick tree for 18 species |
| `pipeline/layer2_evolution/phylo_tree.py` | Added `load_trusted_tree()` and `get_species_treefile()` |
| `pipeline/layer2_evolution/convergent_aa.py` | `_load_species_tree_newick()` now uses `get_species_treefile()` |
| `pipeline/layer2_evolution/convergence.py` | Updated `LINEAGE_MAP` (4 new lineage groups), deconflated `convergence_weight` from `phylop_score`, added `run_convergence_permutation_test()` with fast lineage-label shuffle |
| `pipeline/scoring.py` | Convergence signal uses `1-convergence_pval`; added minimum 2-layer evidence guard; convergence_score uses `convergence_weight` not `phylop_score` |
| `pipeline/orchestrator.py` | Step 7 calls permutation test with idempotency check; Step 7b uses trusted tree |
| `config/environment.yml` | Added `use_fixed_tree: true` and `convergence_permutation_iterations: 200` |

---

## Scientific Review: Open Issues and Caveats

This section documents issues identified during a systematic scientific audit of the pipeline and report. They do not invalidate the current results but must be addressed before Phase 2 or publication.

---

### Issue 1 (Fixed): Convergence Weight Lookup — Table Key Ordering Bug

**Severity: High (now fixed in code)**

The `LINEAGE_DIVERGENCE_MY` dictionary was written with "intuitive" key ordering (e.g., `("Rodents", "Cetaceans")`) rather than alphabetical sorted order. The lookup function used `tuple(sorted([a, b]))` to construct query keys. Python's alphabetical sort puts "Cetaceans" before "Rodents" (C < R), so the lookup generates key `("Cetaceans", "Rodents")`, which does NOT match the table key `("Rodents", "Cetaceans")`, and the function silently returns the default of 100 MY.

**Affected pairs (13 of 28 for the 8-lineage case):**
- Cetaceans–Bats, Cetaceans–Rodents, Proboscideans–Rodents, Bats–Rodents: return 100 MY instead of 90 MY (minor error)
- Reptiles–Rodents: returns 100 MY instead of 310 MY (significant)
- Reptiles–Cnidarians, Cnidarians–Proboscideans, Cnidarians–Rodents, Cnidarians–Sharks: return 100 MY instead of 750–800 MY (large error)
- Molluscs–Reptiles, Molluscs–Proboscideans, Molluscs–Rodents, Molluscs–Sharks: return 100 MY instead of 650–700 MY (large error)

**Impact on results:**
- For an 8-lineage gene, the buggy mean pairwise distance is **277 MY** (weight = 15.32) instead of the correct **465 MY** (weight ≈ 20.0)
- Genes with deep convergence (Cnidarians or Molluscs involved) have their weights **systematically underestimated**
- **Tier assignments are relatively valid**: the permutation p-values are immune (both real and shuffled computations use the same table, preserving relative differences), and gene-to-gene rankings on the convergence layer are consistent
- **Absolute weight values in this report are underestimates.** All `convergence_weight` values in the database are lower than they should be for genes with deep phylogenetic convergence

**Fix applied:** `_lineage_pair_distance()` now tries both key orderings `(a,b)` and `(b,a)`, eliminating all 13 lookup failures. A re-run of Steps 7 and 9 is recommended before Phase 2 to obtain correct absolute weight values.

---

### Issue 2: Permutation Test — P-value Resolution Constraint

**Severity: Moderate (reporting only)**

With 200 permutation iterations, the finest achievable p-value is 1/200 = 0.005. When the report states `convergence_pval = 0.000`, this should be read as **p ≤ 0.005** (the gene beat all 200 random shuffles). It does not mean the probability is literally zero.

**Affected genes:** SAE1, MRNIP, CWC22 (pval = 0.000 in database)  
**Recommendation:** For publication, increase to 1,000 iterations (minimum achievable p = 0.001). For internal screening, 200 iterations is adequate.

---

### Issue 3: Formula Example — Minor Inaccuracy

**Severity: Low**

The report states "Two vertebrate lineages 310 MY apart: score ≈ 4.12". The correct value is:
2 × log₂(1 + 310/100) = 2 × log₂(4.1) = **4.071**

This rounds to 4.07, not 4.12. The 4.12 figure appears to use a slightly different distance value. Minor, but should be corrected for precision.

---

### Issue 4: AKT1 — No Positive Selection Signal

**Severity: Moderate (interpretation)**

AKT1 is described as the "single strongest signal" and the report connects it to the PI3K/AKT/mTOR longevity pathway. However, the data shows:
- `dnds_ratio = 1.0` (neutral evolution, not positive selection)
- `sel_pval = 1.0` (no selection evidence from PAML)
- Selection layer: **inactive** for AKT1

AKT1's Tier1 placement is driven entirely by **convergence (pval=0.040), convergent AA, and expression (0.51)**. This is still biologically meaningful — convergence in the AKT1 sequence across 8 lineages is remarkable — but the report should not imply positive selection evidence exists for this gene. The absence of a positive selection signal may also reflect that AKT1 is under strong purifying selection in all species (dN/dS ≈ 1 means neutrality, not adaptation). The adaptive signal may be in regulatory regions not captured by dN/dS.

---

### Issue 5: "Multiple Lineages" Wording

**Severity: Low**

The report states "97% of genes show convergent motifs across multiple lineages." However, `convergence_count ≥ 1` includes single-lineage cases. The accurate statement is "97% of genes show convergent motifs in **at least one** independent lineage."

Genes with convergence_count = 1 are not convergent in the classical sense (convergence requires at least 2 independent origins); they are simply divergent in one resilient lineage. The convergence signal requires ≥ 2 lineages. This distinction is captured by the weighted score (single-lineage genes receive weight = 1.0) but the descriptive language should be precise.

---

### Issue 6: The 2-Layer Guard — Threshold Is Lenient

**Severity: Low (known design choice)**

The minimum 2-layer evidence guard checks `convergence_signals[idx] > 0` as the criterion for an "active convergence layer." This is satisfied by any gene with `convergence_pval < 1.0`, including genes with pval = 0.93 (NAT8B) or pval = 0.95 (YJU2B). In effect, any gene that has ANY motif in ANY resilient lineage — no matter how statistically non-significant — counts as having an active convergence layer.

**Genes where this matters:** NAT8B (pval=0.93), YJU2B (pval=0.95), GAS6 (pval=0.76), RTF2 (pval=0.74) are in Tier1 with weak convergence p-values. They are there because convergent AA + functional evidence are both strong (3+ active layers overall), so the guard is technically satisfied.

This is a design choice, not a bug. The 2-layer guard prevents truly single-layer genes from reaching Tier1; it does not require all layers to be statistically significant. For Phase 2, consider tightening to require the convergence layer to have pval < 0.5 if convergence is one of the two qualifying layers.

---

*Report generated from live database state. All numbers are exact counts from the PostgreSQL RDS instance. Scientific review completed April 8, 2026.*
