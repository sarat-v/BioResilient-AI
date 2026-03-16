# BioResilient AI — Improvements Backlog

This file tracks scientific and engineering improvements identified during the current
pipeline run. Items are classified by priority and whether they require a full or
partial pipeline rerun.

---

## Status Legend

| Symbol | Meaning |
|--------|---------|
| ✅ | Already implemented in current codebase |
| 🔄 | Partially implemented — needs extension |
| ⏳ | Deferred — implement in v2 before next full run |
| 💡 | Future research direction — evaluate when expanding species panel |

---

## Step 4 — Divergence Detection

### ✅ Hard-coded thresholds moved to config

All sliding-window divergence parameters are now in `config/environment.yml` under
`thresholds:` and read dynamically at runtime via `get_thresholds()`:

```yaml
divergence_min_score: 0.15
divergence_top_n_per_ortholog: 8
divergence_window_size: 15
divergence_window_step: 5
divergence_min_species_per_window: 2
```

No code change needed. Tune these values per run in `environment.yml`.

---

### ✅ Gate 2 — lineage recurrence filter

Implemented in `pipeline/layer1_sequence/divergence.py → filter_by_independent_lineages()`.

Requires a motif to appear in ≥ `convergence_min_lineages` (default 3) independent
lineage groups before it is passed downstream to HyPhy. This already suppresses
orthogroups where divergence is driven by a single taxon.

---

### ✅ Top-N motif cap configurable

`divergence_top_n_per_ortholog` is configurable (currently 8). Previously hard-coded
at 3. No further action needed.

---

### ✅ QC metrics after Step 4

`step_reporter.py → _collect_step4()` logs:

- Motif counts per step
- Orthogroups passing Gate 2
- Top 20 genes ranked by distinct species showing divergence

---

### 🔄 Gate 2 — exclude control species

**Status**: Control species (`is_control: true` in `species_registry.json`) are
excluded from the foreground phenotype comparisons in HyPhy selection tests. However,
Gate 2 currently counts all lineages including controls when checking lineage
recurrence.

**What to fix for v2**: In `filter_by_independent_lineages()`, only count species
where `is_control == False` towards the recurrence threshold. Controls (rat, macaque)
should contribute sequence data for the alignment but not be allowed to satisfy the
convergence gate.

**File**: `pipeline/layer1_sequence/divergence.py`

**Effort**: Low — pass `species_to_is_control: dict[str, bool]` into the function
and skip control species when building `lineages_with_motif`.

---

## Step 4 (Scientific) — Phylogenetic Distance Weighting

**Source**: Suggested by the bioresearcher in `Steps 1-4 (Suggestion).pdf`.

### Background

Very distant species (hydra, ocean quahog clam, sharks) naturally diverge from humans
across most proteins simply due to hundreds of millions of years of evolution, not
because of any phenotype-specific adaptation. This can inflate the raw pool of
divergent motifs in Step 4 before downstream gates have a chance to filter them.

### What is already in place

The convergence scoring layer (`pipeline/layer2_evolution/convergence.py`) already
applies a `phylogenetic_convergence_weight()` function that weights convergence by
**mean pairwise evolutionary distance in millions of years** between lineages. This is
more biologically grounded than the manual weights proposed in the PDF, and it
partially compensates for the inflation at the convergence scoring stage.

Additionally, Gate 2 already requires motifs to appear in ≥ 3 independent lineages,
which suppresses signals driven by a single distant species.

### Why the proposed solution is premature

The PDF proposes:

```
weighted_divergence = raw_divergence × species_weight
```

with manually assigned weights (mammals=1.0, fish=0.6, hydra=0.4). As assessed by
independent review, these weights are arbitrary and not biologically grounded. A
principled implementation should use phylogenetic branch lengths from the species tree
(available after Step 5 IQ-TREE), not hand-tuned constants.

### Recommendation for v2

When the species panel expands (especially adding more cnidarians or molluscs), or
when divergent motif counts grow large enough to slow HyPhy, implement
**branch-length-normalised divergence** using the IQ-TREE species tree output:

```python
# Pseudocode
branch_length = species_tree.get_distance("human", species_id)
normalized_score = raw_divergence_score / branch_length
```

This replaces the manual weight table with a biologically derived value and requires
no arbitrary calibration. The species tree is already computed in Step 5 and available
at `{local_storage_root}/phylo/species.treefile`.

**Files to modify**: `pipeline/layer1_sequence/divergence.py`,
`pipeline/orchestrator.py` (pass tree to step4 call), `config/environment.yml`
(add `divergence_apply_branch_norm: false` flag).

**Prerequisite**: Step 5 (IQ-TREE) must have completed at least once to produce the
species tree. On first run this is a chicken-and-egg — use lineage group weights as
a one-time bootstrap, then switch to branch-length normalisation on reruns.

**Priority**: Low — current protections (Gate 2, selection tests, convergence
weighting) are sufficient for the current 18-species panel.

---

## Step 3 (Technical) — UniProt accession in FASTA headers

**Status**: Fixed in `pipeline/layer1_sequence/download.py → reheader_fasta()` for
future proteome downloads.

**Issue**: The current database (`gene.gene_symbol`, `gene.human_protein`) stores
UniProt **mnemonics** (e.g. `BRCA1_HUMAN`) because `reheader_fasta` was extracting
index `[-1]` instead of index `[1]` from `sp|ACCESSION|MNEMONIC` headers.

**Runtime workaround**: `alphamissense.py` and `pfam.py` now call
`resolve_mnemonics_to_accessions()` at runtime to translate mnemonics to accessions
on the fly, with a disk cache. This is transparent to the rest of the pipeline.

**For v2**: After the current run completes, wipe the `gene`, `ortholog`, and
`divergent_motif` tables and rerun from Step 3 with the corrected `reheader_fasta`.
The accessions will be stored correctly from the start, removing the need for the
runtime mnemonic resolution layer (and reducing Step 4b startup time by ~30s).

---

## Infrastructure

### ⏳ RDS storage — upgrade to gp3

Current RDS instance uses `gp2` storage. Upgrading to `gp3` gives 3x baseline IOPS
at the same cost for the same storage size. Beneficial for Step 6 (HyPhy writes).

**Action**: AWS Console → RDS → Modify instance → Storage type: gp3.

### ⏳ Step 4c (ESM-1v) — GPU instance

Step 4c is deferred to run on a GPU instance (`g4dn.xlarge`). The scoring code
(`pipeline/layer1_sequence/esm1v.py`) is already optimised to perform one GPU forward
pass per unique protein sequence and score all motifs from the resulting log-probability
tensor on CPU.

When running on GPU:
1. Launch `g4dn.xlarge` from the current AMI.
2. Confirm `torch.cuda.is_available()` returns `True`.
3. Run `./run_cancer_resistance_stepwise.sh --from step4c`.

---

*Last updated: 2026-03-11*
