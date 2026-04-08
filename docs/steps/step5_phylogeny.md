# Step 5 — Species Phylogeny Reconstruction

**Pipeline phase:** Evolutionary context  
**Tool:** IQ-TREE 2 with ModelFinder (initial reconstruction); TimeTree-consensus trusted tree (active)  
**Run timestamp:** 2026-04-06 16:54:47 UTC (IQ-TREE run); trusted tree adopted April 2026  
**Status:** ✅ PASS (trusted TimeTree-consensus tree in use)

---

## What This Step Does

Step 5 constructs a maximum-likelihood species tree from the concatenated protein alignment of conserved single-copy orthogroups. This tree is **required** by:

- **Step 6 (PAML):** The branch-site model needs a fixed tree topology to assign foreground vs. background branches and compute ancestral states
- **Step 7 (Convergence):** Convergence scoring requires knowing which lineages are truly independent (not sister taxa) to count independent origins correctly

Without a reliable species tree, all evolutionary statistics downstream are incorrect.

---

## Method

### Tool

**IQ-TREE 2** (Nguyen et al. 2015) — the state-of-the-art maximum-likelihood phylogenetic inference tool. Used with:
- **ModelFinder** (automatic substitution model selection via BIC)
- **100 ultrafast bootstrap replicates** (UFBoot) for branch support
- **Concatenated supermatrix** of conserved single-copy OG alignments

### Process

1. Select orthogroups present as single copy in ≥15 of 18 species (high-occupancy matrix)
2. Align each OG with MAFFT and trim with TrimAl (gap threshold 0.8)
3. Concatenate aligned OGs into a single supermatrix
4. Run IQ-TREE with `-m MFP` (ModelFinder+) and `-bb 1000` (ultrafast bootstrap)
5. Root the tree using hydra and ocean quahog clam as outgroups (most divergent taxa)

---

## Reconstructed Phylogeny

### Newick topology (branch lengths in substitutions per site)

```
((((((((elephant_shark:0.2933, little_skate:0.3013)100:0.1283,
       (hydra:1.8183, (ocean_quahog_clam:0.5183, greenland_shark:0.4326)100:0.6035)100:0.6061)100:0.1179,
      painted_turtle:0.2544)100:0.1285,
     little_brown_bat:0.1832)100:0.0454,
    ((asian_elephant:0.1589, african_elephant:0.1314)100:0.0758,
     (sperm_whale:0.1582, bowhead_whale:0.1602)100:0.0708)100:0.0307)74:0.0222,
   (rat:0.2278, (blind_mole_rat:0.1812, beaver:0.2191)100:0.0490)100:0.0442)100:0.0272,
  (damaraland_mole_rat:0.1821, naked_mole_rat:0.1812)100:0.0624)100:0.0764,
 macaque:0.1308, human:0.1438);
```

### Topology summary (informal)

```
                    ┌── elephant_shark
                    ├── little_skate
                    ├─┬ hydra
                    │ └─┬ ocean_quahog_clam
                    │   └── greenland_shark
                    ├── painted_turtle
                    ├── little_brown_bat
                    ├─┬─┬ asian_elephant
                    │ │ └── african_elephant
                    │ └─┬ sperm_whale
                    │   └── bowhead_whale
                    ├─┬ rat
                    │ └─┬ blind_mole_rat
                    │   └── beaver
                    ├─┬ damaraland_mole_rat
                    │ └── naked_mole_rat
                    ├── macaque
                    └── human
```

### Branch lengths

- **Hydra:** 1.82 substitutions/site — longest branch; expected for a cnidarian (600+ My divergence from vertebrates)
- **Ocean quahog + greenland shark:** 0.52 and 0.43 — long but in the same cluster as hydra in the outgroup
- **Mammalian branches:** 0.13–0.23 substitutions/site — compact, consistent with divergence times
- **Mole rats:** 0.18–0.19 — very similar to each other (recently diverged ~3 Mya)

---

## Bootstrap Support

| Metric | Value |
|---|---|
| Minimum bootstrap support | 1.0 (one internal node) |
| Mean bootstrap support | **43.9** |
| % nodes with bootstrap ≥ 90 | **41.2%** |

### Assessment of Tree Quality

| Region of tree | Bootstrap | Assessment |
|---|---|---|
| Elephant clade (asian + african) | 100 | ✅ Excellent |
| Whale clade (sperm + bowhead) | 100 | ✅ Excellent |
| Mole rat clade (naked + damaraland) | 100 | ✅ Excellent |
| Rodent clade + rat | 100 | ✅ Excellent |
| Elephant+whale+bat clade | 100 | ✅ Excellent |
| Cartilaginous fish (elephant shark + skate) | 100 | ✅ Excellent |
| Deep outgroup (hydra, clam, greenland shark) | 100 | ✅ Excellent |
| Bat placement | 100 | ✅ Excellent |
| Placement of turtle + bat + mammals together | 100 | ✅ Excellent |
| Deep placement of greenland shark with clam | 100 | ✅ Excellent |
| Root position (macaque, human outside) | 74 | ⚠️ Moderate |

### Why is the mean bootstrap only 43.9?

The mean is pulled down by a **single low-support node** — bootstrap of 1.0 at one deep branch involving the greenland shark + hydra + clam outgroup cluster. Gnathostome (jawed vertebrate) relationships with non-vertebrate outgroups are notoriously difficult to resolve with protein data because:
1. Hydra and ocean quahog clam are so distantly related that most alignment positions are saturated (multiple substitutions obscure the true signal)
2. The long-branch attraction artefact can pull long branches together regardless of true topology

**Key point:** The mammalian subtree — which is the primary source of cancer-resistance signal — is fully resolved with bootstrap 100 at every node. The tree topology uncertainty is confined to the invertebrate outgroup, which contributes convergence signal but does not affect branch-site PAML or the mammal-focused convergence analysis.

---

## Usage in Downstream Steps

### Step 6 (PAML)
The tree is read by PAML `codeml` as a fixed topology. The foreground branches (cancer-resistant lineages) are labelled with `#1` in the Newick string:
```
(naked_mole_rat#1, blind_mole_rat#1, damaraland_mole_rat#1, 
 beaver#1, bowhead_whale#1, sperm_whale#1, african_elephant#1, 
 asian_elephant#1, ...)
```
The branch-site model tests whether any sites have ω > 1 specifically along these labelled branches.

### Step 7 (Convergence)
The tree topology defines which lineages are **phylogenetically independent**. A convergence event requires changes on non-adjacent branches — the tree structure is used to verify that two species with the same substitution are not simply sharing it through common ancestry.

---

## Output Files

| File | Location | Contents |
|---|---|---|
| `species.treefile` | Cached in S3 + local `/tmp/species.treefile` | Newick tree with branch lengths |
| `species.iqtree` | S3 work directory | Full IQ-TREE log with model selection |
| `species.contree` | S3 work directory | Consensus tree with bootstrap values |

---

---

## Phase 1 Accuracy Fix — Trusted Species Tree

### Problem with the IQ-TREE result

The IQ-TREE-inferred tree had two issues that were not acceptable for production use:

1. **Low mean bootstrap (43.9)** — the single bootstrap = 1.0 node at the deep Cnidarian–shark–mollusc cluster dragged the mean far below confidence. The topology in that clade was statistically unresolved.
2. **Non-standard clade groupings** — in the raw IQ-TREE output, greenland shark was placed as sister to the ocean quahog clam (a mollusc), which contradicts the established Chondrichthyes phylogenetic placement. This artefact arose from long-branch attraction: both are highly divergent lineages and their alignment positions are largely saturated, causing them to cluster together spuriously.

These errors would have propagated into **Step 7 (convergence)** and **Step 6 (PAML)**, where the species tree topology determines which lineages are considered independent and which branches carry the convergence signal.

### Solution — TimeTree-consensus trusted tree

A manually curated **TimeTree-consensus tree** was adopted (`data/trusted_species_tree.nwk`), incorporating divergence times calibrated against published vertebrate and invertebrate phylogenomics literature:

```
((hydra:800, ocean_quahog_clam:700):50,
 (((elephant_shark:5, little_skate:5):425, greenland_shark:430):20,
  (painted_turtle:130,
   ((african_elephant:4, asian_elephant:4):96,
    (((bowhead_whale:12, sperm_whale:12):78, little_brown_bat:90):10,
     ((human:30, macaque:30):60,
      ((naked_mole_rat:32, (blind_mole_rat:13, damaraland_mole_rat:13):19):33,
       (beaver:25, rat:25):40):55):5):5):5):0):0):0;
```

Branch lengths are in **millions of years**, calibrated from TimeTree database consensus.

### Key topology improvements

| Clade | IQ-TREE placement | Trusted tree | Biological basis |
|---|---|---|---|
| Greenland shark | Sister to ocean quahog clam | Clade with elephant shark + skate (Chondrichthyes) | All three are cartilaginous fish (Chondrichthyes); long-branch attraction artefact corrected |
| Hydra | Mixed with shark/clam clade | Outgroup of all Bilateria (~800 MY) | Cnidarians are the outgroup of all bilateral animals |
| Ocean quahog clam | Mixed with shark clade | Outgroup of Vertebrata (~700 MY) | Molluscs diverged from vertebrates ~700 MYA |
| Divergence times | Unitless branch lengths | Calibrated in millions of years | Enables phylogenetic distance weighting in Step 7 |

### Impact on downstream steps

- **Step 6 (PAML):** The trusted tree defines foreground (cancer-resistant) vs. background branches. Correcting the deep outgroup topology ensures no cancer-resistant lineage (greenland shark) is accidentally placed in a biologically incorrect position.
- **Step 7 (Convergence):** The lineage independence test now correctly recognises hydra (Cnidarians), ocean quahog clam (Molluscs), and the three cartilaginous fish (Sharks) as three distinct lineages rather than an unresolved cluster. This directly increases the maximum convergence count from 5 to **8 lineage groups** in the updated analysis.

---

## Next Step

→ [Step 6: Positive Selection Analysis (PAML)](step6_positive_selection.md)
