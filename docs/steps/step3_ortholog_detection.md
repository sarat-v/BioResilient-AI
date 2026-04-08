# Step 3 — Ortholog Detection & Genomic Region Extraction

**Pipeline phase:** Sequence acquisition  
**Substeps:** 3a (OrthoFinder), 3b (DB load), 3c (nucleotide extraction), 3d (PhyloP conservation)  
**Run timestamps:** 3b: 2026-04-06 16:53:27 | 3c: 2026-04-06 16:53:33 | 3d: 2026-04-06 16:53:39 UTC  
**Status:** ✅ PASS (all substeps)

---

## What This Step Does

Step 3 is the foundation of all downstream evolutionary analysis. It answers: **for each human gene, which proteins in the other 17 species are its evolutionary equivalent?**

This is non-trivial because:
- Genes duplicate, fuse, split, and are lost across 500+ million years of evolution
- Simple BLAST top-hit assignment creates false ortholog calls at scale
- We need the *ortholog* (same gene, different species) not the *paralog* (duplicated copy in the same species)

---

## Step 3a — OrthoFinder Orthogroup Clustering

### Tool

**OrthoFinder v2.5** with **DIAMOND v2** as the pairwise aligner.

### Algorithm

1. **All-vs-all protein alignment:** DIAMOND aligns every protein against every other protein across all 18 proteomes (≈18² = 324 pairwise comparisons run in parallel on AWS Batch)
2. **Bit-score normalisation:** Raw DIAMOND scores are normalised by sequence length and self-score to produce comparable distances across proteins of different lengths
3. **Markov Clustering (MCL):** A graph where proteins are nodes and normalised bit-scores are edge weights is clustered using MCL at inflation parameter I=1.5. This groups proteins into **Orthogroups (OGs)** — sets of proteins sharing common ancestry
4. **Species tree inference:** OrthoFinder reconstructs a gene tree for each OG and reconciles it with the species tree to distinguish ortholog (speciation node) from paralog (duplication node) relationships
5. **1-to-1 ortholog flagging:** Pairs that are single-copy in both species and inferred as orthologous at the speciation event are flagged `is_one_to_one = true`

### Why Orthogroups instead of BLAST hits?

BLAST best-hit strategies systematically confuse orthologs with paralogs, especially in gene families (e.g., the CDK family has many members — you need the gene tree to know which shark CDK corresponds to CDK4 vs CDK6). OrthoFinder's MCL approach correctly groups entire gene families and then uses gene tree reconciliation to distinguish which members are orthologous.

---

## Step 3b — Database Loading Results

### Summary

| Metric | Value |
|---|---|
| Human genes with ≥1 ortholog | **12,795** |
| Total ortholog sequences | **200,979** |
| 1-to-1 ortholog pairs | **199,561** (99.3%) |
| Orthogroups (OGs) | **12,768** |
| OGs with exactly 18 members (all species) | 5,091 |
| OGs with exactly 17 members | 2,226 |
| OGs with exactly 15 members (median) | 1,825 |

The **99.3% one-to-one rate** is exceptionally high and reflects the quality of the species selection — these 18 species span major vertebrate divergence points without excessive recent gene duplications, producing clean 1:1 mappings ideal for comparative genomics.

### OG Size Distribution

| OG size (members) | OG count | Interpretation |
|---|---|---|
| 18 (all species) | 5,091 | Most conserved genes — present in every species |
| 17 | 2,226 | Missing in one species (often hydra or greenland shark) |
| 15–16 | 2,597 | Small losses in phylogenetically distant species |
| ≤14 | 2,854 | Gene loss in multiple species or incomplete genome assemblies |

### Ortholog Coverage by Species (genes with human orthologs)

| Species | Gene count | % of human genes |
|---|---|---|
| Human | 12,795 | 100% |
| Macaque | 12,649 | 98.9% |
| Asian elephant | 12,321 | 96.3% |
| Bowhead whale | 12,293 | 96.1% |
| Blind mole rat | 12,241 | 95.7% |
| Sperm whale | 12,237 | 95.6% |
| Rat | 12,219 | 95.5% |
| African elephant | 12,177 | 95.2% |
| Naked mole rat | 12,162 | 95.1% |
| Beaver | 12,126 | 94.8% |
| Damaraland mole rat | 12,014 | 93.9% |
| Little brown bat | 11,682 | 91.3% |
| Painted turtle | 11,360 | 88.8% |
| Little skate | 10,402 | 81.3% |
| Elephant shark | 10,260 | 80.2% |
| Ocean quahog clam | 7,715 | 60.3% |
| Greenland shark | 7,709 | 60.3% |
| Hydra | 6,617 | 51.7% |

---

## Step 3c — Nucleotide Region Extraction

### What is extracted

For every gene–species ortholog pair, three genomic sequence regions are retrieved from **Ensembl REST API** (and equivalent databases for non-Ensembl species):

| Region | Size | Purpose |
|---|---|---|
| **CDS** | Variable (exons only) | Protein-coding sequence for codon-level analyses (PAML, dN/dS) |
| **Promoter** | 2,000 bp upstream of TSS | Core regulatory region; transcription factor binding sites |
| **Downstream** | 1,000 bp past last exon | Post-transcriptional regulatory elements; polyadenylation signals |

These regions feed into two parallel analyses:
- **CDS** → PAML positive selection (Step 6), codon alignments
- **Promoter + Downstream** → Regulatory divergence scoring (Step 3d), motif analysis (Step 4)

### Results

| Region | Sequences extracted | Genes with sequences |
|---|---|---|
| CDS | 183,622 | 12,784 |
| Promoter | 183,606 | 12,784 |
| Downstream | 183,613 | 12,784 |
| **Total** | **550,841** | **12,784** |

### Nucleotide Score Summary (pairwise vs. human)

Conservation scores computed for each gene–region pair by global pairwise alignment (MUSCLE + identity scoring):

| Region | Genes scored | Avg conservation | Avg % identity |
|---|---|---|---|
| CDS | 5,637 | 7.7% | 83.7% |
| Promoter | 3,583 | 11.1% | 82.2% |
| Downstream | 1,569 | 20.9% | 82.3% |

> The `conservation_score` here measures *variance* in conservation — low average values indicate most CDS positions are highly variable across the 18 species (expected: even conserved proteins vary in faster-evolving lineages like sharks and hydra). The `percent_identity` of ~82–84% across all regions reflects the average pairwise similarity between human and a given species ortholog.

**Regulatory divergence events (pairwise, all species vs. human):**
- Total divergence events: **1,242,496**
- Total regulatory convergence events: **1,241,969**

**Top genes by promoter divergence (most changed across species):**

| Gene | Promoter divergence events | Biological note |
|---|---|---|
| IL32_HUMAN | 2,165 | Interleukin 32 — immune signalling |
| TUG1_HUMAN | 2,080 | lncRNA involved in cell viability |
| OR6M1_HUMAN | 1,454 | Olfactory receptor — rapidly evolving |
| DCDC1_HUMAN | 1,426 | Doublecortin domain — neuronal migration |
| IL26_HUMAN | 1,408 | Interleukin 26 — immune gene |

---

## Step 3d — Phylogenetic Conservation Scoring (PhyloP)

### Tool

**PhyloP** scores from the UCSC 100-vertebrate multiple alignment track. PhyloP uses a phylogenetic model of neutral evolution to test each base for conservation or acceleration relative to the neutral rate.

### Score interpretation

- **PhyloP > 0** (positive): Slower evolution than neutral — purifying selection (site is functionally constrained)
- **PhyloP = 0**: Neutral evolution
- **PhyloP < 0** (negative): Faster evolution than neutral — potentially under positive selection or in a non-functional region

For gene-level summaries, PhyloP scores are summed over all bases in the region.

### Results

| Region | Genes scored | Mean PhyloP score | Interpretation |
|---|---|---|---|
| CDS | 11,400 | 0.01 avg | CDS is broadly conserved but varies widely |
| Promoter | 12,546 | 0.004 | Promoters evolve near-neutrally on average |
| Downstream | 12,546 | 717.4 (summed) | Downstream regions show moderate constraint |

**CDS PhyloP range:** 0.00 – 0.47 (gene-level averages)

> **Why low CDS PhyloP?** The gene-level average PhyloP is low because we are averaging over all sites including fast-evolving surface residues. The important signal is at the *per-site* level — conserved active site residues will have very high PhyloP, while surface loops will be near-zero or negative. The motif analysis in Step 4 works at the per-site level.

**Genes with accelerated promoter evolution (PhyloP < 0 across promoter):** 0 detected at gene level. This likely reflects that single-gene promoter acceleration is difficult to detect with the 100-vertebrate track since the signal is diluted by the many vertebrate species not in our panel.

---

## Database State After Step 3

```
gene table:           12,795 rows
ortholog table:      200,979 rows  (12,768 OGs × avg 15.7 sequences)
nucleotide_region:   550,841 rows  (3 regions × ~12,784 genes × avg 14.4 species)
nucleotide_score:     10,789 rows  (conservation scores per gene-region)
phylo_conservation:  12,546 rows
```

---

## Next Step

→ [Step 4: Protein Divergence Motif Analysis](step4_motif_analysis.md)
