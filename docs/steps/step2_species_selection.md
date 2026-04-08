# Step 2 — Species Selection & Proteome Assembly

**Pipeline phase:** Data preparation  
**Run timestamp:** 2026-03-14 08:54:27 UTC  
**Status:** ✅ PASS

---

## What This Step Does

Step 2 defines the comparative species panel and assembles their reference proteomes. The fundamental hypothesis of BioResilient is:

> *Genes under convergent positive selection in cancer-resistant species are mechanistically linked to the cancer-resistance phenotype.*

This requires a carefully designed species set with three properties:
1. **Independent origins of cancer resistance** — multiple lineages evolved resistance separately, so shared genomic signals are unlikely to be coincidental
2. **Sufficient phylogenetic spread** — species should not all be closely related, otherwise signals reflect shared ancestry rather than convergent adaptation
3. **Control species** — phylogenetically similar but cancer-susceptible species allow the pipeline to subtract out signals that are simply due to normal mammalian evolution

---

## Species Panel — 18 Species

### Experimental Species (cancer-resistant / long-lived)

| Species | Lineage | Max Lifespan | Cancer-Resistance Mechanism |
|---|---|---|---|
| **Naked mole rat** | Rodents | 37 yr | Near-zero cancer incidence; HMW-HA (high-MW hyaluronic acid) induces early contact inhibition; p16/p27 dual checkpoint |
| **Blind mole rat** | Rodents | 21 yr | IFN-β-mediated concerted cell death kills all cells in a clone once colony reaches ~20 cells; no cancer ever reported |
| **Damaraland mole rat** | Rodents | 28 yr | Longest-lived rodent relative to size; negligible senescence |
| **Beaver** | Rodents | 24 yr | Long-lived for a rodent; low documented cancer rates |
| **Bowhead whale** | Cetaceans | 200+ yr | Extreme lifespan and body size (Peto's paradox); somatic mutation rate 10× lower than expected |
| **Sperm whale** | Cetaceans | 70+ yr | Large, long-lived; expanded DNA repair gene families |
| **African elephant** | Proboscideans | 70 yr | 20 copies of *TP53* (humans have 2); enhanced apoptotic response to DNA damage |
| **Asian elephant** | Proboscideans | 70 yr | Same *TP53* amplification mechanism as African elephant |
| **Elephant shark** | Cartilaginous fish | ~100 yr | Slowest-evolving jawed vertebrate genome; ~400 My divergence from human |
| **Little skate** | Sharks | Long-lived | Chondrichthyes outgroup; minimal genome instability |
| **Greenland shark** | Sharks | 400+ yr | Longest-lived vertebrate known; near-zero reported cancer |
| **Little brown bat** | Bats | 34 yr | Exceptional longevity relative to body mass; strong antiviral/DNA damage responses |
| **Painted turtle** | Reptiles | 50+ yr | Long-lived ectotherm; tolerant of extreme hypoxia |
| **Hydra** | Cnidarians | Theoretically immortal | Continuous stem cell renewal; telomere maintenance without telomere shortening |
| **Ocean quahog clam** | Molluscs | 507+ yr | Longest-lived non-colonial animal; Ming the clam (died 2006 during sampling, age 507) |

### Control Species (cancer-susceptible)

| Species | Lineage | Role |
|---|---|---|
| **Rat** | Rodents | Closely related to mole rats and beaver; short-lived (2–3 yr); high cancer rates. Allows isolation of rodent cancer-resistance signals vs. general rodent biology |
| **Macaque** | Primates | Closely related to human; cancer-susceptible. Allows filtering of primate-specific signals |

**Human** is the query species — all evolutionary signals are ultimately mapped back to human gene coordinates.

---

## Proteome Assembly

For each species, reference proteomes are sourced from:
- **UniProt/RefSeq** for well-annotated species (human, macaque, rat, elephants, whales)
- **NCBI GenBank assemblies** for less-annotated species (sharks, clam, hydra)
- Species proteomes are stored in S3 and loaded into the pipeline database

---

## Outputs

### Database — species table

```
18 rows inserted into species table:
  - 16 experimental species (is_control = false)
  - 2 control species (is_control = true): macaque, rat
```

### Ortholog coverage by species (loaded in Step 3b)

| Species | Genes with orthologs | Coverage |
|---|---|---|
| Human | 12,795 | 100% (query) |
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

> The lower coverage in sharks (60%), clams (60%), and hydra (52%) is expected — these are the most divergent species and many mammalian genes either don't have clear orthologs in invertebrates or the genomes are less completely annotated. These species still contribute valuable evolutionary signal for highly conserved genes.

---

## Validation

| Check | Result |
|---|---|
| Species count | ✅ 18 species loaded |
| Proteomes accessible | ✅ All 18 proteome paths resolved in S3 |
| Control species flagged | ✅ macaque (is_control=true), rat (is_control=true) |

---

## Next Step

→ [Step 3: Ortholog Detection & Genomic Region Extraction](step3_ortholog_detection.md)
