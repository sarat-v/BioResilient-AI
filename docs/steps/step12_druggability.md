# Step 12 — Druggability Assessment (fpocket + ChEMBL + P2Rank)

**Pipeline phase:** Phase 2 — Clinical Translation Layer  
**Run timestamp:** 2026-04-11  
**Status:** ✅ PASS — 279 / 280 genes with pocket predictions; 279 with P2Rank ML scores  

---

## What This Step Does

Step 12 assesses whether each candidate gene encodes a protein that could be targeted by a small molecule, antibody, or degrader. Evolutionary evidence tells us *which* proteins differ between cancer-resistant and susceptible species; druggability assessment tells us *whether* those proteins are chemically tractable.

---

## Evidence Sources

### Step 12 — fpocket (Geometric Pocket Detection)

**fpocket** (Le Guilloux et al. 2009) uses Voronoi tessellation of the protein surface to detect cavities suitable for small-molecule binding. Inputs are AlphaFold PDB structures (filtered to pLDDT > 50 to exclude disordered regions).

| Metric | Value |
|---|---|
| Genes with pocket predictions | 279 / 280 |
| Average pocket count per gene | 37.7 |
| Average top pocket score | 0.657 |
| Genes with ChEMBL target ID | 65 |
| Genes with existing drugs | 5 (in ChEMBL database) |
| Druggability tier assigned | 69 genes |

**One gene missing pocket data** is due to no available AlphaFold structure (cannot run fpocket without 3D coordinates).

### Step 12b — P2Rank (Machine Learning Pocket Prediction)

**P2Rank** (Krivák & Hoksza 2018) is a random-forest model trained on known drug-binding sites. It scores protein surface patches using physicochemical and evolutionary features, independent of geometric cavity detection.

Using both fpocket (geometry-based) and P2Rank (ML-based) provides complementary pocket evidence — false positives from one method are unlikely to be false positives in the other.

| Metric | Value |
|---|---|
| Genes with P2Rank predictions | 279 / 280 |
| Average P2Rank score | 0.207 |
| Genes with convergent-pocket-proximal motifs | 152 / 280 (54%) |

**`convergent_pocket_proximal`** flags genes where at least one convergently-selected motif lies within the predicted top druggable pocket — the highest-confidence indicator that evolutionary pressure occurred at a druggable site.

---

## Tractability (from Open Targets, Step 11)

Open Targets Platform annotates each gene's tractability across three drug modalities:

| Modality | Genes tractable | Interpretation |
|---|---|---|
| Small molecule (SM) | 36 | Established chemical matter or druggable active site |
| Antibody (AB) | 64 | Accessible extracellular epitope or validated antibody |
| PROTAC / degrader | 76 | E3 ligase proximity-compatible structure |

PROTAC tractability is notably high (76 genes), consistent with many candidates being nuclear or cytoplasmic proteins where degradation strategies are being actively explored.

---

## Druggability Score Formula

```
druggability_score = 0.20 × pocket_count_norm          [fpocket pocket count, capped at 10]
                   + 0.20 × top_pocket_score_norm       [fpocket best pocket score]
                   + 0.10 × p2rank_score_norm           [P2Rank ML score]
                   + 0.15 × chembl_evidence             [ChEMBL target or drug]
                   + 0.10 × cansar_tier_score           [CanSAR druggability tier]
                   + 0.10 × tractability_bonus          [SM / AB / PROTAC flags]
                   + 0.15 × convergent_pocket_proximal  [convergent motif in top pocket]
```

---

## Top Druggability Scores Among Final Candidates

| Gene | Top Pocket | P2Rank | Conv-Pocket | Tractability | Drug Score |
|---|---|---|---|---|---|
| NR6A1_HUMAN | 0.988 | — | Yes | SM+AB | 0.693 |
| NPSR1_HUMAN | 0.988 | — | Yes | SM+AB+PROTAC | 0.829 |
| ABCB7_HUMAN | 0.986 | — | Yes | SM+AB+PROTAC | 0.823 |
| G6PD_HUMAN | 0.825 | — | Yes | SM+AB+PROTAC | 0.841 |
| CAB39_HUMAN | 0.836 | — | Yes | SM+AB+PROTAC | 0.618 |
| GSTK1_HUMAN | 0.624 | — | No | SM+PROTAC | 0.792 |

---

## Scientific Caveats

- **Druggability score weights (0.20 pocket count, 0.20 pocket score, 0.10 P2Rank, etc.) are heuristic** and not derived from a training set with known drug–target outcomes. They reflect the relative informative value of each evidence source as judged by the pipeline design team. The convergent-pocket-proximal term (0.15) is the most biologically motivated: it directly links evolutionary pressure to a druggable site.
- **ChEMBL coverage (65/280 genes = 23%)** means that 77% of candidates have no small-molecule history in ChEMBL. This is expected for a discovery-focused pipeline: novel genes by definition lack prior chemical matter. ChEMBL absence does not indicate undruggability.
- **P2Rank scores not shown in top druggability table**: the table above shows `—` for P2Rank because P2Rank scores were not cached at the time of the run summary. Scores are stored in `drug_target.p2rank_top_score` and are available in the database. Future table updates should include this column.
- **fpocket + P2Rank agreement** provides orthogonal validation of pocket quality. A gene where both methods identify the same top pocket (convergent_pocket_proximal = True for both) is a stronger druggability call than either method alone.

## Data Quality Assessment

✅ 280 drug_target rows, 0 orphaned  
✅ top_pocket_score in [0, 1]  
✅ convergent_pocket_proximal boolean consistent with CPA pocket_distance_angstrom  
✅ tractability flags from Open Targets API v4 (handles both legacy and current field naming)  
