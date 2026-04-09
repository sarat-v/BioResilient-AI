# BioResilient Phase 2 — Final Check
## Continuity from Phase 1 + Readiness Assessment + Novelty Review

*Prepared April 2026. Reviews all Phase 2 steps against what Phase 1 produced and what the literature says is best practice.*

---

## Part 1: What Phase 1 Hands to Phase 2

Phase 1 produced a richly annotated database. Each of the 14 Tier1 and 37 Tier2 candidates carries the following already-computed signals:

| Evidence Layer | Where Stored | What It Captures |
|---|---|---|
| `EvolutionScore.selection_score` | PAML branch-site ω + LRT p-value | Positive selection pressure in resilient lineages |
| `EvolutionScore.convergence_score` | Permutation-tested convergence weight | Independent molecular convergence across 8 lineages |
| `EvolutionScore.convergence_pval` | 200-iteration label shuffle p-value | Statistical confidence: is convergence above chance? |
| `EvolutionScore.convergence_count` | Raw lineage count (1–8) | How many independent lineages show convergence |
| `DivergentMotif.consequence_score` | AlphaMissense mean score | How functionally significant are the divergent residues |
| `DivergentMotif.domain_name` | Pfam domain annotation | Which functional domain the divergence falls in |
| `DivergentMotif.variant_direction` | GoF / LoF / neutral | Predicted functional direction of change |
| `ConvergentAA.convergent_aa` | Specific amino acid | The exact convergent residue across lineages |
| `ConvergentAA.position` | Protein position (1-indexed) | The exact position of convergence |
| `EvolutionScore.dnds_ratio` | PAML ω value | How strong is positive selection (ω > 1 = selection) |
| `CandidateScore.expression_score` | GTEx/Bgee tissue expression | Baseline functional relevance |
| `CandidateScore.composite_score` | Phase 1 rank-product | Overall Phase 1 ranking |
| `CandidateScore.tier` | Tier1/Tier2/Tier3 | Current classification |

**Phase 2 funnel:** Steps 11–14 run on Tier1+Tier2 only (51 genes). Step 13 (gene therapy) is Tier1-only (14 genes). This is correct and computationally sound.

**Critical observation**: `ConvergentAA.position` contains the exact protein coordinates of the convergent changes. This is a gold mine that Phase 2 does not yet fully exploit. The specific residue-level data flows from Step 7b but is never re-queried in the druggability or disease annotation steps — they treat the gene as a unit, not the specific positions.

---

## Part 2: Step-by-Step Phase 2 Assessment

### Step 10b — AlphaGenome Regulatory Divergence ✅ FULLY CODED

**What it does:** Track B — for genes with high expression divergence but not captured by protein-level Track A, queries promoter sequences (±2 kb from TSS) across resilient species via the AlphaGenome API and scores regulatory effect magnitude.

**Scientific validity:** Sound. Regulatory changes can drive convergent phenotypes without protein-coding changes. This catches a class of targets that protein-only analysis would miss entirely.

**Status:** Fully implemented. Needs `ALPHAGENOME_API_KEY`. Silently skips if absent.

**Gap:** Score formula (`min(max_effect + 0.1 × lineage_count, 1.0)`) is reasonable but could be improved by weighting lineage count more strongly. Right now a gene with 1 highly divergent species scores almost as high as one with 5 moderately divergent lineages.

**Recommendation:** Run this. Set the API key. The output feeds `RegulatoryDivergence` table → `regulatory` weight (5%) in Phase 2 scoring.

---

### Step 11 — Disease Annotation ⚠️ CODED BUT UNDERPERFORMING

**What it does:** Runs OpenTargets, GWAS Catalog, gnomAD, IMPC, Human Protein Atlas, and pathway annotation for all Tier1+Tier2 genes.

**What's good:**
- OpenTargets: queries `associatedDiseases` scores → `DiseaseAnnotation.opentargets_score`
- GWAS Catalog: fetches associations → `gwas_pvalue` stored
- gnomAD: fetches `pLI` (LoF intolerance) → `gnomad_pli`
- IMPC: mouse knockout phenotypes
- HPA: tissue specificity
- Pathway GO annotation

**What's missing — OpenTargets is being used at ~30% of its potential:**

The current query fetches `associatedDiseases` scores (generic phenotype association). It does NOT fetch:
1. **Tractability data** — whether the target is small-molecule/antibody/PROTAC tractable (already in the platform, just not queried)
2. **L2G (Locus-to-Gene) score** — the ML model that identifies if our gene is *causal* at GWAS loci, not just nearby. This is 8.1x enriched for known drug targets.
3. **Safety data** — adverse event flags from OpenTargets (HPA subcellular location for off-target concern, known toxicity annotations)
4. **Known drugs** — whether any approved drug already targets this gene (fastest path to repurposing)

**Recommended query addition for opentargets.py:**

```graphql
query TargetAnnotation($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    tractability {
      label
      modality
      value
    }
    knownDrugs {
      count
      rows {
        drug { name maximumClinicalTrialPhase }
        disease { name }
        phase
      }
    }
    safetyLiabilities { event effects { tissue { label } } }
    associatedDiseases {
      rows { disease { id name } score datatypeScores { id score } }
    }
  }
}
```

Adding tractability + known drugs to the query takes 2 hours of coding and dramatically improves the disease annotation layer. **Priority: HIGH before running.**

---

### Step 11b — Rare Protective Variant Mapping ✅ EXCELLENT — DO NOT CHANGE

**What it does:** For each DivergentMotif, maps its protein positions to hg38 genomic coordinates, queries gnomAD v4 for rare variants (MAF < 1%) at those exact positions, checks whether the rare variant introduces the *same amino acid* seen in the resilient species (or Miyata-equivalent), then queries GWAS Catalog for phenotype associations on those variants.

**Scientific validity:** This is the PCSK9 paradigm applied directly to our convergent animal positions. If a human individual happens to carry a rare variant that mimics what the bowhead whale has at that exact protein position, and that person is healthier — that is as close to experimental validation as you can get computationally.

**Why it's genuinely novel:** No published pipeline does this. The typical approach is to find GWAS associations near a gene. We are finding rare variants at the *exact convergent residues* identified by cross-species evolution. The specificity is an order of magnitude higher.

**Miyata biochemical similarity groups**: Well implemented. Checks not just identical amino acids but biochemically similar ones (same charge, polarity class).

**What could be improved:** Currently upgrades to "Validated" tier if `protective_variant_count >= 1 AND protective_variant_pvalue < 5e-8`. Consider also logging whether the direction of the human rare variant matches the animal change (gain-of-function vs. loss-of-function, from Step 4d `variant_direction`).

---

### Step 11c — Literature Validation (PubMed) ✅ CODED

**What it does:** Searches PubMed for each gene symbol combined with longevity/resilience/cancer-resistance terms.

**Status:** Coded. This is a sanity check, not a primary scoring layer. Good for the report.

---

### Step 11d — Pathway Convergence Scoring ✅ EXCELLENT — DO NOT CHANGE

**What it does:** Maps all Tier1+Tier2 genes to Reactome pathways, computes hypergeometric enrichment (how many of our convergent genes are in the same pathway, above chance?), and weights by the sum of convergence scores of pathway members.

**Scientific validity:** This is the right way to handle it. If 4 of our 14 Tier1 genes are all in the PI3K/AKT pathway, that is not coincidence — it is evidence that the whole pathway is under selection pressure, not just the individual genes. The hypergeometric test on top of evolutionary weighting is scientifically rigorous.

**Why it's novel:** Published convergent evolution papers identify individual genes. We identify *pathways* that are under convergent selection pressure — this is a higher-order signal that strengthens the whole analysis.

**What it outputs:** `PathwayConvergence` table → exposed at `GET /research/pathway-convergence`. This should be prominently featured in the final report.

---

### Step 12 — Druggability Assessment ⚠️ CODED BUT MISSING MODALITY BREADTH

**What it does:**
- `structure.py`: Downloads AlphaFold structures for each gene
- `pockets.py`: Runs fpocket on the structure → `DrugTarget.pocket_count`, `top_pocket_score`
- `chembl.py`: Queries ChEMBL for known drugs/compounds → ChEMBL activity score
- `cansar.py`: Queries CanSAR for druggability tiers
- `peptide.py`: Peptide tractability for divergent motifs

**Step 12b:** P2Rank ML pocket prediction on top of fpocket → improves recall by 14% per benchmark.

**What's good:** fpocket + P2Rank combination is the 2024 benchmark winner. Structure downloads from AlphaFold. ChEMBL integration.

**What's missing:**

1. **PROTAC/degrader tractability**: The PROTACtable genome criteria (protein stability, E3 ligase accessibility, intracellular location) are not assessed. Over 1,000 proteins not tractable by small molecules ARE tractable by PROTACs. For Tier1 genes with poor SM pocket scores, PROTAC should be checked automatically.

2. **Convergent residue proximity to pockets**: fpocket finds pockets on the whole protein. It never checks whether any of our specific convergent positions (from `ConvergentAA` table) are *near* those pockets. A convergent residue within 6Å of the top pocket is far more interesting than a convergent residue in a disordered loop on the opposite face.

3. **Antibody tractability**: Whether the protein has extracellular epitopes (for antibody drugs) is not assessed. Simple UniProt subcellular location lookup.

4. **SGC probe check**: The Structural Genomics Consortium has chemical probes for many proteins — a quick check of the SGC portal API could identify immediate experimental starting points.

**Recommended addition (pocket proximity to convergent residues):**
After fpocket + P2Rank run, for each Tier1 gene, compute the minimum distance between convergent positions (`ConvergentAA.position`) and the top pocket's residue list (from fpocket output). Flag as "convergent-pocket-proximal" if within 6Å. This is the single most valuable addition to druggability and could be done in ~1 day of coding.

---

### Step 13 — Gene Therapy Feasibility ✅ CODED — REASONABLE SCOPE

**What it does:**
- `aav.py`: Gene size from NCBI, checks against AAV 4.7 kb limit, recommends serotype based on GTEx tissue expression
- `crispr.py`: Guide RNA design and off-target assessment (need to verify tool used)

**What's good:** Gene size check is crucial — large genes are immediately flagged as AAV-incompatible. Tissue-to-serotype mapping is scientifically correct.

**What it needs:**
1. This step should use **Step 9b output** (convergent residue structural context) to specify *what edit* is being designed. Currently it assesses delivery feasibility generically for the whole gene rather than for a specific edit at a specific residue.
2. The CRISPR module needs GuideScan2 or similar for off-target scoring — verify what tool `crispr.py` actually uses.

**Overall**: This step is appropriate at Phase 2 scope. A feasibility scorecard (AAV-compatible? CRISPR guide quality? Preferred tissue?) is the right output — not a full therapy design.

---

### Step 14 — Safety Pre-Screen ⚠️ CODED BUT PheWAS NEEDS UPGRADE

**What it does:**
- `phewas.py`: Queries GWAS Catalog for all phenotype associations on the gene → identifies unexpected disease links (safety flags)
- `network.py`: STRING network hub score → flag if gene has many interaction partners (off-target risk)
- `selectivity.py`: Expression breadth / tissue selectivity

**Step 14b:**
- `depmap.py`: DepMap CRISPR essentiality score
- `gtex.py`: GTEx expression breadth (how many tissues express this gene)

**What's good:** The multi-dimensional safety approach is right. STRING network hub analysis is a good proxy for systemic risk. DepMap essentiality is the right question for cancer targets.

**Critical issue — PheWAS is using naive GWAS Catalog gene-level query:**

The current `phewas.py` queries GWAS Catalog for ALL associations on a gene symbol and returns them as safety flags. The problem: many of these associations are driven by LD — a nearby SNP is associated with a trait, not the gene itself. This creates false safety signals for genes near GWAS loci.

**Fix:** Add a simple filter: only flag associations where the GWAS lead SNP is within the gene body (CDS ± 1 kb), not just "gene nearby." This can be done by comparing genomic coordinates. A more rigorous fix would use CoPheScan (R-based, Bayesian LD correction), but the coordinate filter is achievable in a day.

**DepMap essentiality nuance:** Currently using raw essentiality score. The correct metric is **selective essentiality** — essential in cancer cells but not normal cells. The DepMap portal provides context dependency data. The current implementation needs to pull tumor vs. normal context separately.

---

### Step 15 — Final Rescoring ⚠️ WEIGHTS NEED SAFETY CORRECTION

**What it does:** Re-runs `run_scoring(phase="phase2")` which applies the Phase 2 weight configuration, then applies the control species divergence penalty.

**Phase 2 weights (from `config/scoring_weights.json`):**

| Evidence Layer | Weight |
|---|---|
| convergence | 0.22 |
| selection | 0.18 |
| disease | 0.20 |
| druggability | 0.15 |
| expression | 0.10 |
| safety | 0.10 |
| regulatory | 0.05 |

**What's good:** Evolutionary signal remains dominant at 40% combined (convergence + selection). Disease annotation at 20% is appropriately significant. The regulatory layer at 5% picks up AlphaGenome Track B.

**Critical issue — safety scoring is additive when it should be a floor/penalty:**

Currently: `composite = 0.22×conv + 0.18×sel + 0.20×dis + 0.15×drug + 0.10×expr + 0.10×safety + 0.05×reg`

A gene with perfect evolutionary signal (conv=1, sel=1) and maximum druggability (drug=1) but catastrophic safety (safety=0, e.g., expressed in the heart as the primary tissue, hERG channel binder) would score: `0.22 + 0.18 + 0.20 + 0.15 + 0.10 + 0 + 0 = 0.85` — Tier1.

That is wrong. A target with a fatal safety liability should not reach Tier1 regardless of its evolutionary signal.

**Recommended fix in `scoring_weights.json`:**

```json
"phase2": {
  "convergence":  0.25,
  "selection":    0.20,
  "disease":      0.22,
  "druggability": 0.18,
  "expression":   0.10,
  "regulatory":   0.05,
  "safety":       0.00,
  "safety_floor": 0.40
}
```

Then in `scoring.py`, Phase 2 composite = `weighted_sum × safety_multiplier` where `safety_multiplier = safety_score` if `safety_score > safety_floor`, else `0`. A gene with `safety_score < 0.40` gets zero regardless of evolutionary evidence. This is the industry-standard approach.

**Control species penalty timing:** Currently at the end of Step 15. Correct for the first Phase 2 run since all data is fresh. After the first run, this penalty could be pre-applied, but it's fine here for now.

---

## Part 3: Missing Step — Step 9b Structural Annotation

This is the most impactful thing not yet coded. All the evidence we have from Phase 1 is gene-level or motif-level. We know:
- **Gene X** shows convergent evolution in 6 lineages
- **Residue 247** of Gene X is the convergent position
- **AlphaMissense score** at position 247 is 0.82 (likely pathogenic if changed in humans)

What we do NOT know:
- Is residue 247 in the active site, the allosteric cleft, a protein-protein interface, or a disordered loop?
- Is the top fpocket pocket adjacent to residue 247?
- Does the convergent change face into the pocket or away from it?

Without answering this, the druggability step runs on the whole protein and reports generic pockets. Any chemist reviewing the output would immediately ask "yes, but WHERE is the convergent residue relative to the drug binding site?" and we couldn't answer.

**What Step 9b would add:**
- SIFTS residue annotation: is position 247 in an InterPro functional site? A Pfam catalytic residue? A known binding site?
- fpocket residue list: is position 247 in the top pocket's residue set (within 6Å)?
- AlphaMissense contextual classification: is the convergent change classified as pathogenic, benign, or ambiguous?
- Structural label: one of {active_site, allosteric, interface, pocket_proximal, surface, buried, disordered}

**Estimated coding time:** 2–3 days. The data is already available (SIFTS REST API, AlphaMissense TSV already downloaded, fpocket output already computed in Step 12). This is a query and join step, not a new computation.

**When to add it:** Before running Phase 2 for real. Can be added in parallel with Phase 2 infrastructure fixes.

---

## Part 4: Novelty Assessment — Will We Find Something New?

### Evidence Layers That Are Genuinely Novel (Not Published Elsewhere)

**1. The 8-lineage permutation-tested convergence approach (Phase 1, Step 7)**

No published pipeline combines:
- 8 independent evolutionary lineages including Molluscs and Cnidarians (most use 3–5 mammalian lineages)
- Permutation-based null model to assign p-values to convergence (most use raw lineage counts)
- TimeTree-calibrated phylogenetic distance weighting

The BMC Genomics 2023 paper on mTOR found 2 convergent genes in the mTOR network using 6 mammalian species. We found mTOR pathway members (AKT1) with 8-lineage convergence using a broader phylogenetic span. That is an improvement.

**2. Rare protective variant mapping at convergent positions (Step 11b)**

This is the most novel piece. The PCSK9 paradigm is well-known but is never applied at residue-specific resolution from a cross-species convergence signal. We are asking: "does the exact residue that whales and naked mole rats independently evolved ALSO vary in humans, and are people with that variation healthier?"

If even one of the 14 Tier1 genes has a rare human variant at the exact convergent position that associates with cancer resistance or longevity in a GWAS — that is a publication-worthy finding. The specificity of the claim is orders of magnitude higher than "this gene is near a GWAS hit."

**3. Pathway-level convergent selection enrichment (Step 11d)**

Hypergeometric test asking: are our convergent genes enriched in the same pathways above chance? Weighted by their evolutionary scores. This goes beyond gene-level analysis to pathway-level evolutionary pressure. If 4 of our 14 Tier1 genes are in the PI3K/AKT/mTOR pathway — that is not noise. It is the pathway being selected.

**4. Dual-track evolutionary evidence (Protein divergence + Regulatory divergence)**

Track A (Steps 1–9): protein sequence convergence
Track B (Step 10b, AlphaGenome): regulatory/promoter sequence divergence

A gene that shows BOTH protein-level convergence AND promoter-level divergence is extremely compelling. No other published pipeline uses both simultaneously for the same phenotype.

**5. The full evidence stack on the same 14 candidates:**

When Phase 2 is complete, each Tier1 gene will have:
- Evolutionary convergence score with p-value (cross-species)
- PAML positive selection score (within-lineage)
- AlphaMissense consequence score at convergent residues (functional impact)
- ESM-2 structural variant embedding (structural plausibility)
- Pfam domain context (which functional domain)
- OpenTargets disease associations
- Rare protective variant mapping (human genetics validation)
- Druggability score (fpocket + P2Rank)
- Safety profile (PheWAS, DepMap, expression breadth)
- Gene therapy feasibility

No published computational target discovery pipeline generates all of these for the same set of evolution-nominated candidates. This is the novel contribution.

---

## Part 5: Before-You-Run Checklist

Ordered by priority:

### Must-Do Before Running (1–3 days)

**1. Upgrade OpenTargets query to include tractability + known drugs (Step 11)**
- Add tractability (SM/AB/PROTAC), known drugs (fastest path to repurposing), and safety liabilities to the GraphQL query
- File: `pipeline/layer3_disease/opentargets.py`
- Effort: ~4 hours

**2. Fix safety scoring from additive to multiplicative floor (Step 15)**
- Add `safety_floor` parameter to `config/scoring_weights.json`
- Modify `scoring.py` Phase 2 composite to multiply by safety factor rather than add
- File: `pipeline/scoring.py`, `config/scoring_weights.json`
- Effort: ~2 hours

**3. Add convergent-residue-to-pocket proximity check (Step 12)**
- After fpocket + P2Rank run, cross-reference `ConvergentAA.position` against fpocket's pocket residue list
- Flag genes where convergent residues are within 6Å of the top pocket
- File: `pipeline/layer4_druggability/pockets.py`
- Effort: ~1 day

### Should-Do Before Running (3–5 days)

**4. Set ALPHAGENOME_API_KEY**
- Register at https://deepmind.google/technologies/alphagenome/
- Add to environment config
- Effort: 30 minutes (account registration) + testing

**5. Add gnomAD v4.0 pLI/LOEUF upgrade to gnomad.py (Step 11)**
- Current code uses gnomAD v2/v3 endpoint; v4.0 has 6x more samples and updated LOEUF threshold (< 0.6 replaces < 0.35)
- File: `pipeline/layer3_disease/gnomad.py`
- Effort: ~3 hours

**6. Add coordinate filter to PheWAS to reduce LD-driven false safety flags (Step 14)**
- Only flag associations where lead SNP falls within gene body ±1 kb
- File: `pipeline/layer6_safety/phewas.py`
- Effort: ~4 hours

### Plan for Phase 2.1 (After First Phase 2 Run)

**7. Step 9b — Structural Annotation**
- SIFTS API + convergent position functional context + pocket proximity
- Estimated: 2–3 days
- Value: Transforms druggability from gene-level to residue-level assessment

**8. SharePro colocalization for any Tier1 gene with GWAS association**
- Confirms whether our gene is causal at the GWAS locus or just nearby
- Estimated: 1 day (R script, runs on already-collected GWAS data)

---

## Part 6: Summary Verdict

| Dimension | Assessment |
|---|---|
| **Scientific Foundation** | Solid. 14 Tier1 candidates from 12,795 genes with rigorous permutation-tested convergence. |
| **Phase 2 Coding Completeness** | ~80% coded. Key gaps: OpenTargets tractability, safety floor, pocket-proximity check. |
| **Novelty** | High. The rare variant mapping at convergent positions (Step 11b) and pathway convergence enrichment (Step 11d) are not in published pipelines. |
| **Biggest Gap** | Step 9b (structural annotation) — residue-to-structure mapping is missing but feasible. |
| **Most Critical Fix** | Safety scoring must become a multiplicative floor, not an additive term. |
| **Novel Discovery Potential** | High. If one Tier1 gene has a human rare variant at the exact convergent residue with phenotypic association, that is a publishable discovery. |
| **Readiness to Run** | Run after completing the 3 must-do fixes above. Can run current version for a draft result. |

The pipeline is architecturally sound and scientifically ambitious. The combination of evolutionary evidence, human genetic validation, and functional genomics creates a multi-modal evidence stack that is genuinely different from existing target discovery approaches. The three must-do fixes are the difference between a good pipeline and a great one.
