# Step 13 — Gene Therapy Feasibility

**Pipeline phase:** Phase 2 — Clinical Translation Layer  
**Run timestamp:** 2026-04-11  
**Status:** ✅ PASS (informational only — not part of composite scoring)  
**Scope:** Phase 1 Tier 1 genes only (15 genes)  

---

## What This Step Does

Step 13 evaluates whether candidate genes are feasible targets for gene therapy — delivery of a corrective or modifying gene sequence into human cells. This is the most direct path from a cancer-resistance gene to a clinical intervention: if a gene confers protection in longevity model species, could we deliver a modified version of it to human cells?

**Important note on scope:** Step 13 runs on Phase 1 Tier 1 genes (those with strongest evolutionary signal before Phase 2 rescoring). After Phase 2 composite scoring, tier assignments shift; gene therapy data should be re-run on the final Tier 1 / Validated list in a future iteration. The gene therapy assessment does **not** contribute to the composite_score — it is informational only for therapeutic planning.

---

## Evidence Sources

### AAV Compatibility (aav.py)

Adeno-associated virus (AAV) is the most clinically advanced gene delivery vehicle. Key constraint: **AAV cargo capacity is ~4.7 kb**. Gene therapy by AAV is straightforward only for genes whose coding sequence (CDS) fits within this limit.

For each gene:
- Retrieve CDS length from NCBI Gene / RefSeq
- Flag `aav_compatible = True` if CDS ≤ 4,500 bp (conservative threshold, leaving room for promoter + ITRs)
- Identify matching AAV serotypes based on target tissue expression (from GTEx/HPA)

### CRISPR Feasibility (crispr.py)

For larger genes or when overexpression is not desired (base editing, epigenetic reprogramming):
- Count predicted CRISPR cut sites in the gene locus
- Assess off-target risk based on guide RNA uniqueness

---

## Results

| Metric | Value |
|---|---|
| Genes evaluated | 15 (Phase 1 Tier 1) |
| AAV-compatible | 15 / 15 |
| Average CDS length | 1,436 bp |
| CDS range | 750 – 2,334 bp |
| Genes with tissue tropism matched | 14 / 15 |
| Genes with CRISPR sites assessed | 14 / 15 |

All 15 evaluated genes fall well within the AAV size limit. This is not surprising — Phase 1 Tier 1 genes are enriched for well-studied proteins with compact, well-characterized coding sequences.

### Gene Therapy Candidates

| Gene | CDS (bp) | AAV | Tissue Tropism | Composite Score | Tier |
|---|---|---|---|---|---|
| G6PD_HUMAN | 1,845 | ✅ | Available | 0.7091 | Validated |
| FNTA_HUMAN | 1,437 | ✅ | Available | 0.4562 | Validated |
| AKT1_HUMAN | 1,740 | ✅ | Available | 0.4060 | Validated |

**AKT1_HUMAN** is flagged as essential (`is_essential = true`, depmap_score = 0.099) — full AAV delivery/overexpression requires careful dosing to avoid toxicity. Gene therapy approaches for AKT1 would focus on tissue-specific promoters and expression level control.

---

## Scientific Caveats

- **AAV compatibility is necessary but not sufficient**: CDS ≤ 4,500 bp is a prerequisite, but AAV-based gene therapy also requires tissue-specific tropism, manufacturing scalability, immune tolerance (pre-existing anti-AAV antibodies in ~50% of the population for common serotypes), and appropriate promoter design. The pipeline assesses only cargo size and approximate tissue tropism.
- **100% AAV compatibility for 15 evaluated genes**: this is biased towards the Phase 1 Tier 1 gene set, which over-represents well-studied proteins with compact sequences. The full Tier1/2/Validated set (155 genes) likely includes some genes exceeding the 4.5 kb threshold.
- **Step 13 is informational only**: gene therapy scores do **not** feed into the composite_score. They are displayed in the API and dashboard as a separate therapeutic planning layer.

## Limitations and Future Work

1. **Coverage gap**: Only 15 of 155 final Tier 1/2/Validated genes were evaluated (Phase 1 Tier 1 only). A rerun of Step 13 on all 155 Validated/Tier1/Tier2 genes from Phase 2 is **recommended before experimental prioritisation** to ensure gene therapy data reflects the final tier assignments.
2. **Base editing**: The pipeline does not yet model base-editing approaches (e.g., adenine base editors, ABE8e/ABE8.20) that could directly introduce the longevity-model amino acid at convergent positions without delivering a full gene. This is scientifically highly relevant as it would enable single-nucleotide precision to recapitulate the animal adaptation.
3. **mRNA therapy**: No assessment of mRNA delivery feasibility (relevant for transient expression strategies; e.g., liver-tropic LNP delivery for secreted proteins or enzymes).
4. **Epigenetic modulation**: No assessment of whether the gene's expression level (rather than protein sequence) could be therapeutically modulated via HDAC inhibitors, CRISPRa, or dCas9-based systems — relevant for genes where increased expression of the wild-type protein is protective.

---

## Data Quality Assessment

✅ All CDS lengths from NCBI RefSeq (verified against UniProt sequence lengths)  
✅ aav_compatible Boolean consistent with gene_size_bp values  
ℹ️ Not part of composite scoring — no downstream scoring impact  
