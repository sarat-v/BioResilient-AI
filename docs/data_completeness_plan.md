# Data Completeness & Targeted Rerun Plan

## Audit Date: 2026-04-07

---

## 1. Current State of All DB Tables

### evolution_score (the scoring backbone)

| Category | Genes | Situation | In Rank Product |
|---|---|---|---|
| `paml_branch_site` real pval | 940 | PAML ran, detected selection signal | ✅ Real p-value used |
| `paml_branch_site` pval=1.0 | 6,018 | PAML ran, LRT not significant | ✅ Correctly neutral |
| `paml_no_signal` | 4,324 | PAML ran, confirmed no signal | ✅ Correctly neutral |
| `proxy` (HyPhy legacy) | 870 | No PAML run, old HyPhy fake p-values | ✅ Neutralized to 1.0 |
| `NULL` model | 381 | No selection model stored — silent gap | ⚠️ Treated as no signal |
| **Total unique genes** | **12,533** | | |

> **940 genes drive all positive selection signal currently.** The rest compete on convergence + functional evidence only.

### Other Tables — Clean Status

| Table | Rows | Genes | Status |
|---|---|---|---|
| `ortholog` | 200,979 | 12,795 | ✅ Complete |
| `divergent_motif` | 1,692,443 | via orthologs | ✅ Complete (ESM1v scored) |
| `phylo_conservation_score` | 12,546 | 12,546 | ✅ Near-complete |
| `nucleotide_region` | — | 12,795 | ✅ All genes have CDS sequences |
| `nucleotide_score` | 5,680 | 5,680 | ⚠️ Only 44% gene coverage |
| `expression_result` (GEO/DESeq2) | 10,905 | 5,716 | ✅ Step 5 output, not stale |
| `candidate_score` | 12,795 | 12,795 | 🔄 Being computed now (step 9) |
| `gene_therapy_score` | 0 | 0 | 🔴 Step not yet implemented |
| `pathway_convergence` | 0 | 0 | 🔴 Step not yet implemented |
| `regulatory_divergence` | 0 | 0 | 🔴 Step not yet implemented |

---

## 2. Issues Found

### Issue 1 — 381 NULL-model genes (HIGH PRIORITY)
**What:** 381 genes have a row in `evolution_score` with `selection_model IS NULL` and all scoring columns NULL.  
**Why it matters:** These genes HAVE aligned CDS sequences in `nucleotide_region` (avg 12.8 species/gene, 1,206 of 1,251 have 4+ species — sufficient for PAML). They were silently skipped at step 6.  
**Current impact:** Treated as no selection signal. They still score on convergence + functional evidence.  
**Action:** Targeted PAML rerun at step 6.

### Issue 2 — 870 proxy genes without PAML (HIGH PRIORITY)
**What:** 870 genes were processed by the old HyPhy pipeline. They were never submitted to PAML (step 6 in the current pipeline).  
**Why it matters:** Their old HyPhy p-values were artificial (hardcoded floor, not LRT results). These are now neutralized, but they also have full `nucleotide_region` CDS sequences (all 870 confirmed) — PAML can run on them.  
**Current impact:** No selection signal. 857/870 still have convergence scores (phylop avg 3.9) so they can reach Tier 1/2 on convergence alone.  
**Action:** Targeted PAML rerun at step 6 for their OGs. After rerun, replace proxy rows with PAML results.

### Issue 3 — 187 suspicious PAML entries (MEDIUM PRIORITY)
**What:** 187 genes have `selection_model='paml_branch_site'`, `dnds_ratio=99` (PAML's omega cap), and `dnds_pvalue=1.0`.  
**Why it matters:** omega=99 with pval=1.0 typically means PAML's LRT produced a degenerate result — either the null and alternative model converged to the same likelihood, or synonymous sites were too few to compute a valid LRT statistic.  
**Current impact:** Treated as neutral (pval=1.0). If any of these have genuine positive selection, it is being missed.  
**Action:** Rerun PAML for these 187 OGs with stricter convergence settings (or filter at scoring time based on synonymous site count).

### Issue 4 — nucleotide_score covers only 44% of genes (LOW PRIORITY)
**What:** `nucleotide_score` has 5,680 gene rows vs 12,795 total genes. This table stores CDS/promoter/downstream conservation scores from step 7a.  
**Current impact:** Step 7a data is not used as a direct scoring layer in the current rank-product (step 9 uses convergence + selection + functional evidence). So this gap doesn't affect current scores.  
**Action:** Run step 7a for the missing 7,115 genes. Investigate whether step 7a was incomplete or by design (some genes may lack suitable alignment regions).

### Issue 5 — Three scoring layers not yet implemented (FUTURE)
`gene_therapy_score`, `pathway_convergence`, `regulatory_divergence` tables are all empty. These would add additional layers to the rank product and improve scoring resolution.

---

## 3. Remediation Plan

### Phase A — Targeted PAML Rerun (1,251 genes)
**Scope:** 870 proxy + 381 NULL-model = 1,251 genes  
**PAML-eligible:** 870 proxy (all have motifs) + 33 NULL-model (have motifs) = **903 genes' OGs**  
**Not PAML-eligible:** 348 NULL-model genes have no divergent motifs — they can't enter the PAML scatter pipeline and are correctly excluded (no convergence signal = less scientifically meaningful)

**No code changes needed.** The existing `extract_og_ids` logic already:
- Includes OGs that have divergent motifs and are NOT already scored as `paml_branch_site`
- proxy and NULL genes are not `paml_branch_site` → they are included
- `load_selection_scores` does an upsert by `gene_id` → proxy rows are automatically upgraded

**How to run:**
```bash
# Step 1: Bump paml_og_cache_key in nextflow.config (already done: v1 → v2)
# This forces extract_og_ids to re-query DB with fresh state.

# Step 2: Run step 6 only (PAML)
nextflow run nextflow/main.nf \
  --from_step step6 \
  --until step6 \
  -profile aws,seqera \
  -work-dir s3://bioresilient-data/nf-work

# Step 3: After step 6 completes, rerun steps 8 and 9
nextflow run nextflow/main.nf \
  --from_step step8 --until step9 \
  -profile aws,seqera \
  -work-dir s3://bioresilient-data/nf-work
```

**What happens to proxy genes that still fail CDS coverage check:**
- PAML requires ≥60% CDS coverage from NCBI
- Any proxy gene that fails again → aBSREL proxy stored again with `selection_model='proxy'`
- Their p-value is still treated as 1.0 in scoring (code fix in scoring.py is permanent)
- No data degradation — worst case is no change for those genes

**Reference:** Target gene list saved to `docs/paml_rerun_targets.tsv` (1,251 genes with has_motifs flag)

**Expected gain:** Genes with good CDS data will be upgraded from proxy/NULL → real PAML scores. Conservatively 200–600 additional genes gain selection signal.

---

### Phase B — Rerun 187 Suspicious PAML Entries
**After Phase A completes**, rerun PAML for the 187 genes with capped omega (dnds_ratio=99, pval=1.0) using refined parameters:

```sql
-- Identify them
SELECT gene_id FROM evolution_score 
WHERE selection_model = 'paml_branch_site' 
AND dnds_ratio >= 99 AND dnds_pvalue = 1.0;
```

Then targeted step 6 rerun with improved PAML settings (e.g., checking for enough synonymous sites before running the LRT, or using initial omega = 0.5 instead of 1.0 to avoid degenerate convergence).

**Expected gain:** Some of these could flip to real selection signal. Conservatively 20–50 genes.

---

### Phase C — Complete nucleotide_score coverage (Step 7a gap)
Run step 7a for the 7,115 uncovered genes. This enriches the convergence evidence layer.

```
nextflow run nextflow/main.nf \
  --from_step step7a \
  --until step7a \
  -profile aws,seqera \
  -work-dir s3://bioresilient-data/nf-work
```

---

### Phase D — Implement Missing Scoring Layers (Future)
Once Phases A–C are complete, implement:
1. **gene_therapy_score** — gene delivery vehicle compatibility, viral tropism constraints
2. **pathway_convergence** — enrichment of phenotype pathway in convergent gene set
3. **regulatory_divergence** — cis-regulatory changes at convergent loci

These would extend the rank product from 4 to 6–7 layers, significantly improving resolution among Tier 1 genes.

---

## 4. Immediate vs. Deferred

| Action | Blocking current run? | When |
|---|---|---|
| Steps 8+9 with current fixes | No (running now) | ✅ NOW |
| Phase A: targeted PAML rerun (1,251 genes) | No | After step 9 completes |
| Phase B: suspicious PAML 187 genes | No | After Phase A |
| Phase C: nucleotide_score gap | No | After Phase B |
| Phase D: new scoring layers | No | Future sprint |

---

## 5. What is Clean Right Now

- No HyPhy/proxy p-values influence any score (DB + code both enforced)
- No stale `candidate_score` data (cleared before current run)
- No old DepMap probability scores (switched to Chronos + cleared old rows)
- All 12,795 genes compete in scoring — proxy/NULL genes just have neutral selection signal
- Analytical p-values ensure Tier 1/2/3 are meaningfully distributed
