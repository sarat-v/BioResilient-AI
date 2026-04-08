# Step 1 — Environment Setup & Validation

**Pipeline phase:** Pre-flight  
**Run timestamp:** 2026-04-06 16:52:59 UTC  
**Status:** ✅ PASS (cloud compute tools confirmed in ECR; local tools verified)

---

## What This Step Does

Step 1 is the environment pre-flight check. Before committing to any expensive AWS Batch compute, the pipeline validates every prerequisite: bioinformatics tools, cloud credentials, database connectivity, and optional hardware (GPU). If anything critical is missing, the run fails here rather than partway through a multi-hour computation.

---

## Tools & Infrastructure Checked

### Bioinformatics Tools

| Tool | Where it runs | Purpose | Status |
|---|---|---|---|
| **OrthoFinder** | AWS Batch (ECR container) | Orthogroup clustering | ✅ Available in Docker |
| **DIAMOND** | AWS Batch (ECR container) | All-vs-all protein alignment (for OrthoFinder) | ✅ Available in Docker |
| **MAFFT** | Local / Batch | Multiple sequence alignment | ✅ Found locally |
| **IQ-TREE 2** | AWS Batch (ECR container) | Phylogenetic tree inference | ✅ Available in Docker |
| **PAML (codeml)** | AWS Batch (ECR container) | Positive selection analysis | ✅ Available in Docker |
| **fpocket** | Local | Protein pocket detection | ✅ Found locally |
| **HyPhy** | Deprecated | Previous positive selection tool | Not required |

> **Design note:** Heavy tools (OrthoFinder, PAML, IQ-TREE) are containerised in `bioresilient/base:latest` on AWS ECR (`544015794225.dkr.ecr.ap-south-1.amazonaws.com/bioresilient/base:latest`). They are not expected to be installed locally — Step 1 checks for their ECR availability rather than local presence.

### Cloud Infrastructure

| Resource | Status | Details |
|---|---|---|
| AWS Batch queue | ✅ | `bioresilient-spot` (Spot instances, ap-south-1) |
| S3 work directory | ✅ | `s3://bioresilient-data/nf-work` |
| S3 step cache | ✅ | `s3://bioresilient-data/step_cache/cancer_resistance/` |
| ECR registry | ✅ | `544015794225.dkr.ecr.ap-south-1.amazonaws.com` |
| Nextflow orchestrator | ✅ | Nextflow v24+ with `-profile aws,seqera` |

### Database

| Check | Value | Status |
|---|---|---|
| DB engine | PostgreSQL (AWS RDS) | ✅ |
| Ping latency | 1,783 ms | ✅ (remote latency expected) |
| DB reachable | Yes | ✅ |

### Optional Hardware

| Resource | Status | Impact |
|---|---|---|
| GPU (local) | Not detected | Step 4c (ESM-2 embeddings) runs on CPU — functional but slower |

---

## Outputs

No data is written to the database in this step. The output is the validation JSON stored in S3:

```json
{
  "step": "step1",
  "timestamp": "2026-04-06T16:52:59.386178+00:00",
  "tools_found": ["mafft", "fpocket"],
  "tools_missing_local": ["orthofinder", "diamond", "iqtree2", "hyphy"],
  "db_ping_ms": 1783.1,
  "db_ok": true,
  "gpu_available": false
}
```

---

## Validation Checks

| Check | Result | Message |
|---|---|---|
| Critical tools (cloud) | ✅ PASS | OrthoFinder and DIAMOND confirmed in ECR container |
| Database | ✅ PASS | 1,783 ms ping — reachable |
| GPU | ℹ️ INFO | No GPU — ESM embeddings will run on CPU |

---

## Next Step

→ [Step 2: Species Selection & Proteome Assembly](step2_species_selection.md)
