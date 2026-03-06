# BioResilient AI — Claude Code Build Package

This folder contains everything Claude Code needs to start building the BioResilient AI platform. Read these files in order.

---

## 📋 Quick Start

**For Claude Code to begin Phase 1 implementation:**

1. **Start here:** `BioResilient_Technical_Plan.md` (this is what Claude Code actually needs to build)
2. **Reference:** `BioResilient_Master_Build_Plan.docx` (scientific context and six-layer architecture explanation)
3. **Reference (if needed):** `BioResilient_ClaudeCode_Plan.md` (existing patterns from earlier version — may have reusable code/schema)

---

## 📄 File Descriptions

### 1. `BioResilient_Technical_Plan.md` ⭐ START HERE FOR CLAUDE CODE

**What it contains:**
- Hardware-agnostic environment setup (local GPU rig OR AWS cloud)
- Exact repository structure
- Conda environment YAML
- SQLAlchemy database schema (all 6 layers)
- Phase 1 implementation order (10 steps, day by day)
- Critical tool flags (OrthoFinder `--diamond`, HyPhy models, etc.)
- Data formats between each step
- Testing strategy
- All database schemas for layers

**What Claude Code should do:**
- Read this first
- Follow the 10-step Phase 1 sequence exactly
- Do not build Phase 2 (Layer 3–6) until Phase 1 is complete
- Reference the database schema before writing any models

**Key constraint:** Phase 1 ONLY. Section 9 lists what NOT to build yet.

---

### 2. `BioResilient_Master_Build_Plan.docx` (PDF VERSION AT `BioResilient_Master_Build_Plan_Scientific_Context.md`)

**What it contains:**
- Strategic thesis (why this approach)
- The six-layer scientific architecture explained plainly
- Cost breakdown (hardware, cloud, personnel)
- Data gaps and honest limitations
- Four-phase 12-month roadmap
- Species selection guidance

**What Claude Code should do:**
- Read for context on WHY each layer exists
- Reference the six-layer descriptions when designing database schema
- Understand the scientific reasoning, not just the code

**This is NOT code instructions.** This is science + strategy.

---

### 3. `BioResilient_ClaudeCode_Plan.md` (REFERENCE ONLY — PARTIALLY OUTDATED)

**What it contains:**
- Repository structure from earlier version (mostly still valid)
- Database schema from POC (needs expansion for 6 layers)
- Pipeline step descriptions
- Some tool/dependency notes

**Status:** This plan predates the six-layer architecture. It covers the earlier POC direction. Use it only to:
- See if any database patterns are reusable
- Understand the original POC structure
- Reference if Claude Code gets stuck

**What Claude Code should do:**
- Use the Technical Plan (file #1) as the source of truth
- Reference this only if you want to see how the POC was originally structured
- Do NOT follow this in preference to the Technical Plan

---

## 🎯 What Claude Code Will Build

**Phase 1 (months 1–3):** A pipeline that runs end-to-end on 15–20 species, producing Layers 1 + 2 scores (sequence divergence, evolutionary selection, convergence detection). No disease annotation yet.

**Hardware:** Agnostic — works on local 4070 Ti rig or AWS `g4dn.xlarge` + `c5.4xlarge` spot instances.

**Entry point:** `BioResilient_Technical_Plan.md` Section 5, Step 1.

---

## ⚙️ Environment Setup

The Technical Plan has two paths:
- **Local:** WSL2 + CUDA 12.x + Conda (for GPU rig with RTX 3080 Ti / 4070 Ti)
- **Cloud:** AWS with spot instances (g4dn.xlarge for GPU, c5.4xlarge for CPU)

All code is identical. Only config file changes.

---

## 🚀 Next Steps for Claude Code

1. Read `BioResilient_Technical_Plan.md` in full (20 mins)
2. Understand the 10 Phase 1 steps (Section 5)
3. Start Step 1 (environment setup)
4. When building database schema (Step 1), reference the models in Section 4 of the Technical Plan
5. Run tests frequently — fixture dataset provided for each layer
6. Ask questions if tool flags or data flow is unclear

---

## 📞 Questions Claude Code Might Have

**"Should I build all six layers now?"**
No. Phase 1 is Layers 1 + 2 only. Layers 3–6 are Phase 2 and beyond.

**"What if local hardware isn't available?"**
Use cloud (AWS). The code is identical. Section 0 of the Technical Plan has cloud setup commands.

**"What about the existing POC code?"**
It exists and may be reusable. Claude Code should audit it (Technical Plan Step 1). Don't assume it's production-ready.

**"Do I need GPU?"**
Yes, for ESMFold and ESM-2 in Phase 2. Phase 1 uses GPU only optionally. It runs on CPU but slower.

---

## 📊 File Map

```
claude-code/
├── INDEX.md                              ← You are here
├── BioResilient_Technical_Plan.md        ← Claude Code reads this first
├── BioResilient_Master_Build_Plan.docx   ← Reference for context
└── BioResilient_ClaudeCode_Plan.md       ← Historical reference only
```

---

**Last updated:** March 1, 2026
**For:** Claude Code to start Phase 1 implementation
**Status:** Ready for build
