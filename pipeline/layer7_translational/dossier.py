"""Step 18 — Target dossier generation.

For each Tier1/Tier2 gene, generates a structured Markdown target dossier
combining all pipeline evidence into a single researcher-facing document.

Output: outputs/dossiers/{gene_symbol}.md
"""

import logging
from pathlib import Path
from typing import Optional

from db.models import (
    CandidateScore, DiseaseAnnotation, DrugTarget,
    EvolutionScore, Gene, GeneTherapyScore,
    PreclinicalReadiness, SafetyFlag,
)
from db.session import get_session

log = logging.getLogger(__name__)

_SPECIES_SUFFIXES = ("_HUMAN", "_MOUSE", "_RAT", "_DOG", "_CAT", "_WHALE", "_BAT", "_SHARK")

OUTPUT_DIR = Path("outputs/dossiers")


def _hgnc(symbol: str) -> str:
    for suf in _SPECIES_SUFFIXES:
        if symbol.upper().endswith(suf):
            return symbol[: -len(suf)]
    return symbol


def _fmt(val, fmt=".3f", default="—") -> str:
    if val is None:
        return default
    try:
        return format(float(val), fmt)
    except (TypeError, ValueError):
        return str(val)


_PHENOTYPE_ASSAY_CONTEXT: dict[str, dict[str, str]] = {
    "cancer_resistance": {
        "knockdown": "CRISPR/RNAi knockdown in cancer cell lines to confirm growth-inhibition phenotype",
        "repurposing": "Repurposing screen: test {drug} in cancer-resistance model",
        "in_vivo": "AAV-mediated gene overexpression in mouse tumour model to assess cancer resistance",
        "organoid": "Patient-derived tumour organoid (PDO) screen with candidate inhibitors",
    },
    "longevity": {
        "knockdown": "CRISPR/RNAi knockdown in C. elegans or Drosophila to assess lifespan impact",
        "repurposing": "Senescence-reversal assay: test {drug} in aged primary fibroblasts",
        "in_vivo": "AAV-mediated overexpression in aged mice; track metabolic health and survival curves",
        "organoid": "Aged organoid system to test candidate modulators of cellular senescence",
    },
    "viral_immunity": {
        "knockdown": "CRISPR/RNAi knockdown in primary macrophages followed by viral challenge assay",
        "repurposing": "Innate immunity reporter assay: test {drug} in VSV/SARS-CoV-2 infection model",
        "in_vivo": "AAV-mediated overexpression in mice; challenge with relevant virus and track viral load",
        "organoid": "Lung organoid viral infection model with candidate innate-immune modulators",
    },
    "viral_resistance": {
        "knockdown": "CRISPR/RNAi knockdown in primary macrophages followed by viral challenge assay",
        "repurposing": "Innate immunity reporter assay: test {drug} in VSV/SARS-CoV-2 infection model",
        "in_vivo": "AAV-mediated overexpression in mice; challenge with relevant virus and track viral load",
        "organoid": "Lung organoid viral infection model with candidate innate-immune modulators",
    },
    "hypoxia_tolerance": {
        "knockdown": "CRISPR/RNAi knockdown in cardiomyocytes under hypoxic (1% O₂) conditions",
        "repurposing": "Hypoxia-survival assay: test {drug} in HIF-reporter cell line",
        "in_vivo": "AAV-mediated overexpression in mouse chronic hypoxia model (normobaric 10% O₂)",
        "organoid": "Cardiac organoid hypoxia-reoxygenation assay with candidate modulators",
    },
    "dna_repair": {
        "knockdown": "CRISPR/RNAi knockdown followed by ionising radiation survival clonogenic assay",
        "repurposing": "Synthetic lethality screen: test {drug} against BRCA1/2-null background",
        "in_vivo": "AAV-mediated overexpression in mouse DNA-damage (irradiation) model",
        "organoid": "DNA-damage organoid model (IR or cisplatin) to test candidate repair modulators",
    },
}

_DEFAULT_ASSAY_CONTEXT = {
    "knockdown": "CRISPR/RNAi knockdown in relevant primary cells to confirm phenotype",
    "repurposing": "Cell-based assay: test {drug} in appropriate disease model",
    "in_vivo": "AAV-mediated overexpression in mouse model to assess target phenotype",
    "organoid": "Organoid screen with candidate modulators in disease-relevant model",
}


def _experiment_recommendations(
    tier: str,
    therapy: Optional[GeneTherapyScore],
    drug_target: Optional[DrugTarget],
    preclinical: Optional[PreclinicalReadiness],
    phenotype: str = "",
) -> list[str]:
    """Return 3 recommended first experiments based on available evidence and phenotype."""
    ctx = _PHENOTYPE_ASSAY_CONTEXT.get(phenotype, _DEFAULT_ASSAY_CONTEXT)
    recs = []

    # Knock-down / overexpression validation
    recs.append(ctx["knockdown"])

    # Biochemical / compound assay
    if drug_target and drug_target.top_pocket_score and drug_target.top_pocket_score > 0.5:
        recs.append("Fragment-based or virtual screening against the identified druggable pocket")
    elif drug_target and drug_target.existing_drugs:
        recs.append(ctx["repurposing"].format(drug=drug_target.existing_drugs[0]))
    else:
        recs.append("Protein purification + differential scanning fluorimetry (DSF) for binding assay")

    # In-vivo / organoid
    if therapy and therapy.aav_compatible:
        recs.append(ctx["in_vivo"])
    else:
        recs.append(ctx["organoid"])

    return recs[:3]


def generate_dossier(gene_id: str, phenotype: str = "") -> Optional[str]:
    """Generate a Markdown dossier for one gene. Returns the output path or None."""
    with get_session() as session:
        gene = session.get(Gene, gene_id)
        if not gene:
            return None

        ev  = session.get(EvolutionScore, gene_id)
        ann = session.get(DiseaseAnnotation, gene_id)
        dt  = session.get(DrugTarget, gene_id)
        sf  = session.get(SafetyFlag, gene_id)
        gt  = session.get(GeneTherapyScore, gene_id)
        pr  = session.get(PreclinicalReadiness, gene_id)

        # Get best CandidateScore (highest composite)
        cs_rows = (
            session.query(CandidateScore)
            .filter(CandidateScore.gene_id == gene_id)
            .order_by(CandidateScore.composite_score.desc())
            .all()
        )
        cs = cs_rows[0] if cs_rows else None

    symbol = _hgnc(gene.gene_symbol or gene_id)
    tier   = cs.tier if cs else "Unscored"
    score  = cs.composite_score if cs else 0.0

    lines = [
        f"# Target Dossier: {symbol}",
        "",
        f"**Gene ID:** {gene_id}  ",
        f"**Tier:** {tier}  ",
        f"**Composite Score:** {_fmt(score)}",
        "",
        "---",
        "",
        "## 1. Evolutionary Evidence",
        "",
    ]

    if ev:
        lines += [
            f"- **PAML p-value:** {_fmt(ev.dnds_pvalue, '.4e')}  (selection model: {ev.selection_model or '—'})",
            f"- **dN/dS ratio:** {_fmt(ev.dnds_ratio)}",
            f"- **Convergence weight:** {_fmt(ev.convergence_weight)}",
            f"- **Convergence p-value:** {_fmt(ev.convergence_pval, '.4e')}",
            f"- **Convergent lineage count:** {ev.convergence_count or 0}",
            "",
        ]
    else:
        lines += ["*No evolutionary scores available.*", ""]

    lines += ["## 2. Clinical & Disease Evidence", ""]

    if ann:
        lines += [
            f"- **Top disease:** {ann.disease_name or '—'}  (OpenTargets score: {_fmt(ann.ot_association_score)})",
            f"- **GWAS p-value:** {_fmt(ann.gwas_pvalue, '.2e')}",
            f"- **gnomAD pLI:** {_fmt(ann.gnomad_pli)}",
            f"- **Known drug phase:** {ann.known_drug_phase or 'None'}",
            f"- **Literature score:** {_fmt(ann.lit_score)} ({ann.lit_pmid_count or 0} papers)",
            "",
        ]
    else:
        lines += ["*No disease annotation available.*", ""]

    lines += ["## 3. Druggability", ""]

    if dt:
        drugs = ", ".join(dt.existing_drugs[:5]) if dt.existing_drugs else "None identified"
        lines += [
            f"- **Druggability tier:** {dt.druggability_tier or '—'}",
            f"- **Top pocket score:** {_fmt(dt.top_pocket_score)} ({dt.pocket_count or 0} pockets)",
            f"- **Tractability:** SM={dt.tractability_sm}, AB={dt.tractability_ab}, PROTAC={dt.tractability_protac}",
            f"- **Known drugs:** {drugs}",
            f"- **P2Rank score:** {_fmt(dt.p2rank_score)}",
            "",
        ]
    else:
        lines += ["*No druggability data available.*", ""]

    lines += ["## 4. Safety Profile", ""]

    if sf:
        lines += [
            f"- **Broadly essential (DepMap):** {_fmt(sf.depmap_score)} {'⚠ HIGH RISK' if sf.depmap_score and sf.depmap_score > 0.7 else ''}",
            f"- **Hub risk (STRING):** {sf.hub_risk}",
            f"- **Essential (pLI):** {sf.is_essential}",
            f"- **Family size:** {sf.family_size or '—'}",
            "",
        ]
    else:
        lines += ["*No safety data available.*", ""]

    lines += ["## 5. Gene Therapy Potential", ""]

    if gt:
        lines += [
            f"- **AAV-compatible:** {gt.aav_compatible}  (gene size: {gt.gene_size_bp or '—'} bp)",
            f"- **CRISPR sites:** {gt.crispr_sites or '—'}  (risk: {gt.offtarget_risk or '—'})",
            f"- **Guide efficiency:** {_fmt(gt.crispr_efficiency)}",
            f"- **Tissue tropism:** {', '.join(gt.tissue_tropism or []) or '—'}",
            "",
        ]
    else:
        lines += ["*No gene therapy data available.*", ""]

    lines += ["## 6. Preclinical Readiness", ""]

    if pr:
        lines += [
            f"- **Overall score:** {_fmt(pr.overall_readiness_score)} — **{pr.readiness_tier}**",
            f"- **Synthesis score:** {_fmt(pr.synthesis_score)}",
            f"- **Deliverability score:** {_fmt(pr.deliverability_score)}",
            f"- **Model availability:** {_fmt(pr.model_score)}",
            f"- **Assay score:** {_fmt(pr.assay_score)}",
        ]
        if pr.notes:
            lines.append(f"- **Notes:** {pr.notes}")
        lines.append("")
    else:
        lines += ["*Preclinical readiness not yet scored.*", ""]

    lines += ["## 7. Recommended First Experiments", ""]

    recs = _experiment_recommendations(tier, gt, dt, pr, phenotype=phenotype)
    for i, rec in enumerate(recs, 1):
        lines.append(f"{i}. {rec}")
    lines.append("")
    phenotype_label = phenotype or "unspecified"
    lines += [
        "---",
        "",
        f"*Generated by BioResilient Pipeline | Gene: {gene_id} | Phenotype: {phenotype_label}*",
    ]

    content = "\n".join(lines)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_DIR / f"{symbol}.md"
    out_path.write_text(content, encoding="utf-8")
    log.info("Dossier written: %s", out_path)
    return str(out_path)


def run_dossier_pipeline(gene_ids: list[str], phenotype: str = "") -> int:
    """Entry point called from orchestrator step18. Returns number of dossiers written."""
    written = 0
    for gid in gene_ids:
        try:
            path = generate_dossier(gid, phenotype=phenotype)
            if path:
                written += 1
        except Exception as exc:
            log.warning("Dossier generation failed for %s: %s", gid, exc)
    log.info("Dossiers: %d written to %s/ (phenotype=%r)", written, OUTPUT_DIR, phenotype or "unspecified")
    return written
