"""Generate plain-English research summaries for a gene using Claude.

Reads structured data from the DB (scores, evolution, disease, druggability, motifs)
and calls the Anthropic API to produce a short narrative. Result is cached in Gene.narrative.
"""

import json
import logging
from typing import Any, Optional

from db.models import (
    CandidateScore,
    DiseaseAnnotation,
    DivergentMotif,
    DrugTarget,
    EvolutionScore,
    Gene,
    Ortholog,
)
from db.session import get_session
from pipeline.config import get_anthropic_api_key

log = logging.getLogger(__name__)

SYSTEM_PROMPT = """You are a research assistant for a computational biology platform that finds therapeutic target candidates from naturally disease-resistant species (e.g. naked mole rat, bowhead whale). Given structured data about a human gene and its evolutionary/disease/druggability context, write a concise plain-English summary (2–4 short paragraphs) for a researcher. Cover: (1) why this gene is a candidate (convergence across species, selection signal), (2) disease relevance and existing drug links if any, (3) one or two open research questions. Be precise and avoid hype. Do not invent data not present in the input."""


def _gather_gene_context(gene_id: str, session) -> dict[str, Any]:
    """Load all structured data for a gene for the LLM prompt."""
    gene = session.get(Gene, gene_id)
    if not gene:
        return {}
    cs = session.query(CandidateScore).filter_by(gene_id=gene_id, trait_id="").first()
    ev = session.get(EvolutionScore, gene_id)
    da = session.get(DiseaseAnnotation, gene_id)
    dt = session.get(DrugTarget, gene_id)
    orthologs = session.query(Ortholog).filter_by(gene_id=gene_id).all()
    motifs = []
    for o in orthologs:
        for m in o.motifs:
            motifs.append({
                "species_id": o.species_id,
                "start": m.start_pos,
                "end": m.end_pos,
                "human_seq": m.human_seq,
                "animal_seq": m.animal_seq,
                "divergence_score": m.divergence_score,
            })
    return {
        "gene_id": gene.id,
        "gene_symbol": gene.gene_symbol,
        "human_protein": gene.human_protein,
        "human_gene_id": gene.human_gene_id,
        "composite_score": cs.composite_score if cs else None,
        "tier": cs.tier if cs else None,
        "convergence_score": cs.convergence_score if cs else None,
        "selection_score": cs.selection_score if cs else None,
        "expression_score": cs.expression_score if cs else None,
        "dnds_ratio": ev.dnds_ratio if ev else None,
        "dnds_pvalue": ev.dnds_pvalue if ev else None,
        "selection_model": ev.selection_model if ev else None,
        "convergence_count": ev.convergence_count if ev else None,
        "branches_under_selection": ev.branches_under_selection if ev else None,
        "phylop_score": ev.phylop_score if ev else None,
        "disease_name": da.disease_name if da else None,
        "opentargets_score": da.opentargets_score if da else None,
        "gwas_pvalue": da.gwas_pvalue if da else None,
        "gnomad_pli": da.gnomad_pli if da else None,
        "pocket_count": dt.pocket_count if dt else None,
        "top_pocket_score": dt.top_pocket_score if dt else None,
        "existing_drugs": dt.existing_drugs if dt else None,
        "druggability_tier": dt.druggability_tier if dt else None,
        "ortholog_species": [o.species_id for o in orthologs],
        "motifs_sample": motifs[:5],
    }


def generate_narrative(gene_id: str, force: bool = False) -> tuple[str, bool]:
    """Generate or return cached narrative for a gene.

    Returns (narrative_text, cached). cached is True if we returned existing Gene.narrative.
    Raises ValueError if gene not found or API key missing.
    """
    api_key = get_anthropic_api_key()
    if not api_key:
        raise ValueError("ANTHROPIC_API_KEY not set; add to config or environment")

    with get_session() as session:
        gene = session.get(Gene, gene_id)
        if not gene:
            raise ValueError(f"Gene '{gene_id}' not found")
        if not force and getattr(gene, "narrative", None):
            return (gene.narrative or "").strip(), True

        context = _gather_gene_context(gene_id, session)
        if not context:
            raise ValueError(f"No context for gene '{gene_id}'")

    user_content = (
        "Summarize this gene as a therapeutic target candidate for a researcher. "
        "Use only the data below.\n\n"
        f"Structured data (JSON):\n{json.dumps(context, indent=2)}"
    )

    try:
        import anthropic
        client = anthropic.Anthropic(api_key=api_key)
        msg = client.messages.create(
            model="claude-3-5-sonnet-20241022",
            max_tokens=1024,
            system=SYSTEM_PROMPT,
            messages=[{"role": "user", "content": user_content}],
        )
        text = ""
        for block in msg.content:
            if hasattr(block, "text"):
                text += block.text
        narrative = text.strip()
    except Exception as e:
        log.exception("Anthropic API error for gene %s: %s", gene_id, e)
        raise RuntimeError(f"LLM request failed: {e}") from e

    with get_session() as session:
        gene = session.get(Gene, gene_id)
        if gene:
            gene.narrative = narrative
            session.commit()

    return narrative, False
