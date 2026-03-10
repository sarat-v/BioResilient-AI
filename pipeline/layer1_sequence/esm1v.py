"""Step 4c — ESM-1v variant effect scoring for divergent motifs.

ESM-1v (Evolutionary Scale Modeling, variant effect model) computes the
log-likelihood ratio (LLR) of observing a particular amino acid substitution
given the protein's evolutionary context. This is complementary to AlphaMissense:

  - AlphaMissense: predicts pathogenicity (human disease relevance)
  - ESM-1v:        predicts functional impact from evolutionary constraints

For each DivergentMotif:
  1. Reconstruct the full human protein sequence context.
  2. For each divergent position (human_aa ≠ animal_aa), compute:
       LLR = log P(animal_aa | context) - log P(human_aa | context)
  3. Store mean LLR as esm1v_score:
       - Negative LLR: animal change is less likely given protein context
         → adaptive/unusual change, potentially functionally significant
       - Near-zero LLR: biochemically neutral substitution

This gives an ortholog-aware functional prediction beyond simple sequence
comparison, using ESM-1v's 250M-parameter masked language model trained on
250M UniRef90 sequences.

Requirements:
  - `esm` Python package (Meta's ESM library): pip install fair-esm
  - Model is ~10 GB; cached at {storage_root}/esm1v/ on first use
  - Falls back gracefully if ESM not available (non-critical step)

Reference: Meier et al., 2021 "Language models enable zero-shot prediction
of the effects of mutations on protein function" NeurIPS 2021.
"""

import logging
from typing import Optional

from db.models import DivergentMotif, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_storage_root

log = logging.getLogger(__name__)

_ESM1V_MODEL_CACHE = None  # lazy-loaded


def _load_esm1v_model():
    """Load ESM-1v model. Returns (model, alphabet, batch_converter) or None."""
    global _ESM1V_MODEL_CACHE
    if _ESM1V_MODEL_CACHE is not None:
        return _ESM1V_MODEL_CACHE

    try:
        import torch
        import esm as esm_lib

        log.info("Loading ESM-1v model (this may take a few minutes on first use)...")
        model, alphabet = esm_lib.pretrained.esm1v_t33_650M_UR90S_1()
        model.eval()
        batch_converter = alphabet.get_batch_converter()
        _ESM1V_MODEL_CACHE = (model, alphabet, batch_converter)
        log.info("ESM-1v model loaded.")
        return _ESM1V_MODEL_CACHE
    except ImportError:
        log.info("ESM library not installed — skipping ESM-1v scoring (pip install fair-esm to enable).")
        return None
    except Exception as exc:
        log.warning("ESM-1v model load failed: %s", exc)
        return None


def compute_esm1v_llr(
    human_protein_seq: str,
    substitutions: list[tuple[int, str, str]],
    model_tuple,
) -> Optional[float]:
    """Compute mean ESM-1v log-likelihood ratio for a set of substitutions.

    Args:
        human_protein_seq: Full human protein sequence (no gaps).
        substitutions: List of (0-indexed position, human_aa, animal_aa).
        model_tuple: (model, alphabet, batch_converter) from _load_esm1v_model.

    Returns:
        Mean LLR across all substitutions, or None on failure.
    """
    try:
        import torch

        model, alphabet, batch_converter = model_tuple
        device = "cuda" if torch.cuda.is_available() else "cpu"
        model = model.to(device)

        # ESM-1v uses masked marginal scoring:
        # For each position, mask it and compute log P(ref) - log P(alt)
        # LLR(alt) = log P(alt | context) - log P(ref | context)

        llrs = []

        for pos, human_aa, animal_aa in substitutions:
            if pos >= len(human_protein_seq):
                continue

            # Create masked sequence for this position
            masked_seq = human_protein_seq[:pos] + "<mask>" + human_protein_seq[pos + 1:]
            data = [("protein", masked_seq)]
            _, _, batch_tokens = batch_converter(data)
            batch_tokens = batch_tokens.to(device)

            with torch.no_grad():
                logits = model(batch_tokens)["logits"]  # [1, L+2, vocab]

            # Get log-probs at the masked position (offset +1 for <cls> token)
            token_idx = pos + 1
            log_probs = torch.nn.functional.log_softmax(logits[0, token_idx], dim=-1)

            aa_idx_human = alphabet.get_idx(human_aa)
            aa_idx_animal = alphabet.get_idx(animal_aa)

            if aa_idx_human is None or aa_idx_animal is None:
                continue

            llr = (log_probs[aa_idx_animal] - log_probs[aa_idx_human]).item()
            llrs.append(llr)

        return round(sum(llrs) / len(llrs), 4) if llrs else None

    except Exception as exc:
        log.debug("ESM-1v LLR computation failed: %s", exc)
        return None


def annotate_esm1v_scores(gene_ids: Optional[list[str]] = None) -> int:
    """Annotate DivergentMotif rows with ESM-1v variant effect scores.

    For efficiency, loads the ESM-1v model once and processes all motifs.
    Skips gracefully if ESM library is not installed.

    Returns number of motifs scored.
    """
    model_tuple = _load_esm1v_model()
    if model_tuple is None:
        log.info("ESM-1v not available — skipping variant effect scoring.")
        return 0

    with get_session() as session:
        if gene_ids is None:
            gene_ids = [g.id for g in session.query(Gene).all()]
        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()

    log.info("Computing ESM-1v scores for %d genes...", len(genes))
    scored = 0

    for gene in genes:
        # Get the human protein sequence from any human ortholog
        with get_session() as session:
            human_orth = (
                session.query(Ortholog)
                .filter(Ortholog.gene_id == gene.id, Ortholog.species_id == "human")
                .first()
            )
            if not human_orth or not human_orth.protein_seq:
                continue
            human_seq = human_orth.protein_seq.replace("-", "").replace("*", "")

            motifs = (
                session.query(DivergentMotif)
                .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
                .filter(Ortholog.gene_id == gene.id)
                .all()
            )

        for motif in motifs:
            substitutions = []
            for i, (h_aa, a_aa) in enumerate(zip(motif.human_seq, motif.animal_seq)):
                if h_aa == "-" or a_aa == "-" or h_aa == a_aa:
                    continue
                abs_pos = motif.start_pos + i
                if abs_pos < len(human_seq):
                    substitutions.append((abs_pos, h_aa, a_aa))

            if not substitutions:
                continue

            llr = compute_esm1v_llr(human_seq, substitutions, model_tuple)
            if llr is not None:
                with get_session() as session:
                    m = session.get(DivergentMotif, motif.id)
                    if m:
                        m.esm1v_score = llr
                        session.commit()
                scored += 1

    log.info("ESM-1v scoring complete: %d motifs scored.", scored)
    return scored


def run_esm1v_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step4c."""
    return annotate_esm1v_scores(gene_ids)
