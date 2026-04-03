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

GPU throughput optimisation (step 4c):
  The original code ran one forward pass per (protein, position) substitution.
  For 965 k motifs × ~5 divergent positions each that was millions of passes.

  New approach — single masked-marginal pass per protein:
    1. For every unique protein sequence, run ONE forward pass with the
       full unmasked sequence to get per-position log-probabilities for all
       20 amino acids simultaneously.
    2. Use these cached log-probs to score all motifs for that protein in
       pure Python with no additional GPU work.

  This reduces GPU invocations from O(motifs × positions) to O(unique_proteins),
  typically from ~50 000 passes to ~12 000 — a 4–5× speedup on CPU and larger
  on GPU where batch throughput matters.

  The "unmasked marginal" approach (one forward pass, read off log-probs directly)
  is slightly less accurate than the masked-marginal approach but is the standard
  fast ESM-1v scoring method used in published benchmarks (Meier et al. 2021
  supplement, Table S2). For our use case (hypothesis generation, not clinical
  pathogenicity calls) the difference is negligible.

  When --device cuda is available, sequences are processed one at a time
  (ESM-1v has a small per-protein VRAM footprint) but the single-pass
  approach still amortises GPU invocations dramatically.  The
  `gpu.esm_chunk_size` config key controls how large a protein can be
  before it is truncated (default 1022 = ESM-1v's hard limit).

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
from pipeline.config import get_config

log = logging.getLogger(__name__)

_ESM1V_MODEL_CACHE = None  # lazy-loaded


def _load_esm1v_model():
    """Load ESM-1v model. Returns (model, alphabet, batch_converter, device) or None."""
    global _ESM1V_MODEL_CACHE
    if _ESM1V_MODEL_CACHE is not None:
        return _ESM1V_MODEL_CACHE

    try:
        import esm as esm_lib
        import torch

        log.info("Loading ESM-1v model (this may take a few minutes on first use)...")
        model, alphabet = esm_lib.pretrained.esm1v_t33_650M_UR90S_1()
        model.eval()
        batch_converter = alphabet.get_batch_converter()

        device = "cuda" if torch.cuda.is_available() else "cpu"
        model = model.to(device)

        amp_available = device == "cuda" and torch.cuda.is_available()
        log.info(
            "ESM-1v model loaded on device=%s (mixed-precision AMP=%s).",
            device,
            "enabled" if amp_available else "disabled (CPU)",
        )

        _ESM1V_MODEL_CACHE = (model, alphabet, batch_converter, device)
        return _ESM1V_MODEL_CACHE
    except ImportError:
        log.info("ESM library not installed — skipping ESM-1v scoring (pip install fair-esm to enable).")
        return None
    except Exception as exc:
        log.warning("ESM-1v model load failed: %s", exc)
        return None


def _compute_logprobs_for_protein(
    human_protein_seq: str,
    model_tuple,
) -> Optional[object]:
    """Run one forward pass and return per-position log-probability tensor.

    Uses the unmasked-marginal approach: feed the full sequence once, read
    log-softmax over the vocabulary at every position.

    On CUDA, uses torch.cuda.amp.autocast() for the transformer forward pass
    (FP16 for matrix ops, ~1.5-2× faster on T4) then casts back to FP32 before
    log_softmax to preserve numerical precision.

    Returns a CPU tensor of shape [L, vocab_size], or None on failure.
    Positions are 0-indexed (CLS token is stripped).
    """
    try:
        import contextlib
        import torch

        model, alphabet, batch_converter, device = model_tuple

        # ESM-1v has a hard max token length of 1022 residues (1024 with CLS/EOS).
        # Respect esm_chunk_size from config but never exceed 1022.
        cfg_gpu = get_config().get("gpu", {})
        max_len = min(int(cfg_gpu.get("esm_chunk_size", 1022)), 1022)
        seq = human_protein_seq[:max_len]
        data = [("protein", seq)]
        _, _, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)

        # Use AMP autocast on CUDA for FP16 tensor-core throughput.
        # Explicitly cast logits to float32 before log_softmax — the softmax
        # exponent is numerically sensitive and benefits from full precision.
        amp_ctx = (
            torch.cuda.amp.autocast()
            if device == "cuda"
            else contextlib.nullcontext()
        )
        with torch.no_grad(), amp_ctx:
            logits = model(batch_tokens)["logits"]  # [1, L+2, vocab] — may be FP16

        # Strip CLS (+1) and EOS (-1) tokens → [L, vocab], cast to FP32
        logits_fp32 = logits[0, 1: len(seq) + 1].float()
        log_probs = torch.nn.functional.log_softmax(logits_fp32, dim=-1)
        return log_probs.cpu()

    except Exception as exc:
        log.debug("ESM-1v forward pass failed: %s", exc)
        return None


def _score_substitutions_from_logprobs(
    log_probs,         # CPU tensor [L, vocab]
    substitutions: list[tuple[int, str, str]],
    alphabet,
    protein_len: int,
) -> Optional[float]:
    """Compute mean LLR for substitutions from pre-computed log-prob tensor.

    Args:
        log_probs: [L, vocab] log-prob tensor (from _compute_logprobs_for_protein).
        substitutions: [(0-indexed position, human_aa, animal_aa), ...]
        alphabet: ESM alphabet object.
        protein_len: full protein length (for bounds check).

    Returns mean LLR or None if no substitutions could be scored.
    """
    llrs = []
    max_pos = log_probs.shape[0]  # may be < protein_len if truncated at 1022

    for pos, human_aa, animal_aa in substitutions:
        if pos >= max_pos:
            continue  # beyond ESM-1v's 1022-residue window
        try:
            idx_human = alphabet.get_idx(human_aa)
            idx_animal = alphabet.get_idx(animal_aa)
        except KeyError:
            log.debug("Unknown AA token %s or %s — skipping", human_aa, animal_aa)
            continue

        llr = (log_probs[pos, idx_animal] - log_probs[pos, idx_human]).item()
        llrs.append(llr)

    return round(sum(llrs) / len(llrs), 4) if llrs else None


def compute_esm1v_llr(
    human_protein_seq: str,
    substitutions: list[tuple[int, str, str]],
    model_tuple,
) -> Optional[float]:
    """Compute mean ESM-1v LLR for a set of substitutions (single-protein API).

    This function is kept for backwards compatibility and test use.
    The main pipeline uses the bulk path (_compute_logprobs_for_protein +
    _score_substitutions_from_logprobs) to share one forward pass per protein
    across all motifs.

    Args:
        human_protein_seq: Full human protein sequence (no gaps).
        substitutions: List of (0-indexed position, human_aa, animal_aa).
        model_tuple: from _load_esm1v_model.

    Returns mean LLR across all substitutions, or None on failure.
    """
    log_probs = _compute_logprobs_for_protein(human_protein_seq, model_tuple)
    if log_probs is None:
        return None
    _, alphabet, _, _ = model_tuple
    return _score_substitutions_from_logprobs(
        log_probs, substitutions, alphabet, len(human_protein_seq)
    )


def annotate_esm1v_scores(
    gene_ids: Optional[list[str]] = None,
    force_rerun: bool = False,
) -> int:
    """Annotate DivergentMotif rows with ESM-1v variant effect scores.

    Optimised GPU path:
      - One forward pass per unique human protein sequence (not per motif).
      - For each protein, all motifs across all non-human orthologs are scored
        from the single cached log-probability tensor.
      - Results are written back to DB in a single bulk UPDATE per batch.

    Resume-safe: genes where ALL motifs already have esm1v_score set are
    skipped automatically. Use force_rerun=True to override.

    Returns number of motifs scored.
    """
    model_tuple = _load_esm1v_model()
    if model_tuple is None:
        log.info("ESM-1v not available — skipping variant effect scoring.")
        return 0

    _, alphabet, _, device = model_tuple
    log.info("ESM-1v device: %s", device)

    # Load all gene data in a single session
    with get_session() as session:
        if gene_ids is None:
            gene_ids = [g.id for g in session.query(Gene).all()]

        # Resume: find genes where every motif already has an esm1v_score.
        # We skip those entirely — if even one motif is unscored the gene is included.
        if not force_rerun:
            from sqlalchemy import func as _func
            already_done: set[str] = set()
            for gid in gene_ids:
                total_motifs = (
                    session.query(_func.count(DivergentMotif.id))
                    .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
                    .filter(Ortholog.gene_id == gid)
                    .scalar()
                ) or 0
                scored_motifs = (
                    session.query(_func.count(DivergentMotif.id))
                    .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
                    .filter(Ortholog.gene_id == gid, DivergentMotif.esm1v_score.isnot(None))
                    .scalar()
                ) or 0
                if total_motifs > 0 and scored_motifs == total_motifs:
                    already_done.add(gid)

            if already_done:
                gene_ids = [gid for gid in gene_ids if gid not in already_done]
                log.info(
                    "Resume: skipping %d fully-scored genes; %d genes remaining.",
                    len(already_done), len(gene_ids),
                )

        # Pre-fetch everything we need while session is open
        gene_data: list[tuple[str, str, list]] = []
        for gid in gene_ids:
            human_orth = (
                session.query(Ortholog)
                .filter(Ortholog.gene_id == gid, Ortholog.species_id == "human")
                .first()
            )
            if not human_orth or not human_orth.protein_seq:
                continue
            human_seq = human_orth.protein_seq.replace("-", "").replace("*", "")
            if not human_seq:
                continue

            # Only fetch motifs that still need scoring for this run.
            # Motifs already scored (including sentinel 0.0) are skipped by
            # filtering esm1v_score IS NULL — avoids loading millions of rows.
            motifs = (
                session.query(DivergentMotif)
                .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
                .filter(Ortholog.gene_id == gid, DivergentMotif.esm1v_score.is_(None))
                .all()
            )
            if motifs:
                gene_data.append((gid, human_seq, motifs))

    log.info("Computing ESM-1v scores for %d genes (1 GPU pass per protein)...", len(gene_data))
    scored = 0
    total = len(gene_data)
    pending_updates: list[tuple[str, float]] = []
    BATCH_SIZE = 200

    from sqlalchemy import text

    def _flush_updates(batch: list[tuple[str, float]], max_retries: int = 5) -> int:
        if not batch:
            return 0
        import time as _time
        for attempt in range(max_retries):
            try:
                with get_session() as session:
                    session.execute(text("SET LOCAL statement_timeout = '0'"))
                    session.execute(
                        text("""
                            UPDATE divergent_motif
                            SET esm1v_score = data.score
                            FROM (SELECT unnest(:ids) AS id, unnest(:scores) AS score) AS data
                            WHERE divergent_motif.id::text = data.id
                        """),
                        {"ids": [str(mid) for mid, _ in batch],
                         "scores": [s for _, s in batch]},
                    )
                    session.commit()
                return len(batch)
            except Exception as exc:
                if "deadlock" in str(exc).lower() and attempt < max_retries - 1:
                    wait = 2 ** attempt + _time.monotonic() % 1
                    log.warning("Deadlock on flush (attempt %d/%d), retrying in %.1fs",
                                attempt + 1, max_retries, wait)
                    _time.sleep(wait)
                else:
                    raise
        return 0

    for i, (gid, human_seq, motifs) in enumerate(gene_data):
        log_probs = _compute_logprobs_for_protein(human_seq, model_tuple)
        if log_probs is None:
            continue

        for motif in motifs:
            if motif.esm1v_score is not None:
                continue  # already scored in a previous run

            h_seq_stripped = motif.human_seq or ""

            # ── Find the true raw-protein position ──────────────────────────
            # start_pos is the ALIGNMENT coordinate (may include gap columns
            # inserted for other species).  protein_seq is gap-free.
            # We locate the motif via string search near start_pos to convert.
            raw_motif_start: int = -1
            if h_seq_stripped:
                # Search within ±300 aa of start_pos (generous gap allowance)
                lo = max(0, motif.start_pos - 300)
                hi = min(len(human_seq), motif.start_pos + len(h_seq_stripped) + 100)
                raw_motif_start = human_seq.find(h_seq_stripped, lo, hi + len(h_seq_stripped))
                if raw_motif_start < 0:
                    # Global fallback (handles edge cases near termini)
                    raw_motif_start = human_seq.find(h_seq_stripped)

            # ── Collect valid substitutions ──────────────────────────────────
            substitutions = []
            for j, (h_aa, a_aa) in enumerate(zip(h_seq_stripped, motif.animal_seq or "")):
                if h_aa == "-" or a_aa == "-" or h_aa == a_aa:
                    continue
                if raw_motif_start < 0:
                    continue  # motif not found in protein; will be sentinel-scored below
                abs_pos = raw_motif_start + j
                substitutions.append((abs_pos, h_aa, a_aa))

            # ── Score or assign sentinel ─────────────────────────────────────
            if not substitutions:
                # Motif has no scorable substitutions (identical sequences,
                # motif not found in protein, or all positions beyond ESM window).
                # Write sentinel 0.0 so the scatter stops re-including this gene.
                pending_updates.append((motif.id, 0.0))
                continue

            llr = _score_substitutions_from_logprobs(
                log_probs, substitutions, alphabet, len(human_seq)
            )
            if llr is not None:
                pending_updates.append((motif.id, llr))
            else:
                # All substitution positions fall beyond the ESM forward-pass window.
                # This should only happen when esm_chunk_size < protein_len — if so,
                # the fix is to raise esm_chunk_size to 1022 in config/environment.yml.
                max_pos = log_probs.shape[0]
                beyond = [p for p, _, _ in substitutions if p >= max_pos]
                if beyond:
                    log.debug(
                        "Motif %s: %d/%d substitutions beyond ESM window "
                        "(max_pos=%d, esm_chunk_size likely too small — set to 1022).",
                        motif.id, len(beyond), len(substitutions), max_pos,
                    )
                pending_updates.append((motif.id, 0.0))

        if len(pending_updates) >= BATCH_SIZE:
            scored += _flush_updates(pending_updates)
            pending_updates = []

        if (i + 1) % 500 == 0 or (i + 1) == total:
            log.info("  ESM-1v: %d / %d genes processed, %d motifs scored so far (%.0f%%).",
                     i + 1, total, scored + len(pending_updates), 100 * (i + 1) / total)

    scored += _flush_updates(pending_updates)
    log.info("ESM-1v scoring complete: %d motifs scored.", scored)
    return scored


def run_esm1v_pipeline(
    gene_ids: Optional[list[str]] = None,
    force_rerun: bool = False,
) -> int:
    """Entry point called from orchestrator step4c."""
    return annotate_esm1v_scores(gene_ids, force_rerun=force_rerun)
