"""Step 4b — Sliding-window divergence scoring.

For each aligned orthogroup:
  1. Slide a window (size=15, step=5) across the alignment.
  2. Score each window by fraction of differing positions (divergence %).
  3. Optionally compute ESM-2 embedding cosine distance for high-divergence windows.
  4. Load DivergentMotif records into the database.

Threshold: a motif is flagged if sequence identity < 85% in ≥ 2 protective species.
"""

import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import numpy as np

from db.models import DivergentMotif, Ortholog
from db.session import get_session
from pipeline.config import get_thresholds

log = logging.getLogger(__name__)

WINDOW_SIZE = 15
WINDOW_STEP = 5


def score_window(animal_window: str, human_window: str) -> float:
    """Fraction of positions that differ (ignoring gap-only columns).

    Returns a divergence score in [0.0, 1.0].
    """
    if len(animal_window) != len(human_window):
        return 0.0

    diffs = 0
    positions = 0
    for a, h in zip(animal_window, human_window):
        if a == "-" and h == "-":
            continue
        positions += 1
        if a != h:
            diffs += 1

    return diffs / positions if positions > 0 else 0.0


def extract_divergent_motifs(
    og_id: str,
    aligned_seqs: dict[str, str],
    min_divergence: float = 0.15,
) -> list[dict]:
    """Slide a window over each animal sequence vs. the human sequence.

    Args:
        og_id: OrthoGroup identifier.
        aligned_seqs: {label: aligned_sequence} — output from MAFFT.
        min_divergence: Minimum divergence score to flag a window.

    Returns:
        List of motif dicts with keys:
          species_id, protein_id, start_pos, end_pos,
          animal_seq, human_seq, divergence_score.
    """
    human_seq = _get_human_sequence(aligned_seqs)
    if human_seq is None:
        return []

    alignment_len = len(human_seq)
    motifs = []

    for label, animal_seq in aligned_seqs.items():
        if "human" in label:
            continue
        if len(animal_seq) != alignment_len:
            continue

        species_id, protein_id = _parse_label(label)

        for start in range(0, alignment_len - WINDOW_SIZE + 1, WINDOW_STEP):
            end = start + WINDOW_SIZE
            animal_window = animal_seq[start:end]
            human_window = human_seq[start:end]

            # Skip windows that are mostly gaps
            if animal_window.count("-") > WINDOW_SIZE * 0.5:
                continue

            div_score = score_window(animal_window, human_window)
            if div_score < min_divergence:
                continue

            motifs.append({
                "og_id": og_id,
                "species_id": species_id,
                "protein_id": protein_id,
                "start_pos": start,
                "end_pos": end,
                "animal_seq": animal_window.replace("-", ""),
                "human_seq": human_window.replace("-", ""),
                "divergence_score": round(div_score, 4),
            })

    return motifs


def _get_human_sequence(aligned_seqs: dict[str, str]) -> Optional[str]:
    for label, seq in aligned_seqs.items():
        if "human" in label.lower():
            return seq
    return None


def _parse_label(label: str) -> tuple[str, str]:
    """Parse '{species_id}|{protein_id}' label."""
    parts = label.split("|", 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return label, label


def compute_esm_distance(seq_a: str, seq_b: str, model=None) -> Optional[float]:
    """Compute cosine distance between ESM-2 embeddings of two peptide sequences.

    Returns None if ESM-2 is not available (Phase 2 full usage).
    In Phase 1 this is called on small motif windows only.
    """
    try:
        import torch
        import esm

        if model is None:
            model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
            model.eval()
            batch_converter = alphabet.get_batch_converter()
        else:
            model, alphabet, batch_converter = model

        device = "cuda" if torch.cuda.is_available() else "cpu"
        model = model.to(device)

        data = [("seq_a", seq_a), ("seq_b", seq_b)]
        _, _, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33])
            embeddings = results["representations"][33]

        emb_a = embeddings[0, 1:-1].mean(0)
        emb_b = embeddings[1, 1:-1].mean(0)

        cos_sim = torch.nn.functional.cosine_similarity(
            emb_a.unsqueeze(0), emb_b.unsqueeze(0)
        ).item()
        return round(1.0 - cos_sim, 6)

    except ImportError:
        log.debug("ESM-2 not available — skipping embedding distance.")
        return None
    except Exception as exc:
        log.warning("ESM-2 embedding failed: %s", exc)
        return None


def load_motifs_to_db(
    motifs: list[dict],
    ortholog_map: dict = None,  # unused, kept for backwards compatibility
) -> int:
    """Insert DivergentMotif rows into the database.

    Args:
        motifs: Output from extract_divergent_motifs().

    Returns:
        Number of rows inserted.
    """
    thresholds = get_thresholds()
    max_identity = thresholds.get("divergence_identity_max", 0.85)

    inserted = 0
    with get_session() as session:
        for m in motifs:
            # Look up ortholog by species_id + protein_id
            ortholog = (
                session.query(Ortholog)
                .filter_by(species_id=m["species_id"], protein_id=m["protein_id"])
                .first()
            )
            if ortholog is None:
                continue

            # Only keep motifs from species with low overall identity (divergent enough)
            if ortholog.sequence_identity_pct and ortholog.sequence_identity_pct / 100.0 > max_identity:
                continue

            motif = DivergentMotif(
                ortholog_id=ortholog.id,
                start_pos=m["start_pos"],
                end_pos=m["end_pos"],
                animal_seq=m["animal_seq"],
                human_seq=m["human_seq"],
                divergence_score=m["divergence_score"],
                esm_distance=m.get("esm_distance"),
            )
            session.add(motif)
            inserted += 1

    return inserted


def filter_by_independent_lineages(
    motifs_by_og: dict[str, list[dict]],
    species_to_lineage: dict[str, str],
    min_lineages: int = 2,
) -> dict[str, list[dict]]:
    """Gate 2: keep only orthogroups with divergent motifs in ≥ min_lineages independent lineage groups.

    species_to_lineage: {species_id: lineage_group} e.g. {'naked_mole_rat': 'Rodents', ...}
    Expected reduction: ~2k–4k → 500–800 orthogroups passed to HyPhy.
    """
    filtered: dict[str, list[dict]] = {}
    for og_id, motifs in motifs_by_og.items():
        lineages_with_motif = set()
        for m in motifs:
            sid = m.get("species_id")
            if sid and sid in species_to_lineage:
                lineages_with_motif.add(species_to_lineage[sid])
        if len(lineages_with_motif) >= min_lineages:
            filtered[og_id] = motifs
    log.info(
        "  Gate 2: %d orthogroups pass lineage recurrence filter (≥%d independent lineages).",
        len(filtered),
        min_lineages,
    )
    return filtered


def _divergence_worker(args: tuple) -> tuple[str, list[dict]]:
    """Top-level worker for multiprocessing (must be picklable)."""
    og_id, aligned_seqs = args
    return og_id, extract_divergent_motifs(og_id, aligned_seqs)


def run_divergence_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
) -> dict[str, list[dict]]:
    """Run sliding window divergence for all aligned orthogroups in parallel.

    Returns {og_id: [motif_dict, ...]}
    """
    n_cpu = os.cpu_count() or 8
    # Use fewer workers to avoid memory pressure with 12k large orthogroups in flight
    n_workers = max(1, min(n_cpu - 2, 32))
    items = list(aligned_orthogroups.items())
    total = len(items)

    all_motifs: dict[str, list[dict]] = {}
    motif_total = 0
    done = 0

    log.info("Divergence scan: %d orthogroups with %d parallel workers.", total, n_workers)

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        # Submit in batches of 500 to avoid holding 12k futures in memory
        batch_size = 500
        for batch_start in range(0, total, batch_size):
            batch = items[batch_start: batch_start + batch_size]
            futures = {pool.submit(_divergence_worker, item): item[0] for item in batch}
            for future in as_completed(futures, timeout=300):
                try:
                    og_id, motifs = future.result(timeout=60)
                    if motifs:
                        all_motifs[og_id] = motifs
                        motif_total += len(motifs)
                except Exception as exc:
                    og_id = futures[future]
                    log.debug("Divergence worker failed for %s: %s", og_id, exc)
                done += 1
                if done % 1000 == 0 or done == total:
                    log.info("  Divergence: %d / %d orthogroups (%.0f%%).", done, total, 100 * done / total)

    log.info("Divergence scan: %d motifs across %d orthogroups.", motif_total, len(all_motifs))
    return all_motifs


def update_sequence_identities(aligned_orthogroups: dict[str, dict[str, str]]) -> None:
    """Update Ortholog.sequence_identity_pct in the database from alignment data."""
    from pipeline.layer1_sequence.alignment import calculate_sequence_identity

    with get_session() as session:
        for og_id, aligned_seqs in aligned_orthogroups.items():
            human_seq = _get_human_sequence(aligned_seqs)
            if human_seq is None:
                continue

            for label, animal_seq in aligned_seqs.items():
                if "human" in label:
                    continue
                species_id, protein_id = _parse_label(label)

                identity = calculate_sequence_identity(human_seq, animal_seq)

                ortholog = (
                    session.query(Ortholog)
                    .filter_by(species_id=species_id, protein_id=protein_id)
                    .first()
                )
                if ortholog:
                    ortholog.sequence_identity_pct = round(identity, 2)

    log.info("Sequence identities updated.")
