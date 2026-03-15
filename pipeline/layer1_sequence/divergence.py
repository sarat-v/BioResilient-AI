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
from pathlib import Path
from typing import Optional

import numpy as np

from sqlalchemy import text

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
    """Insert DivergentMotif rows into the database using bulk inserts.

    Args:
        motifs: Output from extract_divergent_motifs().

    Returns:
        Number of rows inserted.
    """
    import time
    import uuid as _uuid

    if not motifs:
        return 0

    thresholds = get_thresholds()
    max_identity = thresholds.get("divergence_identity_max", 0.85)

    # Build lookup by fetching all orthologs in a single query — faster than chunked IN()
    unique_pairs = {(m["species_id"], m["protein_id"]) for m in motifs}
    log.info("  Loading %d motifs (%d unique orthologs) to DB...", len(motifs), len(unique_pairs))
    log.info("  Building ortholog lookup (single full-table scan)...")
    t0 = time.time()

    ortholog_lookup: dict[tuple, str] = {}  # (species_id, protein_id) → ortholog.id
    identity_lookup: dict[str, float] = {}  # ortholog.id → sequence_identity_pct

    with get_session() as session:
        all_orthologs = session.query(
            Ortholog.id, Ortholog.species_id, Ortholog.protein_id, Ortholog.sequence_identity_pct
        ).yield_per(10000)
        for r in all_orthologs:
            if (r.species_id, r.protein_id) in unique_pairs:
                ortholog_lookup[(r.species_id, r.protein_id)] = r.id
                identity_lookup[r.id] = r.sequence_identity_pct
        log.info("  Ortholog lookup built: %d entries in %.1fs.", len(ortholog_lookup), time.time() - t0)

    # Pre-build the full insert list in memory (filter first, then insert in committed chunks)
    log.info("  Pre-filtering motifs against identity threshold (max_identity=%.2f)...", max_identity)
    t1 = time.time()
    rows_to_insert = []
    skipped_no_ortholog = 0
    skipped_identity = 0
    for m in motifs:
        key = (m["species_id"], m["protein_id"])
        oid = ortholog_lookup.get(key)
        if oid is None:
            skipped_no_ortholog += 1
            continue
        pct = identity_lookup.get(oid)
        if pct is not None and pct / 100.0 > max_identity:
            skipped_identity += 1
            continue
        rows_to_insert.append({
            "id": str(_uuid.uuid4()),
            "ortholog_id": oid,
            "start_pos": m["start_pos"],
            "end_pos": m["end_pos"],
            "animal_seq": m["animal_seq"],
            "human_seq": m["human_seq"],
            "divergence_score": m["divergence_score"],
            "esm_distance": m.get("esm_distance"),
        })
    log.info(
        "  Pre-filter done in %.1fs: %d to insert, %d skipped (no ortholog), %d skipped (identity).",
        time.time() - t1, len(rows_to_insert), skipped_no_ortholog, skipped_identity,
    )

    # Use PostgreSQL COPY for maximum throughput — ~10-50x faster than bulk_insert_mappings
    # COPY streams all data over a single connection with no per-row overhead
    import io
    from db.session import get_engine

    total_to_insert = len(rows_to_insert)
    log.info("  Inserting %d motif rows via COPY (streaming)...", total_to_insert)
    t2 = time.time()

    # Write rows as TSV into an in-memory buffer, flushing every 500k rows to cap memory
    inserted = 0
    chunk_size = 500_000
    raw_conn = get_engine().raw_connection()
    try:
        cur = raw_conn.cursor()
        for chunk_start in range(0, total_to_insert, chunk_size):
            chunk = rows_to_insert[chunk_start: chunk_start + chunk_size]
            buf = io.StringIO()
            for r in chunk:
                esm = r["esm_distance"] if r["esm_distance"] is not None else "\\N"
                buf.write(
                    f"{r['id']}\t{r['ortholog_id']}\t{r['start_pos']}\t{r['end_pos']}\t"
                    f"{r['animal_seq']}\t{r['human_seq']}\t{r['divergence_score']}\t{esm}\n"
                )
            buf.seek(0)
            cur.copy_from(
                buf, "divergent_motif",
                columns=("id", "ortholog_id", "start_pos", "end_pos",
                         "animal_seq", "human_seq", "divergence_score", "esm_distance"),
            )
            raw_conn.commit()
            inserted += len(chunk)
            elapsed = time.time() - t2
            rate = inserted / elapsed if elapsed > 0 else 0
            log.info(
                "  Motif insert: %d / %d rows copied (%.0f rows/s, ETA ~%.0f min).",
                inserted, total_to_insert, rate,
                (total_to_insert - inserted) / rate / 60 if rate > 0 else 0,
            )
        cur.close()
    finally:
        raw_conn.close()

    log.info("  Inserted %d divergent motif rows in %.1fs.", inserted, time.time() - t2)
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


def _divergence_worker_simple(aligned_seqs: dict[str, str]) -> list[dict]:
    """Worker for pool.map — takes sequences only, returns motifs list."""
    try:
        return extract_divergent_motifs("", aligned_seqs)
    except Exception:
        return []


def run_divergence_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
) -> dict[str, list[dict]]:
    """Run sliding window divergence for all aligned orthogroups.

    Pure Python sliding window — fast enough sequentially (no subprocess overhead).
    Returns {og_id: [motif_dict, ...]}
    """
    items = list(aligned_orthogroups.items())
    total = len(items)

    all_motifs: dict[str, list[dict]] = {}
    motif_total = 0

    log.info("Divergence scan: %d orthogroups (sequential sliding window).", total)

    for done, (og_id, aligned_seqs) in enumerate(items, 1):
        motifs = extract_divergent_motifs(og_id, aligned_seqs)
        if motifs:
            all_motifs[og_id] = motifs
            motif_total += len(motifs)
        if done % 2000 == 0 or done == total:
            log.info("  Divergence: %d / %d orthogroups (%.0f%%).", done, total, 100 * done / total)

    log.info("Divergence scan: %d motifs across %d orthogroups.", motif_total, len(all_motifs))
    return all_motifs


def update_sequence_identities(aligned_orthogroups: dict[str, dict[str, str]]) -> None:
    """Update Ortholog.sequence_identity_pct in the database from alignment data.

    Uses a bulk UPDATE via a single pass: compute all identities in memory,
    then batch-update in chunks to avoid N+1 queries.
    """
    from pipeline.layer1_sequence.alignment import calculate_sequence_identity

    # Build {(species_id, protein_id): identity} map in memory first
    log.info("  Computing sequence identities for %d orthogroups...", len(aligned_orthogroups))
    identity_map: dict[tuple[str, str], float] = {}
    for og_id, aligned_seqs in aligned_orthogroups.items():
        human_seq = _get_human_sequence(aligned_seqs)
        if human_seq is None:
            continue
        for label, animal_seq in aligned_seqs.items():
            if "human" in label:
                continue
            species_id, protein_id = _parse_label(label)
            identity_map[(species_id, protein_id)] = round(
                calculate_sequence_identity(human_seq, animal_seq), 2
            )
    log.info("  Identity map built: %d pairs.", len(identity_map))

    if not identity_map:
        return

    import io
    import time
    from db.session import get_engine

    t0 = time.time()
    log.info("  Bulk-updating %d sequence identity rows via COPY+temp table...", len(identity_map))

    items = list(identity_map.items())
    raw_conn = get_engine().raw_connection()
    try:
        cur = raw_conn.cursor()
        # 1. Create temp table
        cur.execute("""
            CREATE TEMP TABLE _identity_update (
                species_id TEXT,
                protein_id TEXT,
                pct        REAL
            )
        """)
        # 2. Stream all rows via COPY (fastest possible load — no per-row overhead)
        buf = io.StringIO()
        for (sid, pid), pct in items:
            buf.write(f"{sid}\t{pid}\t{pct}\n")
        buf.seek(0)
        cur.copy_from(buf, "_identity_update", columns=("species_id", "protein_id", "pct"))
        # 3. Single UPDATE join
        cur.execute("""
            UPDATE ortholog AS o
            SET sequence_identity_pct = u.pct
            FROM _identity_update AS u
            WHERE o.species_id = u.species_id
              AND o.protein_id = u.protein_id
        """)
        updated = cur.rowcount
        cur.execute("DROP TABLE _identity_update")
        raw_conn.commit()
        cur.close()
    finally:
        raw_conn.close()

    log.info("Sequence identities updated: %d rows in %.1fs.", updated, time.time() - t0)
