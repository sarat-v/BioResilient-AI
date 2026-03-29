"""Step 4b — Sliding-window divergence scoring.

For each aligned orthogroup:
  1. Slide a window (size=15, step=5) across the alignment.
  2. Score each window by fraction of differing positions (divergence %).
  3. Optionally compute ESM-2 embedding cosine distance for high-divergence windows.
  4. Load DivergentMotif records into the database.

Threshold: a motif is flagged if sequence identity < 85% in ≥ 2 protective species.
All thresholds are read from config/environment.yml at runtime — no hardcoded values.
"""

import logging
import os
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Optional

import numpy as np

from sqlalchemy import text

from db.models import DivergentMotif, Ortholog
from db.session import get_session
from pipeline.config import get_thresholds

log = logging.getLogger(__name__)


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
    min_divergence: float = None,
) -> list[dict]:
    """Slide a window over each animal sequence vs. the human sequence.

    Args:
        og_id: OrthoGroup identifier.
        aligned_seqs: {label: aligned_sequence} — output from MAFFT.
        min_divergence: Minimum divergence score to flag a window. Reads from
            config (divergence_min_score) when None.

    Returns:
        List of motif dicts with keys:
          species_id, protein_id, start_pos, end_pos,
          animal_seq, human_seq, divergence_score.
    """
    thresholds = get_thresholds()
    if min_divergence is None:
        min_divergence = float(thresholds.get("divergence_min_score", 0.15))
    window_size = int(thresholds.get("divergence_window_size", 15))
    window_step = int(thresholds.get("divergence_window_step", 5))
    min_species_per_window = int(thresholds.get("divergence_min_species_per_window", 2))

    human_seq = _get_human_sequence(aligned_seqs)
    if human_seq is None:
        return []

    alignment_len = len(human_seq)
    motifs = []

    # Pre-compute which non-human sequences are valid (correct length)
    non_human_seqs = {
        label: seq for label, seq in aligned_seqs.items()
        if "human" not in label and len(seq) == alignment_len
    }

    for label, animal_seq in non_human_seqs.items():
        species_id, protein_id = _parse_label(label)

        for start in range(0, alignment_len - window_size + 1, window_step):
            end = start + window_size
            animal_window = animal_seq[start:end]
            human_window = human_seq[start:end]

            # Skip windows that are mostly gaps in this species
            if animal_window.count("-") > window_size * 0.5:
                continue

            # Species coverage guard: require min_species_per_window to have
            # valid (non-gappy) sequence at this position to avoid alignment
            # artifact signals where only 1 species has residues.
            valid_species_at_window = sum(
                1 for seq in non_human_seqs.values()
                if seq[start:end].count("-") <= window_size * 0.5
            )
            if valid_species_at_window < min_species_per_window:
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
    # Keep only the top-N highest-divergence motifs per ortholog to cap DB volume.
    # Value is read from config (divergence_top_n_per_ortholog, default 8).
    thresholds = get_thresholds()
    TOP_N_PER_ORTHOLOG = int(thresholds.get("divergence_top_n_per_ortholog", 8))
    log.info(
        "  Pre-filtering motifs (max_identity=%.2f, top-%d per ortholog)...",
        max_identity, TOP_N_PER_ORTHOLOG,
    )
    t1 = time.time()

    # Group by ortholog, keep top-N by divergence score
    motifs_by_ortholog: dict[str, list] = defaultdict(list)
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
        motifs_by_ortholog[oid].append(m)

    rows_to_insert = []
    for oid, omotifs in motifs_by_ortholog.items():
        top = sorted(omotifs, key=lambda x: x["divergence_score"], reverse=True)[:TOP_N_PER_ORTHOLOG]
        for m in top:
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
    import io
    from pipeline.config import get_psycopg2_conn

    # Delete existing motifs for the orthologs we're about to re-insert.
    # This makes step4 idempotent — reruns replace data instead of appending duplicates.
    affected_ortholog_ids = list({r["ortholog_id"] for r in rows_to_insert})
    if affected_ortholog_ids:
        log.info("  Deleting existing motifs for %d orthologs (idempotent rerun safety)...", len(affected_ortholog_ids))
        t_del = time.time()
        conn = get_psycopg2_conn()
        try:
            cur = conn.cursor()
            batch_sz = 5000
            total_deleted = 0
            for i in range(0, len(affected_ortholog_ids), batch_sz):
                batch = affected_ortholog_ids[i : i + batch_sz]
                cur.execute("DELETE FROM divergent_motif WHERE ortholog_id = ANY(%s)", (batch,))
                total_deleted += cur.rowcount
            conn.commit()
            cur.close()
        finally:
            conn.close()
        log.info("  Deleted %d stale motif rows in %.1fs.", total_deleted, time.time() - t_del)

    total_to_insert = len(rows_to_insert)
    log.info("  Inserting %d motif rows via COPY (streaming)...", total_to_insert)
    t2 = time.time()

    inserted = 0
    chunk_size = 100_000
    conn = get_psycopg2_conn()
    try:
        cur = conn.cursor()
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
            conn.commit()
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
        conn.close()

    log.info("  Inserted %d divergent motif rows in %.1fs.", inserted, time.time() - t2)
    return inserted


def filter_by_independent_lineages(
    motifs_by_og: dict[str, list[dict]],
    species_to_lineage: dict[str, str],
    min_lineages: int = 2,
) -> dict[str, list[dict]]:
    """Gate 2: keep only orthogroups with divergent motifs in ≥ min_lineages independent lineage groups.

    species_to_lineage: {species_id: lineage_group} — must be pre-filtered by the caller
        to exclude control species (is_control=True) and baseline species.
        Control species' motifs in passing OGs are preserved (they provide background-
        branch sequences for HyPhy) but do NOT contribute to the lineage count.

    Example caller pattern (see orchestrator.step4):
        species_to_lineage = {
            s["id"]: s["lineage_group"]
            for s in species_list
            if not s.get("is_control") and "baseline" not in s.get("phenotypes", [])
        }

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


def branch_length_weights_from_tree(
    treefile: Path,
    reference_species: str = "human",
) -> dict[str, float]:
    """Compute per-species divergence weights from an IQ-TREE Newick species tree.

    Weight = min_distance / distance_to_reference, so the closest species to
    *reference_species* gets weight 1.0 and more distant species get lower weights.
    Closely related species' divergence signals are therefore up-weighted relative
    to background noise from very distant lineages.

    Returns {} (no-op) if the treefile is missing, ete3 is unavailable, or
    the reference species is not a leaf in the tree.
    """
    if not treefile.exists():
        log.info("Branch-length weighting: treefile not found at %s (expected on first run).", treefile)
        return {}

    try:
        from ete3 import Tree  # type: ignore
    except ImportError:
        log.warning("ete3 not installed — branch-length divergence weighting skipped.")
        return {}

    try:
        tree = Tree(str(treefile))
    except Exception as exc:
        log.warning("Could not parse treefile %s: %s — weighting skipped.", treefile, exc)
        return {}

    leaf_names = {n.name for n in tree.get_leaves()}
    if reference_species not in leaf_names:
        log.warning(
            "'%s' not found in tree leaves %s — branch-length weighting skipped.",
            reference_species, sorted(leaf_names),
        )
        return {}

    distances: dict[str, float] = {}
    for name in leaf_names:
        if name == reference_species:
            continue
        try:
            d = tree.get_distance(reference_species, name)
            if d and d > 0:
                distances[name] = d
        except Exception as exc:
            log.debug("Could not compute distance to %s: %s", name, exc)

    if not distances:
        log.warning("No valid distances computed from species tree — weighting skipped.")
        return {}

    min_dist = min(distances.values())
    weights = {sid: round(min_dist / d, 4) for sid, d in distances.items()}

    log.info(
        "Branch-length weights from %s (reference='%s', min_dist=%.4f):",
        treefile.name, reference_species, min_dist,
    )
    for sid, w in sorted(weights.items(), key=lambda x: -x[1]):
        log.info("  %-30s  dist=%.4f  weight=%.4f", sid, distances[sid], w)

    return weights


def apply_branch_length_weighting(
    motifs_by_og: dict[str, list[dict]],
    species_weights: dict[str, float],
) -> dict[str, list[dict]]:
    """Scale each motif's divergence_score by the species' branch-length weight.

    Motifs from evolutionary distant species (e.g., hydra, molluscs) get lower
    divergence scores, reducing background noise in downstream convergence tests.
    Species absent from species_weights are left unchanged (implicit weight = 1.0).
    """
    if not species_weights:
        return motifs_by_og

    updated: dict[str, list[dict]] = {}
    for og_id, motifs in motifs_by_og.items():
        new_motifs = []
        for m in motifs:
            w = species_weights.get(m.get("species_id", ""), 1.0)
            if w != 1.0:
                m = dict(m)
                m["divergence_score"] = round(float(m["divergence_score"]) * w, 4)
            new_motifs.append(m)
        updated[og_id] = new_motifs

    n_weighted = sum(1 for sid in species_weights)
    log.info(
        "Branch-length weighting applied to %d species across %d orthogroups.",
        n_weighted, len(updated),
    )
    return updated


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

    # Step4 QC report — helps detect pipeline drift across runs
    if all_motifs:
        counts = [len(v) for v in all_motifs.values()]
        log.info(
            "Step4 QC: motifs_raw=%d, orthogroups_with_motifs=%d, "
            "median_per_og=%.1f, max_per_og=%d, min_per_og=%d",
            motif_total, len(all_motifs),
            statistics.median(counts), max(counts), min(counts),
        )

    return all_motifs


def update_sequence_identities(aligned_orthogroups: dict[str, dict[str, str]]) -> None:
    """Update Ortholog.sequence_identity_pct in the database from alignment data.

    Strategy: fetch all ortholog (id, species_id, protein_id) in one SELECT,
    resolve identity values in Python, then COPY (id, pct) into a temp table
    and UPDATE by primary key — avoids full table scan on (species_id, protein_id).
    """
    from pipeline.layer1_sequence.alignment import calculate_sequence_identity

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
    from pipeline.config import get_psycopg2_conn

    t0 = time.time()

    # Resolve (species_id, protein_id) → ortholog.id in a single SELECT
    # Then COPY (id, pct) and UPDATE by primary key — uses PK index, no table scan
    conn = get_psycopg2_conn()
    try:
        cur = conn.cursor()

        log.info("  Fetching ortholog id map from DB (single SELECT)...")
        cur.execute("SELECT id, species_id, protein_id FROM ortholog")
        id_rows = cur.fetchall()  # ~184k rows, fast over RDS
        log.info("  Fetched %d ortholog rows in %.1fs.", len(id_rows), time.time() - t0)

        # Resolve to (ortholog_id, pct) pairs in Python
        updates: list[tuple[str, float]] = []
        for oid, sid, pid in id_rows:
            pct = identity_map.get((sid, pid))
            if pct is not None:
                updates.append((oid, pct))
        log.info("  Resolved %d id→pct pairs. Bulk-updating via COPY+PK join...", len(updates))

        cur.execute("""
            CREATE TEMP TABLE _identity_update (
                ortholog_id TEXT,
                pct         REAL
            )
        """)
        buf = io.StringIO()
        for oid, pct in updates:
            buf.write(f"{oid}\t{pct}\n")
        buf.seek(0)
        cur.copy_from(buf, "_identity_update", columns=("ortholog_id", "pct"))

        # UPDATE by primary key — instant with PK index
        cur.execute("""
            UPDATE ortholog AS o
            SET sequence_identity_pct = u.pct
            FROM _identity_update AS u
            WHERE o.id = u.ortholog_id
        """)
        updated = cur.rowcount
        conn.commit()
        cur.close()
    finally:
        conn.close()

    log.info("  Sequence identities updated: %d rows in %.1fs.", updated, time.time() - t0)
