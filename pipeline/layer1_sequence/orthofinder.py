"""Step 3 — OrthoFinder wrapper.

Runs OrthoFinder on all downloaded proteomes simultaneously using the
DIAMOND backend (-S diamond), which is 500x faster than BLAST.

After the run, parses OrthoGroups.tsv and loads into the Ortholog table,
keeping only groups where the human protein is present.
"""

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO

import json
import time

import requests

from db.models import Gene, Ortholog
from db.session import get_session
from pipeline.config import get_local_storage_root, get_tool_config

log = logging.getLogger(__name__)


_HGNC_CACHE_PATH = Path(get_local_storage_root()) / "hgnc_symbol_cache.json" if False else None
_HGNC_CACHE: dict[str, str] = {}


def _load_hgnc_cache() -> None:
    global _HGNC_CACHE_PATH
    _HGNC_CACHE_PATH = Path(get_local_storage_root()) / "hgnc_symbol_cache.json"
    if _HGNC_CACHE_PATH.exists():
        try:
            _HGNC_CACHE.update(json.loads(_HGNC_CACHE_PATH.read_text()))
        except Exception:
            pass


def _save_hgnc_cache() -> None:
    if _HGNC_CACHE_PATH is None:
        return
    try:
        _HGNC_CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
        _HGNC_CACHE_PATH.write_text(json.dumps(_HGNC_CACHE, indent=2))
    except Exception:
        pass


def _resolve_hgnc_symbol(uniprot_accession: str) -> str:
    """Resolve a UniProt accession (e.g. P04637) to an HGNC gene symbol (e.g. TP53).

    Falls back to the raw accession if the API is unavailable or returns no result.
    Results are cached on disk between pipeline runs.
    """
    if not _HGNC_CACHE:
        _load_hgnc_cache()

    key = uniprot_accession.split("|")[-1].strip()
    bare = key.split("_")[0] if "_" in key else key   # AKT1_HUMAN → AKT1

    if bare in _HGNC_CACHE:
        return _HGNC_CACHE[bare]

    # Skip API call for already-readable symbols (e.g. AKT1, TP53)
    if bare.isalpha() and bare.isupper() and 2 <= len(bare) <= 12:
        _HGNC_CACHE[bare] = bare
        return bare

    # Query UniProt REST API for the gene symbol
    try:
        r = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query": f"accession:{bare} AND reviewed:true",
                "fields": "gene_names",
                "format": "json",
                "size": 1,
            },
            timeout=8,
        )
        if r.status_code == 200:
            results = r.json().get("results", [])
            if results:
                names = results[0].get("genes", [])
                if names:
                    symbol = names[0].get("geneName", {}).get("value", bare)
                    _HGNC_CACHE[bare] = symbol
                    _save_hgnc_cache()
                    return symbol
        time.sleep(0.05)
    except Exception as exc:
        log.debug("HGNC symbol lookup %s: %s", bare, exc)

    _HGNC_CACHE[bare] = bare
    return bare


def _orthofinder_output_dir(proteomes_dir: Path) -> Optional[Path]:
    """OrthoFinder creates a dated output folder — locate the most recently modified one."""
    import glob
    results_dir = proteomes_dir.parent / "orthofinder_out" / "Results_*"
    matches = glob.glob(str(results_dir))
    if not matches:
        return None
    # Sort by modification time (most recent last) rather than lexicographic order
    # to correctly handle dates like Results_Jan01 vs Results_Dec31.
    return Path(max(matches, key=lambda p: Path(p).stat().st_mtime))


def _prepare_orthofinder_input(proteomes_dir: Path) -> Path:
    """Build a staging dir with only valid .reheadered.faa (one per species).

    OrthoFinder must see exactly one FASTA per species. The download step leaves
    both raw .faa and .reheadered.faa in proteomes_dir; failed species (e.g. 0 proteins)
    can leave empty/invalid .faa. Using only .reheadered.faa avoids duplicates and
    invalid inputs that break DIAMOND makedb.
    """
    root = Path(get_local_storage_root())
    staging = root / "orthofinder_input"
    staging.mkdir(parents=True, exist_ok=True)
    for p in staging.iterdir():
        p.unlink()

    for src in sorted(proteomes_dir.glob("*.reheadered.faa")):
        if src.stat().st_size < 100:
            log.warning("  Skipping empty or tiny %s", src.name)
            continue
        # First line must be FASTA (">") for DIAMOND
        with open(src) as f:
            first = f.readline().strip()
        if not first.startswith(">"):
            log.warning("  Skipping non-FASTA %s", src.name)
            continue
        species_id = src.stem.replace(".reheadered", "")
        dest = staging / f"{species_id}.fa"
        if dest.exists():
            dest.unlink()
        dest.symlink_to(src.resolve())
    n = len(list(staging.iterdir()))
    log.info("  OrthoFinder input: %d species in %s", n, staging)
    return staging


def run_orthofinder(proteomes_dir: Path) -> Path:
    """Run OrthoFinder with DIAMOND on the proteomes directory.

    If a completed results directory already exists, reuse it (cache).
    Builds a staging dir with only .reheadered.faa (one per species) so OrthoFinder
    does not see raw .faa or invalid files.

    Args:
        proteomes_dir: Directory containing .faa and .reheadered.faa per species.

    Returns:
        Path to the OrthoFinder results directory.
    """
    cfg = get_tool_config()
    n_cpu = os.cpu_count() or 8
    search_threads = cfg.get("orthofinder_threads", n_cpu)
    align_threads  = cfg.get("orthofinder_align_threads", max(1, n_cpu // 4))

    root = Path(get_local_storage_root())
    out_dir = root / "orthofinder_out"

    # Reuse cached results if they exist
    existing = _orthofinder_output_dir(proteomes_dir)
    if existing is not None and (existing / "Orthogroups" / "Orthogroups.tsv").exists():
        log.info("OrthoFinder results already exist — reusing cache: %s", existing)
        return existing

    # Staging dir: only valid reheadered FASTA (one per species)
    input_dir = _prepare_orthofinder_input(proteomes_dir)
    if len(list(input_dir.iterdir())) == 0:
        raise RuntimeError("No valid proteome files for OrthoFinder (check *.reheadered.faa in proteomes dir)")

    # Remove stale/incomplete output dir so OrthoFinder can create a fresh one
    if out_dir.exists():
        log.info("Removing stale OrthoFinder output dir: %s", out_dir)
        shutil.rmtree(out_dir)
    # Do NOT mkdir here — OrthoFinder creates the output dir itself

    cmd = [
        "orthofinder",
        "-f", str(input_dir),
        "-S", "diamond",
        "-M", "dendroblast",
        "-t", str(search_threads),
        "-a", str(align_threads),
        "-o", str(out_dir),
        "-y",
    ]

    log.info("Running OrthoFinder: %s", " ".join(cmd))
    log_file = out_dir.parent / "orthofinder_run.log"
    log.info("OrthoFinder stdout/stderr → %s", log_file)
    with open(log_file, "w") as lf:
        result = subprocess.run(
            cmd,
            stdout=lf,
            stderr=subprocess.STDOUT,
            check=False,
        )
    # Echo the tail of the log to the pipeline logger so errors are visible
    try:
        tail = log_file.read_text().splitlines()[-30:]
        for line in tail:
            log.debug("[orthofinder] %s", line)
    except Exception:
        pass
    if result.returncode != 0:
        raise RuntimeError(
            f"OrthoFinder exited with code {result.returncode} — see {log_file}"
        )

    results_path = _orthofinder_output_dir(proteomes_dir)
    if results_path is None:
        raise FileNotFoundError("OrthoFinder output directory not found after run.")
    if not (results_path / "Orthogroups" / "Orthogroups.tsv").exists():
        raise FileNotFoundError(
            "Orthogroups.tsv missing after OrthoFinder run (DIAMOND or OrthoFinder may have failed; check log above)"
        )

    log.info("OrthoFinder complete → %s", results_path)
    return results_path


def parse_orthogroups(results_dir: Path, proteomes_dir: Path) -> pd.DataFrame:
    """Parse OrthoGroups.tsv into a long-format DataFrame.

    Returns columns: [og_id, species_id, protein_id]
    Only keeps orthogroups where the 'human' column is non-empty.
    """
    og_path = results_dir / "Orthogroups" / "Orthogroups.tsv"
    if not og_path.exists():
        raise FileNotFoundError(f"Orthogroups.tsv not found at {og_path}")

    df = pd.read_csv(og_path, sep="\t", index_col=0)
    df.index.name = "og_id"
    df = df.reset_index()

    # When we feed only .reheadered.faa (as orthofinder_input/species.fa), columns
    # are one per species (e.g. human, naked_mole_rat). When both raw and reheadered
    # exist, keep only reheadered columns.
    reheadered_cols = [c for c in df.columns if "reheadered" in c]
    if reheadered_cols:
        df = df[["og_id"] + reheadered_cols]
    else:
        species_cols = [c for c in df.columns if c != "og_id"]
        df = df[["og_id"] + species_cols]

    # Keep only groups where human has at least one protein
    human_col = _find_species_column(df, "human")
    if human_col is None:
        raise ValueError("Could not find 'human' column in OrthoGroups.tsv")

    df = df[df[human_col].notna() & (df[human_col] != "")]
    log.info("  %d orthogroups with human proteins", len(df))

    # Melt to long format
    species_cols = [c for c in df.columns if c != "og_id"]
    long_rows = []
    for _, row in df.iterrows():
        og_id = row["og_id"]
        for col in species_cols:
            cell = str(row[col]) if pd.notna(row[col]) else ""
            if not cell:
                continue
            for protein_id in cell.split(", "):
                protein_id = protein_id.strip()
                if protein_id:
                    species_id = _col_to_species_id(col, proteomes_dir)
                    long_rows.append({
                        "og_id": og_id,
                        "species_id": species_id,
                        "protein_id": protein_id,
                    })

    return pd.DataFrame(long_rows)


def _find_species_column(df: pd.DataFrame, species_id: str) -> Optional[str]:
    """Match a species_id to a DataFrame column (flexible matching on stem).
    
    Prefer the reheadered column since that's what we loaded into seq_map.
    """
    # Prefer exact reheadered match first
    for col in df.columns:
        if f"{species_id}.reheadered" in col.lower():
            return col
    # Fall back to any match
    for col in df.columns:
        if species_id in col.lower():
            return col
    return None


def _col_to_species_id(col: str, proteomes_dir: Path) -> str:
    """Map an OrthoFinder column name (FASTA stem) to a species_id.

    OrthoFinder uses the filename stem of each input FASTA as the column name.
    The reheadered FASTAs are named {species_id}.reheadered.faa, so the stem
    is {species_id}.reheadered — strip the extra suffix.
    """
    stem = col.replace(".reheadered", "").replace(".faa", "")
    return stem


def load_sequence_map(proteomes_dir: Path) -> dict[str, dict[str, str]]:
    """Build map: {species_id: {protein_id: sequence}} from reheadered FASTA files."""
    seq_map: dict[str, dict[str, str]] = {}
    for faa in proteomes_dir.glob("*.reheadered.faa"):
        species_id = faa.stem.replace(".reheadered", "")
        seq_map[species_id] = {}
        for rec in SeqIO.parse(str(faa), "fasta"):
            seq_map[species_id][rec.id] = str(rec.seq)
    return seq_map


def load_orthologs_to_db(
    long_df: pd.DataFrame,
    seq_map: dict[str, dict[str, str]],
    human_gene_map: dict[str, tuple[str, str]],
    fresh_load: bool = False,
) -> int:
    """Insert orthologs into the database using bulk operations.

    Args:
        long_df: Output from parse_orthogroups().
        seq_map: Protein sequences per species.
        human_gene_map: {protein_id: (gene_id, gene_symbol)} for human proteins.
        fresh_load: If True, truncate Gene and Ortholog tables before loading.
                    Use when rerunning OrthoFinder from scratch on a new set of
                    species — otherwise stale rows from prior runs accumulate.

    Returns:
        Number of rows inserted.
    """
    if fresh_load:
        with get_session() as session:
            log.warning("fresh_load=True: truncating Gene and Ortholog tables before reload.")
            session.query(Ortholog).delete()
            session.query(Gene).delete()
            session.commit()
        log.info("Gene and Ortholog tables cleared.")

    with get_session() as session:
        existing_gene_ids = {g.id for g in session.query(Gene.id).all()}
        existing_orthologs = {
            (o.gene_id, o.species_id)
            for o in session.query(Ortholog.gene_id, Ortholog.species_id).all()
        }

    genes_to_add = []
    orthologs_to_add = []

    for og_id, group in long_df.groupby("og_id"):
        human_rows = group[group["species_id"] == "human"]
        if human_rows.empty:
            continue

        human_protein_id = human_rows.iloc[0]["protein_id"]
        if human_protein_id not in human_gene_map:
            continue

        gene_db_id, gene_symbol = human_gene_map[human_protein_id]

        # Resolve UniProt accession to HGNC gene symbol (e.g. P04637 → TP53).
        # human_gene_map usually stores the OrthoFinder entry name (AKT1_HUMAN) or
        # the bare UniProt accession. Convert to a human-readable HGNC symbol so
        # Gene.gene_symbol is useful in the UI and scoring reports.
        hgnc_symbol = _resolve_hgnc_symbol(gene_symbol or gene_db_id)

        if gene_db_id not in existing_gene_ids:
            genes_to_add.append(Gene(
                id=gene_db_id,
                human_gene_id=human_protein_id,
                gene_symbol=hgnc_symbol,
                human_protein=human_protein_id,
            ))
            existing_gene_ids.add(gene_db_id)

        for _, row in group.iterrows():
            sid = row["species_id"]
            pid = row["protein_id"]
            if (gene_db_id, sid) in existing_orthologs:
                continue
            seq = seq_map.get(sid, {}).get(pid, "")
            orthologs_to_add.append(Ortholog(
                gene_id=gene_db_id,
                species_id=sid,
                protein_id=pid,
                protein_seq=seq or None,
                orthofinder_og=og_id,
            ))
            existing_orthologs.add((gene_db_id, sid))

    inserted = 0
    _BATCH = 5000
    with get_session() as session:
        if genes_to_add:
            session.bulk_save_objects(genes_to_add)
            session.flush()
            log.info("  Bulk-inserted %d Gene rows.", len(genes_to_add))

        for i in range(0, len(orthologs_to_add), _BATCH):
            batch = orthologs_to_add[i : i + _BATCH]
            session.bulk_save_objects(batch)
            session.flush()
            inserted += len(batch)
            if inserted % 10000 == 0 or i + _BATCH >= len(orthologs_to_add):
                log.info("  Ortholog bulk insert: %d / %d", inserted, len(orthologs_to_add))

        session.commit()

    log.info("Inserted %d ortholog rows.", inserted)
    return inserted


def flag_one_to_one_orthogroups() -> int:
    """Mark orthogroups as 1:1 or not based on paralog presence.

    A 1:1 orthogroup has exactly one sequence per species that contributes.
    Any orthogroup where a species contributes >1 sequence is flagged
    is_one_to_one=False so downstream convergence/selection steps can
    filter it out with a simple WHERE clause.

    Returns number of orthogroups flagged as many-to-many.
    """
    many_to_many = 0
    with get_session() as session:
        # Group all orthologs by og_id and count per species
        from sqlalchemy import func as sqlfunc
        rows = (
            session.query(
                Ortholog.orthofinder_og,
                Ortholog.species_id,
                sqlfunc.count(Ortholog.id).label("n"),
            )
            .filter(Ortholog.orthofinder_og.isnot(None))
            .group_by(Ortholog.orthofinder_og, Ortholog.species_id)
            .all()
        )

        # Identify og_ids where any species has >1 sequence
        paralog_ogs: set[str] = set()
        for og_id, species_id, n in rows:
            if n > 1:
                paralog_ogs.add(og_id)

        if paralog_ogs:
            updated = (
                session.query(Ortholog)
                .filter(Ortholog.orthofinder_og.in_(paralog_ogs))
                .update({"is_one_to_one": False}, synchronize_session=False)
            )
            session.commit()
            many_to_many = len(paralog_ogs)
            log.info(
                "  1:1 ortholog filter: %d many-to-many orthogroups flagged (%d rows updated).",
                many_to_many, updated,
            )
        else:
            log.info("  1:1 ortholog filter: all orthogroups are 1:1.")

    return many_to_many
