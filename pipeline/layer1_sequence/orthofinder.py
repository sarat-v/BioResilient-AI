"""Step 3 — OrthoFinder wrapper.

Runs OrthoFinder on all downloaded proteomes simultaneously using the
DIAMOND backend (-S diamond), which is 500x faster than BLAST.

After the run, parses OrthoGroups.tsv and loads into the Ortholog table,
keeping only groups where the human protein is present.
"""

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO

from db.models import Gene, Ortholog
from db.session import get_session
from pipeline.config import get_storage_root, get_tool_config

log = logging.getLogger(__name__)


def _orthofinder_output_dir(proteomes_dir: Path) -> Optional[Path]:
    """OrthoFinder creates a dated output folder — locate the most recent one."""
    results_dir = proteomes_dir.parent / "orthofinder_out" / "Results_*"
    import glob
    matches = sorted(glob.glob(str(results_dir)))
    if matches:
        return Path(matches[-1])
    return None


def run_orthofinder(proteomes_dir: Path) -> Path:
    """Run OrthoFinder with DIAMOND on the proteomes directory.

    If a completed results directory already exists, reuse it (cache).

    Args:
        proteomes_dir: Directory containing one .faa file per species (reheadered).

    Returns:
        Path to the OrthoFinder results directory.
    """
    cfg = get_tool_config()
    search_threads = cfg.get("orthofinder_threads", 8)
    align_threads = cfg.get("orthofinder_align_threads", 4)

    root = Path(get_storage_root())
    out_dir = root / "orthofinder_out"

    # Reuse cached results if they exist
    existing = _orthofinder_output_dir(proteomes_dir)
    if existing is not None and (existing / "Orthogroups" / "Orthogroups.tsv").exists():
        log.info("OrthoFinder results already exist — reusing cache: %s", existing)
        return existing

    # Remove stale/incomplete output dir so OrthoFinder can create a fresh one
    if out_dir.exists():
        log.info("Removing stale OrthoFinder output dir: %s", out_dir)
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "orthofinder",
        "-f", str(proteomes_dir),
        "-S", "diamond",
        "-t", str(search_threads),
        "-a", str(align_threads),
        "-o", str(out_dir),
        "-y",   # split paralogous clades — produces cleaner orthogroups
    ]

    log.info("Running OrthoFinder: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=False, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"OrthoFinder exited with code {result.returncode}")

    results_path = _orthofinder_output_dir(proteomes_dir)
    if results_path is None:
        raise FileNotFoundError("OrthoFinder output directory not found after run.")

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

    # Each column is a species, values are comma-separated protein IDs
    df.index.name = "og_id"
    df = df.reset_index()

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
                    # Column name is the species FASTA filename stem; map to species_id
                    species_id = _col_to_species_id(col, proteomes_dir)
                    long_rows.append({
                        "og_id": og_id,
                        "species_id": species_id,
                        "protein_id": protein_id,
                    })

    return pd.DataFrame(long_rows)


def _find_species_column(df: pd.DataFrame, species_id: str) -> Optional[str]:
    """Match a species_id to a DataFrame column (flexible matching on stem)."""
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
    """Build map: {species_id: {protein_id: sequence}} from FASTA files."""
    seq_map: dict[str, dict[str, str]] = {}
    for faa in proteomes_dir.glob("*.faa"):
        species_id = faa.stem.replace(".reheadered", "")
        seq_map[species_id] = {}
        for rec in SeqIO.parse(str(faa), "fasta"):
            seq_map[species_id][rec.id] = str(rec.seq)
    return seq_map


def load_orthologs_to_db(
    long_df: pd.DataFrame,
    seq_map: dict[str, dict[str, str]],
    human_gene_map: dict[str, tuple[str, str]],
) -> int:
    """Insert orthologs into the database.

    Args:
        long_df: Output from parse_orthogroups().
        seq_map: Protein sequences per species.
        human_gene_map: {protein_id: (gene_id, gene_symbol)} for human proteins.

    Returns:
        Number of rows inserted.
    """
    inserted = 0

    with get_session() as session:
        # Group by og_id to find human anchor
        for og_id, group in long_df.groupby("og_id"):
            human_rows = group[group["species_id"] == "human"]
            if human_rows.empty:
                continue

            # Use the first human protein as anchor gene
            human_protein_id = human_rows.iloc[0]["protein_id"]
            if human_protein_id not in human_gene_map:
                continue

            gene_db_id, gene_symbol = human_gene_map[human_protein_id]

            # Ensure Gene row exists
            gene = session.get(Gene, gene_db_id)
            if gene is None:
                gene = Gene(
                    id=gene_db_id,
                    human_gene_id=human_protein_id,
                    gene_symbol=gene_symbol,
                    human_protein=human_protein_id,
                )
                session.add(gene)
                session.flush()

            # Insert orthologs for all species
            for _, row in group.iterrows():
                sid = row["species_id"]
                pid = row["protein_id"]
                seq = seq_map.get(sid, {}).get(pid, "")

                existing = (
                    session.query(Ortholog)
                    .filter_by(gene_id=gene_db_id, species_id=sid)
                    .first()
                )
                if existing:
                    continue

                ortholog = Ortholog(
                    gene_id=gene_db_id,
                    species_id=sid,
                    protein_id=pid,
                    protein_seq=seq or None,
                    orthofinder_og=og_id,
                )
                session.add(ortholog)
                inserted += 1

    log.info("Inserted %d ortholog rows.", inserted)
    return inserted
