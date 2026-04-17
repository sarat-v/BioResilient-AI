"""Step 12b — P2Rank ML pocket prediction for druggability.

P2Rank (Krivak & Hoksza, 2018) is a machine-learning-based tool for predicting
protein-ligand binding pockets from 3D structure. Unlike fpocket (geometry-based),
P2Rank uses a random forest trained on known binding sites, giving more reliable
druggability predictions for non-obvious pockets.

For each Tier1/Tier2 gene with an AlphaFold structure:
  1. Run P2Rank on the PDB file.
  2. Parse the output CSV to extract pocket scores.
  3. Store p2rank_score (best pocket probability) and p2rank_pocket_count on DrugTarget.

Requirements:
  - P2Rank JAR: download from https://github.com/rdk/p2rank/releases
  - Java 11+ runtime
  - Configure tool path in environment.yml: p2rank_jar: /path/to/p2rank.jar
  - AlphaFold PDB files must exist (Step 12a downloads them)

P2Rank scores [0–1]: closer to 1 = higher probability of druggable pocket.
"""

import csv
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from db.models import DrugTarget, Gene
from db.session import get_session
from pipeline.config import get_local_storage_root, get_tool_config

log = logging.getLogger(__name__)


def _p2rank_output_dir(pdb_path: Path) -> Path:
    root = Path(get_local_storage_root()) / "p2rank"
    root.mkdir(parents=True, exist_ok=True)
    return root / pdb_path.stem


def run_p2rank(pdb_path: Path) -> Optional[list[dict]]:
    """Run P2Rank on a PDB file and return parsed pocket predictions.

    Returns list of pocket dicts: [{"rank": int, "score": float, "probability": float}, ...]
    sorted by rank (best first). Returns None if P2Rank is unavailable or fails.
    """
    import os as _os
    tools = get_tool_config()
    # Use the prank shell script (sets full classpath with Groovy + lib/*).
    # Check YAML config first, then P2RANK_BIN env var (set in clinical Dockerfile).
    prank_bin = tools.get("prank_bin") or _os.environ.get("P2RANK_BIN", "")
    if not prank_bin:
        # Fallback: derive prank script path from P2RANK_JAR env var
        jar = tools.get("p2rank_jar") or _os.environ.get("P2RANK_JAR", "")
        if jar:
            prank_bin = str(Path(jar).parent.parent / "prank")

    log.info("P2Rank prank script: %r  (P2RANK_BIN env=%r)", prank_bin, _os.environ.get("P2RANK_BIN"))
    if not prank_bin:
        log.warning("prank_bin not configured — skipping P2Rank.")
        return None
    if not Path(prank_bin).exists():
        log.warning("prank script not found at %s — skipping.", prank_bin)
        return None

    out_dir = _p2rank_output_dir(pdb_path)
    predictions_csv = out_dir / f"{pdb_path.stem}.pdb_predictions.csv"

    # Use cached result if available
    if predictions_csv.exists():
        return _parse_p2rank_predictions(predictions_csv)

    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        cmd = [
            str(prank_bin),
            "predict", "-f", str(pdb_path),
            "-o", str(out_dir),
        ]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0:
            log.warning("P2Rank failed for %s (rc=%d):\nstdout: %s\nstderr: %s",
                        pdb_path.name, result.returncode, result.stdout[:200], result.stderr[:300])
            return None

        if not predictions_csv.exists():
            # P2Rank may use a slightly different output filename
            csvs = list(out_dir.glob("*.csv"))
            log.info("P2Rank output CSVs for %s: %s", pdb_path.name, [c.name for c in csvs])
            if csvs:
                predictions_csv = csvs[0]
            else:
                log.warning("P2Rank produced no CSV for %s; out_dir contents: %s",
                            pdb_path.name, [f.name for f in out_dir.iterdir()] if out_dir.exists() else "dir missing")
                return None

        pockets = _parse_p2rank_predictions(predictions_csv)
        log.info("P2Rank parsed %d pockets for %s", len(pockets), pdb_path.name)
        return pockets

    except subprocess.TimeoutExpired:
        log.warning("P2Rank timed out for %s", pdb_path.name)
        return None
    except FileNotFoundError:
        log.warning("Java not found — P2Rank requires Java 11+.")
        return None
    except Exception as exc:
        log.warning("P2Rank error for %s: %s", pdb_path.name, exc)
        return None


def _parse_p2rank_predictions(csv_path: Path) -> list[dict]:
    """Parse P2Rank _predictions.csv output file.

    P2Rank uses space-padded column headers (e.g., " probability", " rank").
    We strip whitespace from all keys and values before reading.

    Columns: name, rank, score, probability, sas_points, surf_atoms, ...
    """
    pockets = []
    try:
        with open(csv_path, newline="") as f:
            reader = csv.DictReader(f)
            # Normalise header names: strip surrounding whitespace
            if reader.fieldnames:
                reader.fieldnames = [h.strip() for h in reader.fieldnames]
            _header_logged = False
            for row in reader:
                # Strip whitespace from values too
                row = {k.strip(): v.strip() if isinstance(v, str) else v for k, v in row.items()}
                if not _header_logged:
                    log.info("P2Rank CSV columns: %s", list(row.keys()))
                    _header_logged = True
                try:
                    pockets.append({
                        "rank": int(row.get("rank") or 0),
                        "score": float(row.get("score") or 0),
                        "probability": float(row.get("probability") or 0),
                    })
                except (ValueError, TypeError):
                    continue
        pockets.sort(key=lambda p: p["rank"])
    except Exception as exc:
        log.warning("P2Rank CSV parse error for %s: %s", csv_path.name, exc)
    return pockets


def annotate_p2rank(
    gene_to_pdb: dict[str, Path],
    gene_ids: Optional[list[str]] = None,
) -> int:
    """Run P2Rank on AlphaFold structures and update DrugTarget rows.

    Args:
        gene_to_pdb: {gene_id: pdb_path} — from the structure step (step12a).
        gene_ids: Optional filter to specific gene IDs.

    Returns:
        Number of genes with P2Rank scores stored.
    """
    if gene_ids:
        gene_to_pdb = {g: p for g, p in gene_to_pdb.items() if g in gene_ids}

    log.info("Running P2Rank on %d structures...", len(gene_to_pdb))
    scored = 0

    with get_session() as session:
        for gene_id, pdb_path in gene_to_pdb.items():
            if not pdb_path or not Path(pdb_path).exists():
                continue

            pockets = run_p2rank(Path(pdb_path))
            if pockets is None:
                continue

            dt = session.get(DrugTarget, gene_id)
            if dt is None:
                dt = DrugTarget(gene_id=gene_id)
                session.add(dt)

            dt.p2rank_pocket_count = len(pockets)
            dt.p2rank_score = pockets[0]["probability"] if pockets else 0.0
            scored += 1

        session.commit()

    log.info("P2Rank annotation complete: %d genes scored.", scored)
    return scored


def run_p2rank_pipeline(gene_to_pdb: dict, gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step12b."""
    return annotate_p2rank(gene_to_pdb, gene_ids)
