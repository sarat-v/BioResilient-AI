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
from pipeline.config import get_storage_root, get_tool_config

log = logging.getLogger(__name__)


def _p2rank_output_dir(pdb_path: Path) -> Path:
    root = Path(get_storage_root()) / "p2rank"
    root.mkdir(parents=True, exist_ok=True)
    return root / pdb_path.stem


def run_p2rank(pdb_path: Path) -> Optional[list[dict]]:
    """Run P2Rank on a PDB file and return parsed pocket predictions.

    Returns list of pocket dicts: [{"rank": int, "score": float, "probability": float}, ...]
    sorted by rank (best first). Returns None if P2Rank is unavailable or fails.
    """
    tools = get_tool_config()
    p2rank_jar = tools.get("p2rank_jar")

    if not p2rank_jar:
        log.debug("p2rank_jar not configured — skipping P2Rank.")
        return None
    if not Path(p2rank_jar).exists():
        log.debug("P2Rank JAR not found at %s — skipping.", p2rank_jar)
        return None

    out_dir = _p2rank_output_dir(pdb_path)
    predictions_csv = out_dir / f"{pdb_path.stem}.pdb_predictions.csv"

    # Use cached result if available
    if predictions_csv.exists():
        return _parse_p2rank_predictions(predictions_csv)

    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        cmd = [
            "java", "-jar", p2rank_jar,
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
            log.debug("P2Rank failed for %s: %s", pdb_path.name, result.stderr[:300])
            return None

        if not predictions_csv.exists():
            # P2Rank may use a slightly different output filename
            csvs = list(out_dir.glob("*.csv"))
            if csvs:
                predictions_csv = csvs[0]
            else:
                return None

        return _parse_p2rank_predictions(predictions_csv)

    except subprocess.TimeoutExpired:
        log.debug("P2Rank timed out for %s", pdb_path.name)
        return None
    except FileNotFoundError:
        log.debug("Java not found — P2Rank requires Java 11+.")
        return None
    except Exception as exc:
        log.debug("P2Rank error: %s", exc)
        return None


def _parse_p2rank_predictions(csv_path: Path) -> list[dict]:
    """Parse P2Rank _predictions.csv output file.

    Columns: name, rank, score, probability, sas_points, surf_atoms, ...
    """
    pockets = []
    try:
        with open(csv_path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    pockets.append({
                        "rank": int(row.get("rank", 0)),
                        "score": float(row.get("score", 0)),
                        "probability": float(row.get("probability", 0)),
                    })
                except (ValueError, TypeError):
                    continue
        pockets.sort(key=lambda p: p["rank"])
    except Exception as exc:
        log.debug("P2Rank CSV parse error: %s", exc)
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
