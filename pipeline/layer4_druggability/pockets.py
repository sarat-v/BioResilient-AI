"""fpocket — pocket detection on PDB structures."""

import logging
import re
import subprocess
from pathlib import Path
from typing import Optional

from db.models import DrugTarget, Gene
from db.session import get_session

log = logging.getLogger(__name__)


def run_fpocket(pdb_path: Path, min_alpha_spheres: int = 3) -> Optional[dict]:
    """Run fpocket on a PDB file. Returns dict with pocket_count and top_pocket_score."""
    if not pdb_path.exists() or pdb_path.stat().st_size < 100:
        return None
    out_dir = pdb_path.parent / f"{pdb_path.stem}_out"
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["fpocket", "-f", str(pdb_path), "-m", str(min_alpha_spheres)],
            cwd=str(pdb_path.parent),
            capture_output=True,
            check=True,
            timeout=120,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        log.debug("fpocket %s: %s", pdb_path.name, exc)
        return None

    # Parse fpocket output: {stem}_out/pocket*_info.txt or {stem}_out/pockets/pocket*_info.txt
    pocket_count = 0
    top_score = 0.0
    for search_dir in [out_dir, out_dir / "pockets"]:
        if not search_dir.exists():
            continue
        for f in search_dir.glob("pocket*_info.txt"):
            pocket_count += 1
            try:
                text = f.read_text()
                m = re.search(r"Druggability\s+Score[:\s]+([\d.]+)", text, re.I)
                if m:
                    s = float(m.group(1))
                    if s > top_score:
                        top_score = s
            except Exception:
                pass
        break

    return {"pocket_count": pocket_count, "top_pocket_score": round(top_score, 4) if top_score else None}


def annotate_pockets(gene_to_pdb: dict[str, Path]) -> int:
    """Run fpocket for each gene's structure and populate DrugTarget."""
    updated = 0
    with get_session() as session:
        for gene_id, pdb_path in gene_to_pdb.items():
            result = run_fpocket(pdb_path)
            if result is None:
                continue
            dt = session.get(DrugTarget, gene_id)
            if dt is None:
                dt = DrugTarget(gene_id=gene_id)
                session.add(dt)
            dt.pocket_count = result["pocket_count"]
            if result.get("top_pocket_score") is not None:
                dt.top_pocket_score = result["top_pocket_score"]
            updated += 1
    log.info("fpocket: updated %d genes.", updated)
    return updated
