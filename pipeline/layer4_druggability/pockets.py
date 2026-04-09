"""fpocket — pocket detection on PDB structures + convergent residue proximity.

Two functions:
  run_fpocket()          — runs fpocket, returns pocket_count + top_pocket_score.
  annotate_pockets()     — runs fpocket for each gene, writes DrugTarget, then
                           computes convergent_pocket_proximal by cross-referencing
                           fpocket's top-pocket residue list against ConvergentAA
                           positions for that gene.

Convergent pocket proximity logic:
  fpocket writes a {stem}_out/pockets/pocket1_atm.pdb (or pocket1_vert.pdb for
  alpha-sphere centroids).  We parse ATOM/HETATM records from the top-ranked
  pocket PDB to get residue numbers.  If any residue number matches a known
  convergent AA position (± PROXIMITY_RESIDUE_WINDOW), the gene is flagged
  convergent_pocket_proximal=True.

  We use a residue number window rather than an Å distance threshold because we
  are working in 1D sequence space here — a proper Å-distance check would require
  loading the full PDB structure (heavy dependency).  A window of ±5 residues is
  a conservative proxy; in real secondary structure elements this is 7–15 Å.
  If exact distance is needed later, mdanalysis / BioPython can be added.
"""

import logging
import re
import subprocess
from pathlib import Path
from typing import Optional

from db.models import DivergentMotif, DrugTarget, Gene, Ortholog
from db.session import get_session

log = logging.getLogger(__name__)

# Residue sequence window treated as "proximal" when we lack 3D coordinates.
# ±5 residues in a helix is ~7.5 Å; in a loop it may be up to ~20 Å.
# Conservative but computationally free.
PROXIMITY_RESIDUE_WINDOW = 5


def run_fpocket(pdb_path: Path, min_alpha_spheres: int = 3) -> Optional[dict]:
    """Run fpocket on a PDB file.

    Returns dict with:
      pocket_count        — number of pockets found
      top_pocket_score    — highest fpocket druggability score
      top_pocket_residues — set of residue numbers (int) in the top-scoring pocket
    """
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

    pocket_count = 0
    top_score = 0.0
    top_pocket_idx = None

    for search_dir in [out_dir, out_dir / "pockets"]:
        if not search_dir.exists():
            continue
        for f in sorted(search_dir.glob("pocket*_info.txt")):
            pocket_count += 1
            try:
                text = f.read_text()
                m = re.search(r"Druggability\s+Score[:\s]+([\d.]+)", text, re.I)
                if m:
                    s = float(m.group(1))
                    if s > top_score:
                        top_score = s
                        # Extract pocket index from filename: pocket3_info.txt → 3
                        idx_m = re.search(r"pocket(\d+)_info", f.name)
                        if idx_m:
                            top_pocket_idx = int(idx_m.group(1))
            except Exception:
                pass
        break

    # Parse residue numbers from the top pocket's atom PDB file
    top_pocket_residues: set[int] = set()
    if top_pocket_idx is not None:
        for search_dir in [out_dir / "pockets", out_dir]:
            if not search_dir.exists():
                continue
            # fpocket names pocket atom files: pocket{N}_atm.pdb
            atm_file = search_dir / f"pocket{top_pocket_idx}_atm.pdb"
            if not atm_file.exists():
                # Some versions use _vert.pdb (alpha-sphere vertices)
                atm_file = search_dir / f"pocket{top_pocket_idx}_vert.pdb"
            if atm_file.exists():
                top_pocket_residues = _parse_pdb_residue_numbers(atm_file)
                break

    return {
        "pocket_count": pocket_count,
        "top_pocket_score": round(top_score, 4) if top_score else None,
        "top_pocket_residues": top_pocket_residues,
    }


def _parse_pdb_residue_numbers(pdb_path: Path) -> set[int]:
    """Extract unique residue sequence numbers from ATOM/HETATM records in a PDB file."""
    residues: set[int] = set()
    try:
        for line in pdb_path.read_text(errors="replace").splitlines():
            if not line.startswith(("ATOM", "HETATM", "STP")):  # STP = fpocket alpha spheres
                continue
            # PDB columns 23–26 (1-indexed) = residue sequence number
            try:
                res_num = int(line[22:26].strip())
                residues.add(res_num)
            except (ValueError, IndexError):
                pass
    except Exception as exc:
        log.debug("PDB residue parse %s: %s", pdb_path.name, exc)
    return residues


def _convergent_positions_for_gene(gene_id: str) -> set[int]:
    """Return 1-indexed protein positions of all convergent motifs for this gene.

    Convergent amino acid positions are the midpoints of DivergentMotif records
    where convergent_aa_count > 0 (i.e. the same substitution appears in multiple
    independent lineages). We use the midpoint of (start_pos, end_pos) as the
    representative residue number. This is an approximation — the exact convergent
    position within the motif window is not individually stored.
    """
    with get_session() as session:
        rows = (
            session.query(DivergentMotif.start_pos, DivergentMotif.end_pos)
            .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
            .filter(
                Ortholog.gene_id == gene_id,
                DivergentMotif.convergent_aa_count > 0,
                DivergentMotif.start_pos.isnot(None),
            )
            .all()
        )
    positions: set[int] = set()
    for start, end in rows:
        if start is not None:
            midpoint = (start + end) // 2 if end is not None else start
            positions.add(int(midpoint))
    return positions


def _is_pocket_proximal(convergent_positions: set[int], pocket_residues: set[int]) -> bool:
    """Return True if any convergent position is within PROXIMITY_RESIDUE_WINDOW of a pocket residue."""
    if not convergent_positions or not pocket_residues:
        return False
    for conv_pos in convergent_positions:
        for pocket_res in pocket_residues:
            if abs(conv_pos - pocket_res) <= PROXIMITY_RESIDUE_WINDOW:
                return True
    return False


def annotate_pockets(gene_to_pdb: dict[str, Path]) -> int:
    """Run fpocket for each gene's structure, populate DrugTarget, and flag convergent pocket proximity."""
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

            # Convergent residue proximity to top pocket
            top_residues = result.get("top_pocket_residues", set())
            conv_positions = _convergent_positions_for_gene(gene_id)
            proximal = _is_pocket_proximal(conv_positions, top_residues)
            dt.convergent_pocket_proximal = proximal

            if conv_positions:
                log.debug(
                    "Pocket proximity %s: conv_positions=%s pocket_residues=%d proximal=%s",
                    gene_id, sorted(conv_positions), len(top_residues), proximal,
                )

            updated += 1

    log.info("fpocket: updated %d genes (pocket_count + convergent_pocket_proximal).", updated)
    return updated
