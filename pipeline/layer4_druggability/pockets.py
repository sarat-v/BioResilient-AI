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

import json
import logging
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from db.models import DivergentMotif, DrugTarget, Gene, Ortholog
from db.session import get_session

log = logging.getLogger(__name__)


def _fpocket_cache_path(pdb_path: Path) -> Path:
    """Sidecar JSON cache alongside the PDB file."""
    return pdb_path.parent / f"{pdb_path.stem}_fpocket_cache.json"


def _load_fpocket_cache(pdb_path: Path) -> Optional[dict]:
    """Return cached fpocket result dict or None if cache is absent / stale."""
    cache_path = _fpocket_cache_path(pdb_path)
    if not cache_path.exists():
        return None
    try:
        cache_mtime = cache_path.stat().st_mtime
        pdb_mtime   = pdb_path.stat().st_mtime
        if pdb_mtime > cache_mtime:
            return None   # PDB newer than cache — re-run
        data = json.loads(cache_path.read_text())
        # JSON stores residue list; convert back to set[int]
        data["top_pocket_residues"] = set(data.get("top_pocket_residues") or [])
        return data
    except Exception:
        return None


def _save_fpocket_cache(pdb_path: Path, result: dict) -> None:
    """Write fpocket result to sidecar cache (residues serialised as sorted list)."""
    cache_path = _fpocket_cache_path(pdb_path)
    try:
        serialisable = dict(result)
        serialisable["top_pocket_residues"] = sorted(result.get("top_pocket_residues") or [])
        cache_path.write_text(json.dumps(serialisable))
    except Exception as exc:
        log.debug("fpocket cache write failed for %s: %s", pdb_path.name, exc)

# Residue sequence window treated as "proximal" when we lack 3D coordinates.
# ±5 residues in a helix is ~7.5 Å; in a loop it may be up to ~20 Å.
# Conservative but computationally free.
PROXIMITY_RESIDUE_WINDOW = 5

# AlphaFold stores pLDDT in the B-factor column.  Residues below this threshold
# are disordered and create noise that prevents fpocket from detecting pockets.
MIN_PLDDT = 50.0


def _filter_alphafold_pdb(pdb_path: Path, min_plddt: float = MIN_PLDDT) -> Path:
    """Return a pLDDT-filtered copy of an AlphaFold PDB file.

    Writes a temporary file containing only ATOM/HETATM lines where the
    B-factor (= pLDDT) is >= min_plddt.  This prevents disordered tails
    from blocking fpocket pocket detection.  Non-ATOM lines (HEADER, REMARK,
    TER, END …) are always kept.

    Returns the original path unchanged if no filtering was needed or if the
    file doesn't look like an AlphaFold structure.
    """
    try:
        lines = pdb_path.read_text(errors="replace").splitlines(keepends=True)
    except Exception:
        return pdb_path

    filtered = []
    removed = 0
    for line in lines:
        record = line[:6].strip()
        if record in ("ATOM", "HETATM"):
            try:
                bfactor = float(line[60:66])
            except (ValueError, IndexError):
                bfactor = 100.0  # keep if B-factor unreadable
            if bfactor < min_plddt:
                removed += 1
                continue
        filtered.append(line)

    if removed == 0:
        return pdb_path  # nothing to filter

    tmp = tempfile.NamedTemporaryFile(
        suffix=".pdb", dir=str(pdb_path.parent), delete=False
    )
    tmp.write("".join(filtered).encode())
    tmp.close()
    log.debug(
        "pLDDT filter %s: removed %d low-confidence residues (pLDDT < %g)",
        pdb_path.name, removed, min_plddt,
    )
    return Path(tmp.name)


def run_fpocket(pdb_path: Path, min_alpha_spheres: int = 3) -> Optional[dict]:
    """Run fpocket on a PDB file.

    Returns dict with:
      pocket_count        — number of pockets found
      top_pocket_score    — highest fpocket druggability score
      top_pocket_residues — set of residue numbers (int) in the top-scoring pocket

    Results are cached in a sidecar JSON file so that step 12 can re-use
    step 9b's computation without re-running fpocket on the same PDB.
    """
    if not pdb_path.exists() or pdb_path.stat().st_size < 100:
        return None

    # Return cached result when PDB hasn't changed since last run
    cached = _load_fpocket_cache(pdb_path)
    if cached is not None:
        log.debug("fpocket cache hit: %s", pdb_path.name)
        return cached

    # Filter out low-pLDDT residues (disordered tails block pocket detection).
    filtered_path = _filter_alphafold_pdb(pdb_path)
    tmp_created = filtered_path != pdb_path  # track whether we need to clean up

    # Do NOT pre-create out_dir — fpocket creates it itself.
    out_dir = filtered_path.parent / f"{filtered_path.stem}_out"

    try:
        proc = subprocess.run(
            ["fpocket", "-f", str(filtered_path), "-m", str(min_alpha_spheres)],
            cwd=str(filtered_path.parent),
            capture_output=True,
            check=False,      # don't raise — some builds exit 1 on warnings but still write output
            timeout=120,
        )
        if proc.returncode != 0:
            stderr_snippet = (proc.stderr or b"")[-400:].decode(errors="replace").strip()
            log.debug("fpocket %s exit %d: %s", filtered_path.name, proc.returncode, stderr_snippet)
    except FileNotFoundError as exc:
        log.warning("fpocket not found: %s", exc)
        if tmp_created:
            filtered_path.unlink(missing_ok=True)
        return None
    except subprocess.TimeoutExpired:
        log.debug("fpocket %s timed out", filtered_path.name)
        if tmp_created:
            filtered_path.unlink(missing_ok=True)
        return None

    import os as _os

    # ------------------------------------------------------------------ #
    # fpocket 3.x output format (bioconda):                               #
    #   {stem}_info.txt          — summary with all pocket scores         #
    #   pockets/pocket{N}_atm.pdb — atom coords for pocket N              #
    #   pockets/pocket{N}_vert.pqr — alpha-sphere centroids               #
    #                                                                      #
    # Older format (< 3.x):                                               #
    #   pocket{N}_info.txt        — per-pocket info file                  #
    #   pocket{N}_atm.pdb                                                  #
    # ------------------------------------------------------------------ #

    pocket_count = 0
    top_score = 0.0
    top_pocket_idx = None

    if out_dir.exists():
        # --- Strategy 1: fpocket 3.x single summary file ---
        summary_file = out_dir / f"{filtered_path.stem}_info.txt"
        if summary_file.exists():
            try:
                text = summary_file.read_text()
                # Count "Pocket N :" sections
                pocket_headers = re.findall(r"^Pocket\s+(\d+)\s*:", text, re.M)
                pocket_count = len(pocket_headers)
                # Per-pocket druggability scores: "Pocket N :\n\t...\n\tDruggability Score : X"
                for m in re.finditer(
                    r"Pocket\s+(\d+)\s*:.*?Druggability\s+Score\s*:\s*([\d.]+)",
                    text, re.S | re.I
                ):
                    idx, score = int(m.group(1)), float(m.group(2))
                    if score > top_score:
                        top_score = score
                        top_pocket_idx = idx
            except Exception:
                pass

        # --- Strategy 2: older fpocket — per-pocket info files ---
        if pocket_count == 0:
            for search_dir in [out_dir / "pockets", out_dir]:
                if not search_dir.exists():
                    continue
                per_pocket_files = sorted(search_dir.glob("pocket*_info.txt"))
                if not per_pocket_files:
                    continue
                for f in per_pocket_files:
                    pocket_count += 1
                    try:
                        t = f.read_text()
                        m = re.search(r"Druggability\s+Score[:\s]+([\d.]+)", t, re.I)
                        if m:
                            s = float(m.group(1))
                            if s > top_score:
                                top_score = s
                                idx_m = re.search(r"pocket(\d+)_info", f.name)
                                if idx_m:
                                    top_pocket_idx = int(idx_m.group(1))
                    except Exception:
                        pass
                break

        # --- Fallback: count pocket*_atm.pdb files if info files missing ---
        if pocket_count == 0:
            for search_dir in [out_dir / "pockets", out_dir]:
                atm_files = sorted(search_dir.glob("pocket*_atm.pdb")) if search_dir.exists() else []
                if atm_files:
                    pocket_count = len(atm_files)
                    top_pocket_idx = 1  # default to pocket 1
                    break

    log.debug("fpocket %s: pocket_count=%d top_score=%s", filtered_path.name, pocket_count, top_score)

    # Parse residue numbers from the top pocket's atom PDB file
    top_pocket_residues: set[int] = set()
    if top_pocket_idx is not None and out_dir.exists():
        for search_dir in [out_dir / "pockets", out_dir]:
            if not search_dir.exists():
                continue
            atm_file = search_dir / f"pocket{top_pocket_idx}_atm.pdb"
            if atm_file.exists():
                top_pocket_residues = _parse_pdb_residue_numbers(atm_file)
                break

    if tmp_created:
        filtered_path.unlink(missing_ok=True)

    result = {
        "pocket_count": pocket_count,
        "top_pocket_score": round(top_score, 4) if top_score else None,
        "top_pocket_residues": top_pocket_residues,
    }
    _save_fpocket_cache(pdb_path, result)
    return result


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
    for gene_id, pdb_path in gene_to_pdb.items():
        result = run_fpocket(pdb_path)
        if result is None:
            continue

        conv_positions = _convergent_positions_for_gene(gene_id)
        top_residues = result.get("top_pocket_residues", set())
        proximal = _is_pocket_proximal(conv_positions, top_residues)

        if conv_positions:
            log.debug(
                "Pocket proximity %s: conv_positions=%s pocket_residues=%d proximal=%s",
                gene_id, sorted(conv_positions), len(top_residues), proximal,
            )

        with get_session() as session:
            dt = session.get(DrugTarget, gene_id)
            if dt is None:
                dt = DrugTarget(gene_id=gene_id)
                session.add(dt)

            dt.pocket_count = result["pocket_count"]
            if result.get("top_pocket_score") is not None:
                dt.top_pocket_score = result["top_pocket_score"]
            dt.convergent_pocket_proximal = proximal
            session.commit()

        updated += 1
        log.debug("fpocket %s: %d pockets, score=%s, proximal=%s",
                  pdb_path.name, result["pocket_count"], result.get("top_pocket_score"), proximal)

    log.info("fpocket: updated %d genes (pocket_count + convergent_pocket_proximal).", updated)
    return updated
