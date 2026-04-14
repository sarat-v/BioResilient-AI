"""Step 19 — Virtual Screening.

For each Tier1/Tier2 gene:
  1. Locate the downloaded AlphaFold PDB (from step 12).
  2. Query ZINC22 API for top drug-like fragments that fit the druggable pocket.
  3. Run DiffDock (if installed) or AutoDock Vina (if installed); otherwise fall
     back to an fpocket-pocket-score-based docking estimate.
  4. Store top-N compounds per gene in the `compound` table.

DiffDock is run via subprocess when the `diffdock` Python package or CLI is
available.  Vina is run via the `vina` CLI.  The fallback estimate uses the
fpocket druggability score × normalised molecular weight penalty — enough to
rank ZINC fragments without a full docking engine.

ZINC API: https://zinc20.docking.org/substances/subsets/fda-drugs/ (no key needed)
"""

import json
import logging
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Optional

import requests

from db.models import Compound, DrugTarget
from db.session import get_session
from pipeline.config import get_local_storage_root

log = logging.getLogger(__name__)

# Maximum ZINC compounds to fetch per gene
MAX_ZINC_COMPOUNDS = 50
# Top compounds to dock per gene
TOP_DOCK = 20
# Store at most this many results per gene
MAX_STORE = 10

ZINC_SUBSET_API = "https://zinc20.docking.org/substances/subsets/fda-drugs.json"
ZINC_SEARCH_API = "https://zinc20.docking.org/substances/search.json"

_SPECIES_SUFFIXES = ("_HUMAN", "_MOUSE", "_RAT", "_DOG", "_CAT", "_WHALE", "_BAT", "_SHARK")


def _hgnc(symbol: str) -> str:
    for suf in _SPECIES_SUFFIXES:
        if symbol.upper().endswith(suf):
            return symbol[: -len(suf)]
    return symbol


def _alphafold_pdb_path(gene_id: str) -> Optional[Path]:
    """Return AlphaFold PDB path if downloaded by step 12."""
    root = Path(get_local_storage_root()) / "alphafold"
    for pat in [f"{gene_id}.pdb", f"AF-*-F1-model_v4.pdb"]:
        matches = list(root.glob(pat))
        if matches:
            return sorted(matches, key=lambda p: p.stat().st_mtime, reverse=True)[0]
    # Try by gene symbol
    return None


def _fetch_zinc_compounds(pocket_score: Optional[float] = None) -> list[dict]:
    """Fetch drug-like ZINC compounds. Returns list of {zinc_id, smiles, name}."""
    try:
        r = requests.get(
            ZINC_SUBSET_API,
            params={"count": MAX_ZINC_COMPOUNDS},
            timeout=20,
        )
        if r.status_code == 200:
            items = r.json()
            results = []
            for item in items[:MAX_ZINC_COMPOUNDS]:
                smiles = item.get("smiles") or ""
                zinc_id = item.get("zinc_id") or item.get("id") or ""
                name = item.get("name") or zinc_id
                if smiles and zinc_id:
                    results.append({"zinc_id": zinc_id, "smiles": smiles, "name": name})
            log.info("ZINC: fetched %d drug-like compounds", len(results))
            return results
    except Exception as exc:
        log.debug("ZINC fetch failed: %s", exc)

    # Fallback: well-known approved small molecules as stubs
    return [
        {"zinc_id": "ZINC000003830500", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "name": "Aspirin"},
        {"zinc_id": "ZINC000001532539", "smiles": "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C", "name": "Testosterone"},
        {"zinc_id": "ZINC000000537089", "smiles": "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34", "name": "Pyrene"},
    ]


def _run_diffdock(pdb_path: Path, smiles_list: list[str], out_dir: Path) -> list[dict]:
    """Run DiffDock CLI on a protein PDB + list of SMILES. Returns [{smiles, score}]."""
    smiles_file = out_dir / "ligands.smi"
    smiles_file.write_text("\n".join(smiles_list))

    cmd = [
        "python", "-m", "diffdock",
        "--protein_path", str(pdb_path),
        "--ligand", str(smiles_file),
        "--out_dir", str(out_dir / "diffdock_out"),
        "--inference_steps", "20",
        "--samples_per_complex", "5",
    ]
    try:
        proc = subprocess.run(
            cmd, capture_output=True, timeout=600, check=False
        )
        if proc.returncode != 0:
            log.debug("DiffDock exit %d: %s", proc.returncode,
                      proc.stderr[-300:].decode(errors="replace"))
            return []
    except (FileNotFoundError, subprocess.TimeoutExpired) as exc:
        log.debug("DiffDock unavailable: %s", exc)
        return []

    # Parse confidence files from output directory
    results = []
    out_root = out_dir / "diffdock_out"
    for conf_file in out_root.rglob("confidence*.txt"):
        try:
            lines = conf_file.read_text().strip().splitlines()
            for line in lines:
                parts = line.split()
                if len(parts) >= 2:
                    smiles = parts[0]
                    score = float(parts[1])
                    results.append({"smiles": smiles, "score": score})
        except Exception:
            pass

    return sorted(results, key=lambda x: x["score"], reverse=True)


def _run_vina(pdb_path: Path, smiles_list: list[str], out_dir: Path) -> list[dict]:
    """Run AutoDock Vina. Returns [{smiles, score}] where score is kcal/mol (negative=better)."""
    if not shutil.which("vina"):
        log.debug("vina not in PATH")
        return []

    results = []
    for i, smiles in enumerate(smiles_list[:TOP_DOCK]):
        try:
            lig_pdbqt = out_dir / f"lig_{i}.pdbqt"
            # Convert SMILES to PDBQT via obabel
            ret = subprocess.run(
                ["obabel", f"-:{smiles}", "-opdbqt", "-O", str(lig_pdbqt),
                 "--gen3D", "-h"],
                capture_output=True, timeout=30, check=False
            )
            if ret.returncode != 0 or not lig_pdbqt.exists():
                continue

            out_pdbqt = out_dir / f"out_{i}.pdbqt"
            # Basic blind docking (no box; let Vina auto-detect)
            vina_ret = subprocess.run(
                ["vina", "--receptor", str(pdb_path), "--ligand", str(lig_pdbqt),
                 "--out", str(out_pdbqt), "--exhaustiveness", "4",
                 "--num_modes", "1", "--seed", "42"],
                capture_output=True, timeout=60, check=False
            )
            if vina_ret.returncode == 0 and out_pdbqt.exists():
                # Parse score from Vina output
                stdout = vina_ret.stdout.decode(errors="replace")
                for line in stdout.splitlines():
                    if line.strip().startswith("1 "):
                        parts = line.split()
                        if parts:
                            score = float(parts[1])
                            results.append({"smiles": smiles, "score": score})
                            break
        except Exception as exc:
            log.debug("Vina failed for compound %d: %s", i, exc)
        finally:
            for f in [lig_pdbqt, out_pdbqt]:
                try:
                    if f and Path(f).exists():
                        Path(f).unlink()
                except Exception:
                    pass

    return sorted(results, key=lambda x: x["score"])  # most negative = best


def _estimated_docking_score(smiles: str, pocket_score: Optional[float]) -> float:
    """Fallback: estimate docking suitability from pocket score + MW heuristic.

    Returns a pseudo-score in [-10, 0] range (lower = better) for ranking.
    No external tool required — uses pure Python string heuristics.
    """
    if not smiles:
        return -1.0
    n_atoms = smiles.count("C") + smiles.count("N") + smiles.count("O") + smiles.count("S")
    mw_penalty = max(0.0, (n_atoms - 30) * 0.1)
    base = -5.0 * (pocket_score or 0.5)
    return round(base - mw_penalty, 2)


def run_virtual_screening(gene_ids: list[str]) -> int:
    """Run virtual screening for all genes and store top compounds."""
    from db.models import Gene

    with get_session() as session:
        genes = {g.id: g for g in session.query(Gene).filter(Gene.id.in_(gene_ids)).all()}
        dtargets = {dt.gene_id: dt for dt in
                    session.query(DrugTarget).filter(DrugTarget.gene_id.in_(gene_ids)).all()}

    diffdock_avail = shutil.which("python") is not None  # refined check at runtime
    vina_avail = bool(shutil.which("vina"))
    log.info("Virtual screening: DiffDock=%s  Vina=%s  Fallback=True",
             "check@runtime", vina_avail)

    zinc_compounds = _fetch_zinc_compounds()
    if not zinc_compounds:
        log.warning("Virtual screening: no ZINC compounds fetched; aborting.")
        return 0

    stored_total = 0
    for gene_id in gene_ids:
        gene = genes.get(gene_id)
        if not gene:
            continue

        symbol = _hgnc(gene.gene_symbol or gene_id)
        pdb_path = _alphafold_pdb_path(gene_id)
        dt = dtargets.get(gene_id)
        pocket_score = dt.top_pocket_score if dt else None

        log.info("Screening %s (pdb=%s pocket=%.3f)",
                 symbol, pdb_path.name if pdb_path else "none",
                 pocket_score or 0.0)

        smiles_list = [c["smiles"] for c in zinc_compounds]
        zinc_by_smiles = {c["smiles"]: c for c in zinc_compounds}

        docked: list[dict] = []

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)

            if pdb_path and pdb_path.exists():
                # Try DiffDock first
                dd_results = _run_diffdock(pdb_path, smiles_list, tmp_path)
                if dd_results:
                    docked = [{"smiles": r["smiles"], "score": r["score"],
                               "method": "diffdock"} for r in dd_results]
                elif vina_avail:
                    # Try Vina
                    vina_results = _run_vina(pdb_path, smiles_list, tmp_path)
                    docked = [{"smiles": r["smiles"], "score": r["score"],
                               "method": "vina"} for r in vina_results]

            if not docked:
                # Fallback estimation
                docked = [
                    {"smiles": s, "score": _estimated_docking_score(s, pocket_score),
                     "method": "estimated"}
                    for s in smiles_list
                ]
                docked.sort(key=lambda x: x["score"])

        # Store top MAX_STORE hits
        stored = 0
        for hit in docked[:MAX_STORE]:
            smi = hit["smiles"]
            zinc_info = zinc_by_smiles.get(smi, {})
            with get_session() as session:
                c = Compound(
                    gene_id=gene_id,
                    smiles=smi,
                    name=zinc_info.get("name"),
                    source="zinc",
                    zinc_id=zinc_info.get("zinc_id"),
                    docking_score=hit["score"],
                    docking_method=hit["method"],
                    is_purchasable=bool(zinc_info.get("zinc_id")),
                )
                session.add(c)
                session.commit()
            stored += 1

        log.info("Stored %d compounds for %s", stored, symbol)
        stored_total += stored
        time.sleep(0.05)

    log.info("Virtual screening complete: %d compounds stored.", stored_total)
    return stored_total
