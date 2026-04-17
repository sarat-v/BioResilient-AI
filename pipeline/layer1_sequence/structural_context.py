"""Step 9b — Structural context annotation for convergent motif regions.

For each Tier1/Tier2 gene this module characterises the structural environment
of every convergently-evolved motif (DivergentMotif with convergent_aa_count
>= MIN_CONVERGENT_LINEAGES).  The analysis answers four questions per motif:

  1. Is this change at a functionally annotated site?
     → UniProt REST API feature table (active_site, binding_site, PTM).

  2. Is the specific convergent substitution functionally significant?
     → AlphaMissense pathogenicity score for the exact human→animal change.
     The full AlphaMissense index (built in Step 4b) is reused here; we make
     a targeted query for convergent positions instead of averaging all motif
     positions.

  3. How confident is the AlphaFold model at this residue?
     → pLDDT from the B-factor column of the AlphaFold PDB file.
     Disordered regions (pLDDT < 70) have low structural reliability and
     therefore low functional significance for drug targeting.

  4. How close is the change to the best druggable pocket?
     → Euclidean 3D distance (Å) from the Cα centroid of the convergent
     motif residues to the centroid of the top-ranked fpocket pocket.
     A distance ≤ 8 Å is flagged as pocket-adjacent.

Results are stored in ConvergentPositionAnnotation (one row per motif).
A per-gene structural_score [0,1] aggregates all rows; this score is picked
up by scoring._run_scoring_weighted() in Phase 2.

AlphaMissense reference: Cheng et al., Science 2023, 10.1126/science.adg7492
UniProt REST API: https://rest.uniprot.org/uniprotkb/{accession}.json
"""

import logging
import math
import time
import uuid
from pathlib import Path
from typing import Optional

import requests

from db.models import (
    CandidateScore,
    ConvergentPositionAnnotation,
    DivergentMotif,
    Gene,
    Ortholog,
)
from db.session import get_session
from pipeline.layer4_druggability.pockets import run_fpocket
from pipeline.layer4_druggability.structure import download_alphafold_structure
from pipeline.layer1_sequence.alphamissense import (
    _get_gene_accession_map,
    _am_path,
    build_am_index,
)

log = logging.getLogger(__name__)

# AlphaMissense classification thresholds (Cheng et al. 2023, Table S2)
AM_PATHOGENIC_THRESHOLD = 0.564
AM_BENIGN_THRESHOLD     = 0.340

# Minimum number of independent lineages sharing a convergent change
MIN_CONVERGENT_LINEAGES = 3

# Pocket adjacency: residue Cα centroid within this distance of pocket centroid
POCKET_ADJACENT_ANGSTROM = 8.0

# AlphaFold pLDDT below this threshold = structurally disordered
DISORDERED_PLDDT = 70.0

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
REQUEST_TIMEOUT = 20

# UniProt feature types that indicate functional importance
_FUNCTIONAL_FEATURE_TYPES = frozenset({
    "Active site",
    "Binding site",
    "Metal binding",
    "Site",
    "Modified residue",
    "Disulfide bond",
})


# ---------------------------------------------------------------------------
# 1. UniProt functional site annotation
# ---------------------------------------------------------------------------


def fetch_uniprot_features(accession: str) -> dict[int, dict]:
    """Fetch per-position functional site annotations from UniProt.

    Returns {position_1indexed: {'type': str, 'description': str}} for all
    positions covered by active/binding/metal/site/modified features.
    If multiple features overlap a position, the highest-priority one is kept
    (Active site > Binding site > others).
    """
    priority = {
        "Active site": 0,
        "Binding site": 1,
        "Metal binding": 2,
        "Site": 3,
        "Modified residue": 4,
        "Disulfide bond": 5,
    }
    try:
        r = requests.get(
            f"{UNIPROT_API}/{accession}.json",
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code != 200:
            log.warning("UniProt features %s: HTTP %d", accession, r.status_code)
            return {}
        data = r.json()
        features = data.get("features", [])
        result: dict[int, dict] = {}
        for feat in features:
            ft_type = feat.get("type", "")
            if ft_type not in _FUNCTIONAL_FEATURE_TYPES:
                continue
            location = feat.get("location", {})
            start = location.get("start", {}).get("value")
            end = location.get("end", {}).get("value")
            if start is None:
                continue
            end = end or start
            description = feat.get("description", "") or ft_type.lower()
            ft_norm = ft_type.lower().replace(" ", "_")
            prio = priority.get(ft_type, 99)
            for pos in range(start, end + 1):
                existing = result.get(pos)
                if existing is None or prio < priority.get(
                    existing["type"].replace("_", " ").title(), 99
                ):
                    result[pos] = {
                        "type": ft_norm,
                        "description": description,
                    }
        return result
    except Exception as exc:
        log.warning("UniProt features %s: %s", accession, exc)
        return {}


def _best_feature_for_range(
    feature_map: dict[int, dict], start: int, end: int
) -> Optional[dict]:
    """Return the highest-priority UniProt feature in a residue range."""
    priority = {
        "active_site": 0,
        "binding_site": 1,
        "metal_binding": 2,
        "site": 3,
        "modified_residue": 4,
        "disulfide_bond": 5,
    }
    best: Optional[dict] = None
    best_p = 999
    for pos in range(start, end + 1):
        feat = feature_map.get(pos)
        if feat:
            p = priority.get(feat["type"], 99)
            if p < best_p:
                best_p = p
                best = feat
    return best


# ---------------------------------------------------------------------------
# 2. AlphaMissense per-motif score
# ---------------------------------------------------------------------------


def _score_motif_am(
    human_seq: str,
    animal_seq: str,
    human_protein: str,
    am_protein_index: dict[int, dict[str, float]],
) -> tuple[Optional[float], Optional[str]]:
    """Compute mean AlphaMissense score for divergent positions in one motif.

    Returns (score, am_class) where score is mean AM pathogenicity [0–1] for
    the specific convergent (human→animal) substitutions in this motif, and
    am_class is 'likely_pathogenic', 'ambiguous', or 'likely_benign'.

    Reuses the same gap-corrected position mapping as _consequence_for_motif
    in alphamissense.py — see that module for the full rationale.
    """
    h_stripped = human_seq.replace("-", "")
    if not h_stripped:
        return None, None
    protein_start = human_protein.find(h_stripped)
    if protein_start < 0:
        return None, None

    scores: list[float] = []
    res_pos = protein_start  # 0-indexed
    for h_aa, a_aa in zip(human_seq, animal_seq):
        if h_aa == "-":
            continue
        if h_aa != a_aa and a_aa != "-":
            pos_1 = res_pos + 1  # 1-indexed for AM lookup
            aa_scores = am_protein_index.get(pos_1, {})
            if a_aa in aa_scores:
                scores.append(aa_scores[a_aa])
        res_pos += 1

    if not scores:
        return None, None
    mean_score = round(sum(scores) / len(scores), 4)
    if mean_score >= AM_PATHOGENIC_THRESHOLD:
        am_class = "likely_pathogenic"
    elif mean_score <= AM_BENIGN_THRESHOLD:
        am_class = "likely_benign"
    else:
        am_class = "ambiguous"
    return mean_score, am_class


def _map_motif_to_protein_positions(
    human_seq: str, human_protein: str
) -> Optional[tuple[int, int]]:
    """Map a motif's human_seq to 1-indexed (start, end) in the canonical protein.

    Returns None if the motif cannot be located in the human protein sequence.
    """
    h_stripped = human_seq.replace("-", "")
    if not h_stripped:
        return None
    idx = human_protein.find(h_stripped)
    if idx < 0:
        return None
    start_1 = idx + 1
    end_1 = idx + len(h_stripped)
    return start_1, end_1


# ---------------------------------------------------------------------------
# 3. AlphaFold structure: pLDDT and Cα coordinates
# ---------------------------------------------------------------------------


def _parse_pdb_ca(pdb_path: Path) -> dict[int, dict]:
    """Parse an AlphaFold PDB and return per-residue Cα info.

    Returns {res_num_1indexed: {'plddt': float, 'x': float, 'y': float, 'z': float}}.
    pLDDT is stored in the B-factor column (standard for AlphaFold v2+ models).
    """
    result: dict[int, dict] = {}
    try:
        for line in pdb_path.read_text(errors="replace").splitlines():
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            try:
                res_num = int(line[22:26].strip())
                bfactor = float(line[60:66].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                result[res_num] = {"plddt": bfactor, "x": x, "y": y, "z": z}
            except (ValueError, IndexError):
                pass
    except Exception as exc:
        log.debug("PDB Cα parse %s: %s", pdb_path, exc)
    return result


def _centroid(points: list[tuple[float, float, float]]) -> Optional[tuple[float, float, float]]:
    """Geometric centroid of a list of 3D points."""
    if not points:
        return None
    n = len(points)
    return (
        sum(p[0] for p in points) / n,
        sum(p[1] for p in points) / n,
        sum(p[2] for p in points) / n,
    )


def _distance(p1: tuple[float, float, float], p2: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


# ---------------------------------------------------------------------------
# 4. fpocket pocket centroid (from Cα coordinates of pocket residues)
# ---------------------------------------------------------------------------


def _pocket_centroid_from_ca(
    pdb_path: Path, ca_coords: dict[int, dict]
) -> Optional[tuple[float, float, float]]:
    """Run fpocket on pdb_path and return the Cα centroid of the top-pocket residues.

    Uses run_fpocket() which handles pLDDT filtering internally. The pocket
    centroid is computed from the Cα coordinates of residues in the top-ranked
    pocket (a physically meaningful reference point for distance calculations).

    Returns None if fpocket fails, finds no pockets, or no Cα coords are
    available for the pocket residues (e.g. the protein has no AlphaFold model).
    """
    result = run_fpocket(pdb_path)
    if result is None:
        return None
    pocket_residues: set[int] = result.get("top_pocket_residues", set())
    if not pocket_residues or not ca_coords:
        return None
    pts = [
        (ca_coords[r]["x"], ca_coords[r]["y"], ca_coords[r]["z"])
        for r in pocket_residues
        if r in ca_coords
    ]
    return _centroid(pts) if pts else None


# ---------------------------------------------------------------------------
# 5. Structural context classifier
# ---------------------------------------------------------------------------


def _classify_context(
    uniprot_feature: Optional[str],
    is_pocket_adjacent: bool,
    is_disordered: bool,
    pocket_distance: Optional[float],
) -> str:
    """Assign a structural context label to a convergent motif region.

    Priority order (most to least clinically relevant):
      active_site      — catalytic / active site residues
      binding_site     — ligand / cofactor / metal binding
      disordered       — pLDDT < 70; model unreliable
      pocket_adjacent  — within 8 Å of top druggable pocket centroid
      surface_near_pocket — within 15 Å; potentially accessible
      distal           — far from any known pocket
    """
    if uniprot_feature in ("active_site",):
        return "active_site"
    if uniprot_feature in ("binding_site", "metal_binding"):
        return "binding_site"
    if is_disordered:
        return "disordered"
    if is_pocket_adjacent:
        return "pocket_adjacent"
    if pocket_distance is not None and pocket_distance < 15.0:
        return "surface_near_pocket"
    return "distal"


# ---------------------------------------------------------------------------
# 6. Per-gene structural score
# ---------------------------------------------------------------------------


def compute_gene_structural_score(gene_id: str) -> float:
    """Aggregate ConvergentPositionAnnotation rows into a [0,1] structural score.

    Weights:
      0.40 — fraction of convergent motifs in functionally important context
              (active_site, binding_site, or pocket_adjacent)
      0.30 — mean AlphaMissense score (higher = changes at functionally critical positions)
      0.20 — fraction with UniProt functional annotation (non-null uniprot_feature)
      0.10 — mean pLDDT / 100 (penalises disordered convergence; max confidence = 1.0)

    Returns 0.0 if no annotation rows exist for this gene (step 9b not yet run
    or gene had no convergent motifs above the lineage threshold).
    """
    with get_session() as session:
        rows = (
            session.query(ConvergentPositionAnnotation)
            .filter_by(gene_id=gene_id)
            .all()
        )
    if not rows:
        return 0.0

    n = len(rows)
    functional_contexts = {"active_site", "binding_site", "pocket_adjacent"}
    functional_frac = sum(
        1 for r in rows if r.structural_context in functional_contexts
    ) / n

    am_scores = [r.am_score for r in rows if r.am_score is not None]
    mean_am = sum(am_scores) / len(am_scores) if am_scores else 0.0

    feature_frac = sum(1 for r in rows if r.uniprot_feature is not None) / n

    plddt_vals = [r.plddt_mean for r in rows if r.plddt_mean is not None]
    plddt_norm = (sum(plddt_vals) / len(plddt_vals)) / 100.0 if plddt_vals else 0.5

    score = (
        0.40 * functional_frac
        + 0.30 * mean_am
        + 0.20 * feature_frac
        + 0.10 * plddt_norm
    )
    return round(min(score, 1.0), 4)


# ---------------------------------------------------------------------------
# 7. Main annotation entry point
# ---------------------------------------------------------------------------


def annotate_convergent_structural_context(gene_ids: Optional[list[str]] = None) -> int:
    """Annotate convergent motif regions with structural context for Step 9b.

    For each gene in gene_ids (default: all Tier1/Tier2 CandidateScore genes):
      1. Collect DivergentMotif rows with convergent_aa_count >= MIN_CONVERGENT_LINEAGES.
      2. For each motif, map to protein coordinates, fetch UniProt features,
         AlphaMissense score, pLDDT, and pocket proximity.
      3. Upsert ConvergentPositionAnnotation row.
      4. Update CandidateScore.structural_score for the gene.

    Returns number of genes with at least one annotated motif.
    """
    # ---- resolve gene set ------------------------------------------------
    if gene_ids is None:
        with get_session() as session:
            gene_ids = [
                r.gene_id
                for r in session.query(CandidateScore)
                .filter(CandidateScore.tier.in_(["Tier1", "Tier2", "Validated"]))
                .all()
            ]
    if not gene_ids:
        log.info("Step 9b: no Tier1/Tier2 genes found — nothing to annotate.")
        return 0

    log.info("Step 9b: annotating structural context for %d genes.", len(gene_ids))

    # ---- build AlphaMissense in-memory index for target proteins ------------
    am_tsv = _am_path()
    if am_tsv.exists():
        gene_acc_map = _get_gene_accession_map()
        target_accessions = set(gene_acc_map.values()) if gene_acc_map else None
        am_index = build_am_index(am_tsv, target_accessions=target_accessions)
        log.info("Step 9b: AlphaMissense index loaded (%d proteins).", len(am_index))
    else:
        am_index = {}
        log.warning("Step 9b: AlphaMissense TSV not found — AM scores will be NULL.")

    # ---- collect gene metadata in one DB query ----------------------------
    with get_session() as session:
        genes = session.query(Gene).filter(Gene.id.in_(gene_ids)).all()
        gene_map = {g.id: g for g in genes}

        # Pre-fetch human protein sequences (from Ortholog where species_id='human')
        human_seqs = {
            row.gene_id: row.protein_seq.replace("-", "").replace("*", "")
            for row in (
                session.query(Ortholog.gene_id, Ortholog.protein_seq)
                .filter(
                    Ortholog.species_id == "human",
                    Ortholog.gene_id.in_(gene_ids),
                    Ortholog.protein_seq.isnot(None),
                )
                .all()
            )
            if row.protein_seq
        }

    # ---- pre-resolve UniProt accessions -----------------------------------
    gene_acc_map = _get_gene_accession_map()

    annotated_genes = 0

    for gene_id in gene_ids:
        gene = gene_map.get(gene_id)
        if not gene:
            continue

        accession = gene_acc_map.get(gene_id)
        human_protein_seq = human_seqs.get(gene_id, "")
        symbol = (gene.gene_symbol or "").split("_")[0].split("|")[-1]

        log.info("Step 9b: processing %s (accession=%s)", symbol, accession)

        # ---- fetch UniProt features once per gene -------------------------
        uniprot_feature_map: dict[int, dict] = {}
        if accession:
            uniprot_feature_map = fetch_uniprot_features(accession)
            time.sleep(0.15)   # polite rate-limiting for UniProt REST

        # ---- resolve AlphaFold structure + fpocket -----------------------
        pdb_path: Optional[Path] = None
        pdb_stem: Optional[str] = None
        ca_coords: dict[int, dict] = {}
        pocket_centroid: Optional[tuple[float, float, float]] = None

        if gene.human_protein:
            pdb_path = download_alphafold_structure(gene.human_protein)
        if pdb_path and pdb_path.exists():
            pdb_stem = pdb_path.stem
            ca_coords = _parse_pdb_ca(pdb_path)
            # Compute pocket centroid from Cα coords of top-pocket residues.
            # run_fpocket handles pLDDT filtering internally; the Cα centroid is
            # more meaningful than alpha-sphere centroids for distance calculations.
            pocket_centroid = _pocket_centroid_from_ca(pdb_path, ca_coords)
        else:
            log.info("Step 9b: %s — no AlphaFold structure available.", symbol)

        # ---- get AlphaMissense index for this protein --------------------
        am_protein_index: dict[int, dict[str, float]] = {}
        if accession and am_index:
            am_protein_index = am_index.get(accession.upper(), {})

        # ---- fetch convergent motifs for this gene -----------------------
        with get_session() as session:
            convergent_motifs = (
                session.query(DivergentMotif)
                .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
                .filter(
                    Ortholog.gene_id == gene_id,
                    DivergentMotif.convergent_aa_count >= MIN_CONVERGENT_LINEAGES,
                )
                .all()
            )

        if not convergent_motifs:
            log.info("Step 9b: %s — no convergent motifs (count >= %d).",
                     symbol, MIN_CONVERGENT_LINEAGES)
            continue

        log.info("Step 9b: %s — %d convergent motifs to annotate.",
                 symbol, len(convergent_motifs))

        gene_got_annotation = False
        motif_annotations: list[dict] = []

        for motif in convergent_motifs:
            # ---- map motif to protein positions --------------------------
            pos_range = None
            if human_protein_seq and motif.human_seq:
                pos_range = _map_motif_to_protein_positions(
                    motif.human_seq, human_protein_seq
                )

            prot_start = pos_range[0] if pos_range else None
            prot_end   = pos_range[1] if pos_range else None

            # ---- AlphaMissense -------------------------------------------
            am_score: Optional[float] = None
            am_class: Optional[str] = None
            if human_protein_seq and motif.human_seq and motif.animal_seq and am_protein_index:
                am_score, am_class = _score_motif_am(
                    motif.human_seq,
                    motif.animal_seq,
                    human_protein_seq,
                    am_protein_index,
                )
            # Fallback: use pre-computed consequence_score from step 4b
            if am_score is None and motif.consequence_score is not None:
                am_score = motif.consequence_score
                if am_score >= AM_PATHOGENIC_THRESHOLD:
                    am_class = "likely_pathogenic"
                elif am_score <= AM_BENIGN_THRESHOLD:
                    am_class = "likely_benign"
                else:
                    am_class = "ambiguous"

            # ---- UniProt feature for this range --------------------------
            best_feat: Optional[dict] = None
            if prot_start and uniprot_feature_map:
                best_feat = _best_feature_for_range(
                    uniprot_feature_map, prot_start, prot_end or prot_start
                )

            # ---- pLDDT for the motif residues ----------------------------
            plddt_mean: Optional[float] = None
            is_disordered: Optional[bool] = None
            if prot_start and ca_coords:
                residues_in_range = [
                    ca_coords[r]["plddt"]
                    for r in range(prot_start, (prot_end or prot_start) + 1)
                    if r in ca_coords
                ]
                if residues_in_range:
                    plddt_mean = round(
                        sum(residues_in_range) / len(residues_in_range), 2
                    )
                    is_disordered = plddt_mean < DISORDERED_PLDDT

            # ---- pocket proximity (3D) -----------------------------------
            pocket_dist: Optional[float] = None
            is_pocket_adj: Optional[bool] = None
            if prot_start and ca_coords and pocket_centroid:
                motif_ca_pts = [
                    (ca_coords[r]["x"], ca_coords[r]["y"], ca_coords[r]["z"])
                    for r in range(prot_start, (prot_end or prot_start) + 1)
                    if r in ca_coords
                ]
                if motif_ca_pts:
                    motif_cent = _centroid(motif_ca_pts)
                    if motif_cent:
                        pocket_dist = round(_distance(motif_cent, pocket_centroid), 2)
                        is_pocket_adj = pocket_dist <= POCKET_ADJACENT_ANGSTROM

            # ---- structural context classifier ---------------------------
            # Treat None pLDDT as 'unknown' rather than 'ordered' (bool(None)==False)
            is_disordered_known: Optional[bool] = is_disordered  # may be None
            structural_ctx = _classify_context(
                uniprot_feature=best_feat["type"] if best_feat else None,
                is_pocket_adjacent=bool(is_pocket_adj),
                is_disordered=bool(is_disordered_known) if is_disordered_known is not None else False,
                pocket_distance=pocket_dist,
            )

            log.info(
                "  %s motif@%s-%s: AM=%.3f(%s) feat=%s plddt=%.1f dist=%s ctx=%s",
                symbol,
                prot_start, prot_end,
                am_score or 0, am_class,
                best_feat["type"] if best_feat else "—",
                plddt_mean or 0,
                f"{pocket_dist:.1f}Å" if pocket_dist is not None else "—",
                structural_ctx,
            )

            # ---- collect annotation data for batched DB write ------------
            motif_annotations.append({
                "motif_id": str(motif.id),
                "accession": accession,
                "protein_start_pos": prot_start,
                "protein_end_pos": prot_end,
                "convergent_lineage_count": motif.convergent_aa_count,
                "am_score": am_score,
                "am_class": am_class,
                "uniprot_feature": best_feat["type"] if best_feat else None,
                "uniprot_feature_description": best_feat["description"] if best_feat else None,
                "plddt_mean": plddt_mean,
                "is_disordered": is_disordered_known,
                "pocket_distance_angstrom": pocket_dist,
                "is_pocket_adjacent": is_pocket_adj,
                "structural_context": structural_ctx,
            })
            gene_got_annotation = True

        # ---- batched DB write for all motifs of this gene ---------------
        if motif_annotations:
            with get_session() as session:
                for ann in motif_annotations:
                    row = (
                        session.query(ConvergentPositionAnnotation)
                        .filter_by(gene_id=gene_id, divergent_motif_id=ann["motif_id"])
                        .first()
                    )
                    if row is None:
                        row = ConvergentPositionAnnotation(
                            id=str(uuid.uuid4()),
                            gene_id=gene_id,
                            divergent_motif_id=ann["motif_id"],
                        )
                        session.add(row)
                    row.uniprot_accession = ann["accession"]
                    row.protein_start_pos = ann["protein_start_pos"]
                    row.protein_end_pos   = ann["protein_end_pos"]
                    row.convergent_lineage_count = ann["convergent_lineage_count"]
                    row.am_score          = ann["am_score"]
                    row.am_class          = ann["am_class"]
                    row.uniprot_feature   = ann["uniprot_feature"]
                    row.uniprot_feature_description = ann["uniprot_feature_description"]
                    row.plddt_mean        = ann["plddt_mean"]
                    row.is_disordered     = ann["is_disordered"]
                    row.pocket_distance_angstrom = ann["pocket_distance_angstrom"]
                    row.is_pocket_adjacent = ann["is_pocket_adjacent"]
                    row.structural_context = ann["structural_context"]
                # Single commit for entire gene's motifs
                session.commit()

        # ---- update CandidateScore.structural_score ----------------------
        if gene_got_annotation:
            s_score = compute_gene_structural_score(gene_id)
            with get_session() as session:
                for cs in session.query(CandidateScore).filter_by(gene_id=gene_id).all():
                    cs.structural_score = s_score
                session.commit()
            log.info("Step 9b: %s structural_score = %.4f", symbol, s_score)
            annotated_genes += 1

    log.info(
        "Step 9b complete: %d / %d genes annotated with structural context.",
        annotated_genes, len(gene_ids),
    )
    return annotated_genes
