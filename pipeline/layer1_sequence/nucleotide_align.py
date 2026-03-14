"""Step 3c (part 2) — Nucleotide alignment, conservation scoring, and regulatory divergence.

For each gene, takes the NucleotideRegion sequences (written by nucleotide_scan.py) and:

  1. Aligns each resilient species' sequence against the human reference using
     minimap2 -x asm5 (assembly-to-assembly). Falls back to LASTZ, then BLASTN.

  2. Computes per-gene per-region conservation statistics:
       conservation_score  = (pct_identity × alignment_length) / region_length
       percent_identity    = mean across resilient species
       alignment_length    = mean
       gap_fraction        = mean

  3. Detects regulatory divergence in promoter regions:
       regulatory_divergence_count:  position mutated in ≥3 resilient species,
                                     absent in human AND absent in all controls
       regulatory_convergence_count: same mutation independently present in ≥2
                                     distinct evolutionary lineage clusters

  4. Upserts NucleotideScore rows.

Entry point: run_nucleotide_alignment(gene_ids=None) -> int
"""

import logging
import os
import re
import shutil
import subprocess
import tempfile
import uuid
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from db.models import Gene, NucleotideRegion, NucleotideScore, Species
from db.session import get_session

log = logging.getLogger(__name__)

# Lineage clusters used to determine regulatory convergence.
# Each species' lineage_group (stored in Species.lineage_group) is checked
# against this set to assign a cluster label.
_KNOWN_LINEAGE_CLUSTERS = {
    "Rodents", "Cetaceans", "Proboscideans", "Bats",
    "Sharks", "Turtles", "Cnidarians", "Molluscs",
    "Primates",
}

# Minimum number of resilient species that must share a mutation for it to count
# as regulatory divergence (absent in human + controls).
_MIN_RESILIENT_SHARING = 3

# Minimum number of distinct lineage clusters that must share the same mutation
# for it to count as regulatory convergence.
_MIN_CONVERGENT_CLUSTERS = 2


# ─────────────────────────────────────────────────────────────────────────────
# Tool detection
# ─────────────────────────────────────────────────────────────────────────────

def _detect_aligner() -> str:
    """Return the first available nucleotide aligner: blastn → minimap2 → lastz.

    blastn is preferred for cross-species comparisons since minimap2 presets
    are optimised for same-species or high-identity assemblies.
    """
    for tool in ("blastn", "minimap2", "lastz"):
        if shutil.which(tool):
            return tool
    log.warning("No nucleotide aligner found (blastn/minimap2/lastz). Alignment step will be skipped.")
    return ""


# ─────────────────────────────────────────────────────────────────────────────
# Alignment functions
# ─────────────────────────────────────────────────────────────────────────────

def _write_fasta(path: Path, seqid: str, seq: str) -> None:
    path.write_text(f">{seqid}\n{seq}\n")


def _run_minimap2(query_seq: str, target_seq: str, query_id: str = "query", target_id: str = "target", threads: int = 1) -> Optional[dict]:
    """Align query → target with minimap2. Uses asm5 for sequences ≥10 kb,
    falls back to sr (short-read) preset for shorter sequences."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        query_fa  = tmp / "query.fa"
        target_fa = tmp / "target.fa"
        _write_fasta(query_fa,  query_id,  query_seq)
        _write_fasta(target_fa, target_id, target_seq)
        # asm5 requires long sequences; use sr preset for short regions
        preset = "asm5" if min(len(query_seq), len(target_seq)) >= 10_000 else "sr"
        try:
            result = subprocess.run(
                ["minimap2", "-x", preset, "-t", str(threads), "--cs", "-c", str(target_fa), str(query_fa)],
                capture_output=True, text=True, timeout=120,
            )
            parsed = _parse_paf(result.stdout, len(query_seq))
            if parsed:
                return parsed
            # If sr preset also failed, try with no preset (general purpose)
            if preset == "sr":
                result2 = subprocess.run(
                    ["minimap2", "-t", str(threads), "--cs", "-c", str(target_fa), str(query_fa)],
                    capture_output=True, text=True, timeout=120,
                )
                return _parse_paf(result2.stdout, len(query_seq))
            return None
        except Exception as exc:
            log.debug("minimap2 failed: %s", exc)
            return None


def _parse_paf(paf_output: str, query_len: int) -> Optional[dict]:
    """Parse PAF format output from minimap2 and return alignment stats."""
    best: Optional[dict] = None
    best_aln = 0
    for line in paf_output.strip().splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 12:
            continue
        try:
            aln_len   = int(cols[10])
            n_matches = int(cols[9])
            n_blocks  = int(cols[11])
        except (ValueError, IndexError):
            continue
        if aln_len < 1:
            continue
        pct_id = round(n_matches / aln_len * 100.0, 2)
        # Count gaps from CIGAR-like cs tag if present
        gap_count = 0
        for col in cols[12:]:
            if col.startswith("cg:Z:"):
                cigar = col[5:]
                for m in re.finditer(r"(\d+)([ID])", cigar):
                    gap_count += int(m.group(1))
        gap_frac = round(gap_count / aln_len, 4) if aln_len > 0 else 0.0
        if aln_len > best_aln:
            best_aln = aln_len
            best = {
                "pct_identity":    pct_id,
                "alignment_length": aln_len,
                "gap_fraction":     gap_frac,
            }
    return best


def _run_blastn(query_seq: str, target_seq: str, threads: int = 1) -> Optional[dict]:
    """Fallback alignment using BLASTN when minimap2/lastz are unavailable."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        query_fa  = tmp / "query.fa"
        target_fa = tmp / "target.fa"
        db_path   = tmp / "blastdb"
        _write_fasta(query_fa,  "query",  query_seq)
        _write_fasta(target_fa, "target", target_seq)
        try:
            subprocess.run(
                ["makeblastdb", "-in", str(target_fa), "-dbtype", "nucl", "-out", str(db_path)],
                capture_output=True, timeout=60, check=True,
            )
            res = subprocess.run(
                ["blastn", "-query", str(query_fa), "-db", str(db_path),
                 "-outfmt", "6 pident length gaps",
                 "-num_alignments", "1", "-dust", "no",
                 "-num_threads", str(threads)],
                capture_output=True, text=True, timeout=120,
            )
            for line in res.stdout.strip().splitlines():
                cols = line.split("\t")
                if len(cols) >= 3:
                    pct_id = float(cols[0])
                    aln    = int(cols[1])
                    gaps   = int(cols[2])
                    return {
                        "pct_identity":    round(pct_id, 2),
                        "alignment_length": aln,
                        "gap_fraction":     round(gaps / aln, 4) if aln > 0 else 0.0,
                    }
        except Exception as exc:
            log.debug("BLASTN fallback failed: %s", exc)
        return None


def _run_lastz(query_seq: str, target_seq: str) -> Optional[dict]:
    """Fallback alignment using LASTZ."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        query_fa  = tmp / "query.fa"
        target_fa = tmp / "target.fa"
        _write_fasta(query_fa,  "query",  query_seq)
        _write_fasta(target_fa, "target", target_seq)
        try:
            res = subprocess.run(
                ["lastz", f"{target_fa}[multiple]", str(query_fa),
                 "--format=general:identity,length,gaps"],
                capture_output=True, text=True, timeout=120,
            )
            for line in res.stdout.strip().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.split("\t")
                if len(cols) >= 3:
                    try:
                        pct_id = float(cols[0].replace("%", ""))
                        aln    = int(cols[1])
                        gaps   = int(cols[2])
                        return {
                            "pct_identity":    round(pct_id, 2),
                            "alignment_length": aln,
                            "gap_fraction":     round(gaps / aln, 4) if aln > 0 else 0.0,
                        }
                    except ValueError:
                        continue
        except Exception as exc:
            log.debug("LASTZ fallback failed: %s", exc)
        return None


def align_sequences(query_seq: str, target_seq: str, aligner: str) -> Optional[dict]:
    """Run alignment with the specified tool and return stats dict.
    Always falls back to blastn if primary aligner returns nothing."""
    if not query_seq or not target_seq:
        return None
    result = None
    if aligner == "minimap2":
        result = _run_minimap2(query_seq, target_seq)
        if result is None:
            result = _run_blastn(query_seq, target_seq)
    elif aligner == "lastz":
        result = _run_lastz(query_seq, target_seq)
        if result is None:
            result = _run_blastn(query_seq, target_seq)
    elif aligner == "blastn":
        result = _run_blastn(query_seq, target_seq)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Conservation score
# ─────────────────────────────────────────────────────────────────────────────

def compute_conservation_score(pct_identity: float, alignment_length: int, region_length: int) -> float:
    """Normalised conservation score: (pct_identity/100 × alignment_length) / region_length."""
    if region_length <= 0:
        return 0.0
    return round(min((pct_identity / 100.0) * alignment_length / region_length, 1.0), 4)


# ─────────────────────────────────────────────────────────────────────────────
# Regulatory divergence detection
# ─────────────────────────────────────────────────────────────────────────────

def _extract_variant_positions(query_seq: str, target_seq: str) -> set[int]:
    """Return 0-indexed positions where query differs from target using a simple
    column-by-column comparison of equal-length sequences.

    When sequences have different lengths (due to indels), we do a rough
    positional comparison on the shorter sequence length.
    """
    positions: set[int] = set()
    min_len = min(len(query_seq), len(target_seq))
    for i in range(min_len):
        if query_seq[i].upper() != target_seq[i].upper():
            positions.add(i)
    return positions


def compute_regulatory_divergence(
    gene_id: str,
    region_type: str,
    regions_by_species: dict[str, NucleotideRegion],
    species_meta: dict[str, Species],
) -> tuple[int, int]:
    """Compute regulatory_divergence_count and regulatory_convergence_count.

    regulatory_divergence_count:
      Number of base positions mutated (vs. human) in ≥3 resilient species
      AND absent in human (by definition) AND absent in all control species.

    regulatory_convergence_count:
      Of those divergent positions, count how many appear independently in
      ≥2 distinct evolutionary lineage clusters.

    Returns (divergence_count, convergence_count).
    """
    human_region = regions_by_species.get("human")
    if human_region is None or not human_region.sequence:
        return 0, 0

    human_seq = human_region.sequence

    # Separate resilient and control species
    resilient_seqs: dict[str, str]   = {}   # species_id → sequence
    control_seqs: dict[str, str]     = {}
    lineage_map: dict[str, str]      = {}   # species_id → lineage_group

    for sid, region in regions_by_species.items():
        if sid == "human" or not region.sequence:
            continue
        sp = species_meta.get(sid)
        if sp is None:
            continue
        if sp.is_control:
            control_seqs[sid] = region.sequence
        else:
            resilient_seqs[sid] = region.sequence
            lineage_map[sid] = sp.lineage_group or "Unknown"

    if not resilient_seqs:
        return 0, 0

    # Find positions mutated vs. human for each resilient species
    per_species_variants: dict[str, set[int]] = {}
    for sid, seq in resilient_seqs.items():
        per_species_variants[sid] = _extract_variant_positions(seq, human_seq)

    # Positions mutated in control species (to exclude them)
    control_variant_positions: set[int] = set()
    for seq in control_seqs.values():
        control_variant_positions |= _extract_variant_positions(seq, human_seq)

    # Regulatory divergence: mutated in ≥3 resilient, absent in controls
    position_sharing: dict[int, list[str]] = defaultdict(list)
    for sid, positions in per_species_variants.items():
        for pos in positions:
            position_sharing[pos].append(sid)

    divergent_positions: set[int] = set()
    for pos, sharing_species in position_sharing.items():
        if len(sharing_species) >= _MIN_RESILIENT_SHARING and pos not in control_variant_positions:
            divergent_positions.add(pos)

    divergence_count = len(divergent_positions)

    # Regulatory convergence: same divergent position in ≥2 distinct lineage clusters
    convergence_count = 0
    for pos in divergent_positions:
        sharing_species = position_sharing[pos]
        clusters_represented = {lineage_map.get(sid, "Unknown") for sid in sharing_species}
        if len(clusters_represented) >= _MIN_CONVERGENT_CLUSTERS:
            convergence_count += 1

    return divergence_count, convergence_count


# ─────────────────────────────────────────────────────────────────────────────
# DB upsert
# ─────────────────────────────────────────────────────────────────────────────

def _upsert_nucleotide_score(session, gene_id: str, region_type: str, stats: dict) -> None:
    existing = (
        session.query(NucleotideScore)
        .filter_by(gene_id=gene_id, region_type=region_type)
        .first()
    )
    if existing is None:
        existing = NucleotideScore(id=str(uuid.uuid4()), gene_id=gene_id, region_type=region_type)
        session.add(existing)
    for k, v in stats.items():
        setattr(existing, k, v)


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def _score_gene(args: tuple) -> tuple[str, dict[str, dict] | None]:
    """Worker: align all region_types for a single gene and return scores.

    Runs in a subprocess worker — DB writes are done in the main process.
    Returns (gene_id, {region_type: stats_dict}) or (gene_id, None) on skip.
    """
    gene_id, gene_symbol, gene_regions_data, species_meta_data, aligner = args

    # gene_regions_data: [{region_type, species_id, sequence}, ...]
    by_type: dict[str, dict] = defaultdict(dict)
    for row in gene_regions_data:
        by_type[row["region_type"]][row["species_id"]] = row["sequence"]

    results: dict[str, dict] = {}
    for region_type, regions_by_species in by_type.items():
        human_seq = regions_by_species.get("human")
        if not human_seq:
            continue

        region_len = len(human_seq)
        pct_ids, aln_lens, gap_fracs = [], [], []

        for sid, seq in regions_by_species.items():
            if sid == "human" or not seq:
                continue
            sp = species_meta_data.get(sid)
            if sp is None or sp.get("is_control"):
                continue
            result = align_sequences(seq, human_seq, aligner)
            if result:
                pct_ids.append(result["pct_identity"])
                aln_lens.append(result["alignment_length"])
                gap_fracs.append(result["gap_fraction"])

        if not pct_ids:
            log.debug("  No alignment results for gene %s region %s", gene_symbol, region_type)
            continue

        mean_pct_id = round(sum(pct_ids) / len(pct_ids), 2)
        mean_aln    = round(sum(aln_lens) / len(aln_lens))
        mean_gap    = round(sum(gap_fracs) / len(gap_fracs), 4)
        cons_score  = compute_conservation_score(mean_pct_id, mean_aln, region_len)

        results[region_type] = {
            "conservation_score":  cons_score,
            "percent_identity":    mean_pct_id,
            "alignment_length":    mean_aln,
            "gap_fraction":        mean_gap,
            # Regulatory divergence requires full Species ORM objects — computed in main process
            "regulatory_divergence_count":  0,
            "regulatory_convergence_count": 0,
        }

    return gene_id, results if results else None


def run_nucleotide_alignment(gene_ids: list[str] | None = None) -> int:
    """Align nucleotide regions and compute conservation + regulatory divergence scores.

    Reads NucleotideRegion rows, aligns each species against human in parallel
    (one worker per gene), aggregates conservation stats and regulatory
    divergence/convergence counts, then upserts NucleotideScore rows.

    Returns count of genes scored.
    """
    aligner = _detect_aligner()
    if not aligner:
        log.warning("Nucleotide alignment skipped — no aligner available. Install minimap2.")
        return 0

    n_cpu = os.cpu_count() or 8
    n_workers = max(1, n_cpu - 2)
    log.info("Nucleotide alignment using: %s  |  %d parallel workers", aligner, n_workers)

    with get_session() as session:
        species_meta: dict[str, Species] = {sp.id: sp for sp in session.query(Species).all()}
        # Serialisable form for workers (no ORM objects across process boundaries)
        species_meta_data: dict[str, dict] = {
            sid: {"is_control": sp.is_control, "lineage_group": sp.lineage_group}
            for sid, sp in species_meta.items()
        }

        q = session.query(Gene)
        if gene_ids:
            q = q.filter(Gene.id.in_(gene_ids))
        genes = q.all()
        log.info("Aligning nucleotide regions for %d genes...", len(genes))

        # Pre-fetch all NucleotideRegion rows in a single query and group by gene_id
        all_regions = session.query(NucleotideRegion)
        if gene_ids:
            all_regions = all_regions.filter(NucleotideRegion.gene_id.in_(gene_ids))
        regions_by_gene: dict[str, list[dict]] = defaultdict(list)
        for nr in all_regions:
            if nr.sequence:
                regions_by_gene[nr.gene_id].append({
                    "region_type": nr.region_type,
                    "species_id":  nr.species_id,
                    "sequence":    nr.sequence,
                })

        # Prepare worker args
        work_items = [
            (g.id, g.gene_symbol, regions_by_gene.get(g.id, []), species_meta_data, aligner)
            for g in genes
            if regions_by_gene.get(g.id)
        ]
        total = len(work_items)
        genes_scored = 0
        done = 0

        # Run alignment in parallel; collect results and write to DB in main process
        # (avoids SQLAlchemy session sharing across processes)
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {pool.submit(_score_gene, item): item[0] for item in work_items}
            for future in as_completed(futures):
                gene_id, region_scores = future.result()
                done += 1

                if region_scores:
                    # Compute regulatory divergence (needs ORM objects) in main process
                    gene_regions = regions_by_gene.get(gene_id, [])
                    by_type_full: dict[str, dict[str, NucleotideRegion]] = defaultdict(dict)
                    for row in gene_regions:
                        # We need real NucleotideRegion-like objects; reconstruct minimally
                        class _R:
                            pass
                        r = _R()
                        r.sequence = row["sequence"]
                        r.species_id = row["species_id"]
                        by_type_full[row["region_type"]][row["species_id"]] = r

                    for region_type, stats in region_scores.items():
                        if region_type == "promoter":
                            div_count, conv_count = compute_regulatory_divergence(
                                gene_id, region_type, by_type_full.get(region_type, {}), species_meta
                            )
                            stats["regulatory_divergence_count"]  = div_count
                            stats["regulatory_convergence_count"] = conv_count

                        _upsert_nucleotide_score(session, gene_id, region_type, stats)
                        genes_scored += 1

                if done % 500 == 0 or done == total:
                    log.info("  Aligned %d / %d genes (%.0f%%).", done, total, 100 * done / total)
                    session.commit()

        session.commit()

    log.info("Nucleotide alignment complete: %d genes scored", genes_scored)
    return genes_scored
