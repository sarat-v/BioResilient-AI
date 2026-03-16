"""Step 4b-ii — AlphaMissense functional consequence scoring for divergent motifs.

AlphaMissense (Cheng et al., 2023) provides precomputed pathogenicity scores
for all ~216M possible human missense variants. A score near 1 = likely
pathogenic; near 0 = likely benign.

For each DivergentMotif, this module:
  1. Identifies divergent positions (animal_seq[i] != human_seq[i], no gaps).
  2. Looks up the AlphaMissense score for substituting the human residue to the
     animal residue at each divergent protein position.
  3. Computes consequence_score = mean AlphaMissense score across all divergent
     positions in the motif window.

A high consequence_score means the divergence is at positions where substitutions
tend to be functionally significant in humans — much more informative than raw
divergence % alone.

Data source:
  AlphaMissense TSV (~1.6 GB compressed) from:
  https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
  Cached locally at {storage_root}/alphamissense/AlphaMissense_hg38.tsv.gz
  after first download.

The TSV format:
  #CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id  protein_variant  am_pathogenicity  am_class
  (protein_variant format: e.g. "A123V" — position is 1-indexed in the canonical transcript)
"""

import gzip
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Optional

import requests

from db.models import DivergentMotif, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_storage_root

log = logging.getLogger(__name__)

AM_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
AM_CHUNK = 1024 * 1024 * 8   # 8 MB download chunks


def _am_path() -> Path:
    root = Path(get_storage_root()) / "alphamissense"
    root.mkdir(parents=True, exist_ok=True)
    return root / "AlphaMissense_hg38.tsv.gz"


def _am_index_path() -> Path:
    return _am_path().parent / "am_index.db"


def download_alphamissense(force: bool = False) -> Path:
    """Download AlphaMissense TSV if not already cached.

    Returns path to the gzipped TSV file.
    Skips download if file already exists and force=False.
    """
    path = _am_path()
    if path.exists() and not force:
        log.info("AlphaMissense TSV already cached at %s", path)
        return path

    log.info("Downloading AlphaMissense scores (~1.6 GB)...")
    try:
        with requests.get(AM_URL, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(path, "wb") as f:
                downloaded = 0
                for chunk in r.iter_content(chunk_size=AM_CHUNK):
                    f.write(chunk)
                    downloaded += len(chunk)
                    if downloaded % (100 * 1024 * 1024) == 0:
                        log.info("  Downloaded %.0f MB...", downloaded / 1e6)
    except Exception as exc:
        log.warning("AlphaMissense download failed: %s — skipping consequence scoring.", exc)
        if path.exists():
            path.unlink()
        return path

    log.info("AlphaMissense download complete: %s", path)
    return path


def _get_target_accessions() -> set[str]:
    """Return the set of UniProt accessions present in our Gene table.

    Used to filter the AlphaMissense TSV during index build so we only keep
    entries for proteins we actually have motifs for, instead of loading all
    ~20 000 human proteins into memory.
    """
    from db.models import Gene
    from db.session import get_session

    accessions: set[str] = set()
    with get_session() as session:
        for gene in session.query(Gene).all():
            acc = gene.human_protein
            if not acc:
                continue
            if "|" in acc:
                acc = acc.split("|")[-1]
            acc = acc.split("-")[0].strip()  # strip isoform suffix
            if acc:
                accessions.add(acc.upper())
    return accessions


def build_am_index(
    tsv_gz_path: Path,
    target_accessions: Optional[set[str]] = None,
) -> dict[str, dict[int, dict[str, float]]]:
    """Build an in-memory index: {uniprot_id: {position: {alt_aa: score}}}.

    Args:
        tsv_gz_path: Path to AlphaMissense_hg38.tsv.gz.
        target_accessions: Optional set of UniProt IDs to keep. When provided,
            only variants for those proteins are loaded — reducing RAM from ~3 GB
            to typically <100 MB and cutting parse time from ~5 min to ~30 s.
            If None, all ~216 M variants are loaded (legacy behaviour).

    Returns empty dict if file does not exist.
    """
    if not tsv_gz_path.exists():
        return {}

    if target_accessions is None:
        target_accessions = _get_target_accessions()

    n_target = len(target_accessions)
    log.info(
        "Building AlphaMissense index for %d target proteins from %s...",
        n_target, tsv_gz_path,
    )
    index: dict[str, dict[int, dict[str, float]]] = defaultdict(lambda: defaultdict(dict))
    n_loaded = 0
    n_skipped = 0

    try:
        with gzip.open(tsv_gz_path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                # Columns: CHROM POS REF ALT genome uniprot_id transcript_id protein_variant am_pathogenicity am_class
                uniprot_id = parts[5].upper()

                # Skip proteins we don't need — this is the key optimisation.
                # Cuts the number of parsed lines from ~216M to a few million.
                if target_accessions and uniprot_id not in target_accessions:
                    n_skipped += 1
                    continue

                protein_variant = parts[7]    # e.g. "A123V"
                try:
                    am_score = float(parts[8])
                except ValueError:
                    continue

                m = re.match(r"^([A-Z])(\d+)([A-Z])$", protein_variant)
                if not m:
                    continue
                pos = int(m.group(2))
                alt_aa = m.group(3)

                index[uniprot_id][pos][alt_aa] = am_score
                n_loaded += 1

    except Exception as exc:
        log.warning("AlphaMissense index build failed: %s", exc)
        return {}

    log.info(
        "AlphaMissense index built: %d variants across %d proteins "
        "(%d variants skipped as off-target).",
        n_loaded, len(index), n_skipped,
    )
    return dict(index)


def _consequence_for_motif(
    human_seq: str,
    animal_seq: str,
    motif_start_pos: int,
    am_index: dict[int, dict[str, float]],
) -> Optional[float]:
    """Compute mean AlphaMissense score for all divergent positions in a motif.

    Args:
        human_seq: Human residues in the motif window (may contain gaps '-').
        animal_seq: Animal residues in the motif window (aligned, same length).
        motif_start_pos: 0-indexed start position of the motif in the full alignment.
        am_index: {position_1indexed: {alt_aa: am_score}} for the human protein.

    Returns mean score, or None if no variants could be looked up.
    """
    scores = []
    human_residue_pos = motif_start_pos  # running 0-indexed alignment position

    for i, (h_aa, a_aa) in enumerate(zip(human_seq, animal_seq)):
        if h_aa == "-":
            continue  # gap in human — skip
        human_residue_pos_1 = human_residue_pos + 1  # convert to 1-indexed
        if h_aa != a_aa and a_aa != "-":
            # Divergent position — look up AM score
            pos_scores = am_index.get(human_residue_pos_1, {})
            if a_aa in pos_scores:
                scores.append(pos_scores[a_aa])
        human_residue_pos += 1

    return round(sum(scores) / len(scores), 4) if scores else None


def annotate_motif_consequences(
    am_index: dict[str, dict[int, dict[str, float]]],
    gene_ids: Optional[list[str]] = None,
) -> int:
    """Annotate DivergentMotif rows with AlphaMissense consequence_score.

    Args:
        am_index: Full in-memory AM index (from build_am_index()).
        gene_ids: Optional filter — only process these gene IDs.

    Returns number of motifs that received a consequence_score.
    """
    scored = 0

    with get_session() as session:
        q = session.query(Gene)
        if gene_ids:
            q = q.filter(Gene.id.in_(gene_ids))
        genes = q.all()

        log.info("Scoring AlphaMissense consequences for %d genes...", len(genes))

        for gene in genes:
            accession = gene.human_protein
            if not accession:
                continue
            if "|" in accession:
                accession = accession.split("|")[-1]

            protein_am = am_index.get(accession)
            if not protein_am:
                continue

            motifs = (
                session.query(DivergentMotif)
                .join(Ortholog, DivergentMotif.ortholog_id == Ortholog.id)
                .filter(Ortholog.gene_id == gene.id)
                .all()
            )

            for motif in motifs:
                score = _consequence_for_motif(
                    motif.human_seq,
                    motif.animal_seq,
                    motif.start_pos,
                    protein_am,
                )
                if score is not None:
                    motif.consequence_score = score
                    scored += 1

        session.commit()

    log.info("AlphaMissense annotation complete: %d motifs scored.", scored)
    return scored


def run_alphamissense_pipeline(gene_ids: Optional[list[str]] = None) -> int:
    """Entry point called from orchestrator step4b.

    Downloads AlphaMissense data if needed, builds a filtered index (only the
    UniProt accessions present in our Gene table), annotates motifs.
    Gracefully skips if download fails (non-critical step).
    """
    tsv_path = download_alphamissense()
    if not tsv_path.exists():
        log.warning("AlphaMissense TSV not available — skipping consequence scoring.")
        return 0

    # Pre-compute target accessions so build_am_index skips irrelevant variants
    target_accessions = _get_target_accessions()
    am_index = build_am_index(tsv_path, target_accessions=target_accessions)
    if not am_index:
        log.warning("AlphaMissense index is empty — skipping consequence scoring.")
        return 0

    return annotate_motif_consequences(am_index, gene_ids)
