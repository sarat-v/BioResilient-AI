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
import json
import logging
import re
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import requests

from db.models import DivergentMotif, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_local_storage_root, get_tool_config

log = logging.getLogger(__name__)

AM_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
AM_CHUNK = 1024 * 1024 * 8   # 8 MB download chunks


def _am_path() -> Path:
    root = Path(get_local_storage_root()) / "alphamissense"
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


def _extract_uniprot_accession(raw: str) -> Optional[str]:
    """Extract a bare UniProt accession from various formats.

    Handles:
      - "human|P12345"          → "P12345"    (reheadered FASTA, accession preserved)
      - "human|BRCA1_HUMAN"     → None        (reheadered FASTA with mnemonic — not an accession)
      - "sp|P12345|BRCA1_HUMAN" → "P12345"    (raw UniProt FASTA header)
      - "P12345"                → "P12345"    (bare accession)
      - "P12345-2"              → "P12345"    (isoform)

    UniProt accessions match the pattern: [OPQ][0-9][A-Z0-9]{3}[0-9]
    or [A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
    A simpler heuristic: 6–10 uppercase alphanumeric chars with no underscore.
    Mnemonics look like GENE_SPECIES (contains underscore).

    Returns None if no valid accession can be extracted.
    """
    if not raw:
        return None

    # Split on pipe and examine all parts
    parts = raw.strip().split("|")

    # "sp|P12345|BRCA1_HUMAN" → check index 1 first (standard UniProt format)
    if len(parts) >= 3:
        candidate = parts[1].split("-")[0].strip()
        if re.match(r"^[A-Z][A-Z0-9]{4,9}$", candidate) and "_" not in candidate:
            return candidate.upper()

    # Try each part: take the first one that looks like an accession (no underscore, 5-10 chars)
    for part in parts:
        candidate = part.split("-")[0].strip()
        if re.match(r"^[A-Z][A-Z0-9]{4,9}$", candidate) and "_" not in candidate:
            return candidate.upper()

    return None


def _get_gene_accession_map() -> dict[str, str]:
    """Return {gene_id: uniprot_accession} for all genes in the Gene table.

    Handles both databases that store bare accessions (P12345) and those that
    store mnemonics (BRCA1_HUMAN). Mnemonics are resolved via UniProt API
    and cached on disk.
    """
    accession_map: dict[str, str] = {}
    mnemonic_genes: dict[str, str] = {}  # gene_id → mnemonic

    with get_session() as session:
        for gene in session.query(Gene).all():
            sym = (gene.gene_symbol or "").strip().upper()
            if not sym:
                continue
            # Already a bare accession?
            if re.match(r"^[A-Z][A-Z0-9]{4,9}$", sym) and "_" not in sym:
                accession_map[gene.id] = sym
            elif _is_mnemonic(sym):
                mnemonic_genes[gene.id] = sym
            else:
                acc = _extract_uniprot_accession(gene.human_protein or "")
                if acc:
                    accession_map[gene.id] = acc
                else:
                    mnemonic_genes[gene.id] = sym  # try as mnemonic anyway

    if mnemonic_genes:
        log.info("Resolving %d UniProt mnemonics to accessions...", len(mnemonic_genes))
        unique_mnemonics = list(set(mnemonic_genes.values()))
        resolved = resolve_mnemonics_to_accessions(unique_mnemonics)
        for gene_id, mnemonic in mnemonic_genes.items():
            if mnemonic in resolved:
                accession_map[gene_id] = resolved[mnemonic]

    log.info("Resolved %d / %d genes to UniProt accessions.",
             len(accession_map), len(accession_map) + len(mnemonic_genes) - sum(
                 1 for gid in mnemonic_genes if gid in accession_map))
    return accession_map


def _get_target_accessions() -> Optional[set[str]]:
    """Return the set of UniProt accessions for all genes, resolving mnemonics.

    Returns None only if the Gene table is empty or all lookups fail.
    """
    acc_map = _get_gene_accession_map()
    if not acc_map:
        log.warning("No UniProt accessions found — loading full AlphaMissense index.")
        return None
    return set(acc_map.values())


_UNSET = object()  # sentinel — distinct from None which means "load all variants"

UNIPROT_SEARCH_API = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ACCESSIONS_API = "https://rest.uniprot.org/uniprotkb/accessions"
_MNEMONIC_CACHE_FILE = "mnemonic_to_accession.json"
_MNEMONIC_BATCH = 50   # /search with OR queries; keep batches small to stay under URL limits
_MNEMONIC_TIMEOUT = 30


def _mnemonic_cache_path() -> Path:
    root = Path(get_local_storage_root()) / "alphamissense"
    root.mkdir(parents=True, exist_ok=True)
    return root / _MNEMONIC_CACHE_FILE


def resolve_mnemonics_to_accessions(mnemonics: list[str]) -> dict[str, str]:
    """Resolve UniProt mnemonics (e.g. 'BRCA1_HUMAN') to accessions ('P38398').

    Uses the UniProt /search endpoint with id: queries
    (e.g. id:BRCA1_HUMAN OR id:TP53_HUMAN) — the correct API for mnemonic lookup.
    Results are cached on disk so subsequent calls are free.

    Returns {mnemonic: accession} for every mnemonic that could be resolved.
    Unresolvable mnemonics are absent from the result dict.
    """
    cache_path = _mnemonic_cache_path()
    cache: dict[str, str] = {}
    if cache_path.exists():
        try:
            cache = json.loads(cache_path.read_text())
        except Exception:
            cache = {}

    to_fetch = [m for m in mnemonics if m not in cache]
    if not to_fetch:
        return {m: cache[m] for m in mnemonics if m in cache}

    log.info("Resolving %d UniProt mnemonics to accessions (%d cached)...",
             len(to_fetch), len(mnemonics) - len(to_fetch))

    def _fetch_batch(batch: list[str]) -> list[tuple[str, str]]:
        """Returns [(mnemonic, accession), ...] for the batch.

        Uses /search?query=id:BRCA1_HUMAN OR id:TP53_HUMAN&fields=accession,id
        which is the correct UniProt endpoint for mnemonic lookup.
        """
        query = " OR ".join(f"id:{m}" for m in batch)
        params = {
            "query": query,
            "fields": "accession,id",
            "format": "json",
            "size": len(batch),
        }
        for attempt in range(3):
            try:
                r = requests.get(UNIPROT_SEARCH_API, params=params,
                                 timeout=_MNEMONIC_TIMEOUT)
                if r.status_code == 200:
                    pairs = []
                    for entry in r.json().get("results", []):
                        acc = entry.get("primaryAccession", "").upper()
                        mnemonic = entry.get("uniProtkbId", "").upper()
                        if acc and mnemonic:
                            pairs.append((mnemonic, acc))
                    return pairs
                if r.status_code == 429:
                    time.sleep(2 * (attempt + 1))
                else:
                    log.debug("UniProt mnemonic batch HTTP %d: %s", r.status_code, r.text[:100])
                    break
            except Exception as exc:
                log.debug("UniProt mnemonic fetch error: %s", exc)
                time.sleep(attempt + 1)
        return []

    workers = int(get_tool_config().get("uniprot_batch_workers", 4))
    batches = [to_fetch[i:i + _MNEMONIC_BATCH]
               for i in range(0, len(to_fetch), _MNEMONIC_BATCH)]

    new_pairs: list[tuple[str, str]] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_fetch_batch, b): b for b in batches}
        done = 0
        for future in as_completed(futures):
            new_pairs.extend(future.result())
            done += 1
            if done % 10 == 0 or done == len(batches):
                log.info("  Mnemonic resolution: %d / %d batches done.", done, len(batches))

    for mnemonic, accession in new_pairs:
        cache[mnemonic] = accession

    try:
        cache_path.write_text(json.dumps(cache))
    except Exception as exc:
        log.debug("Mnemonic cache write failed: %s", exc)

    log.info("Resolved %d / %d mnemonics to UniProt accessions.",
             len(new_pairs), len(to_fetch))
    return {m: cache[m] for m in mnemonics if m in cache}


def _is_mnemonic(s: str) -> bool:
    """Return True if s looks like a UniProt mnemonic (GENE_SPECIES format)."""
    return bool(s and "_" in s and re.match(r"^[A-Z0-9]+_[A-Z]+$", s))


def build_am_index(
    tsv_gz_path: Path,
    target_accessions: Optional[set[str]] = _UNSET,  # type: ignore[assignment]
) -> dict[str, dict[int, dict[str, float]]]:
    """Build an in-memory index: {uniprot_id: {position: {alt_aa: score}}}.

    Args:
        tsv_gz_path: Path to AlphaMissense_hg38.tsv.gz.
        target_accessions: Set of UniProt IDs to keep. When provided, only
            variants for those proteins are loaded — reducing RAM from ~3 GB
            to typically <100 MB and cutting parse time from ~5 min to ~30 s.
            Pass None explicitly to load all ~216 M variants.
            When omitted (default), calls _get_target_accessions() which may
            return None if IDs cannot be reliably extracted, triggering the
            full load.

    Returns empty dict if file does not exist.
    """
    if not tsv_gz_path.exists():
        return {}

    if target_accessions is _UNSET:
        target_accessions = _get_target_accessions()  # may return None → load all

    if target_accessions is not None:
        n_target = len(target_accessions)
        log.info(
            "Building AlphaMissense index for %d target proteins (filtered) from %s...",
            n_target, tsv_gz_path,
        )
    else:
        log.info("Building AlphaMissense index (full, unfiltered) from %s...", tsv_gz_path)
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
    human_protein: str,
    am_index: dict[int, dict[str, float]],
) -> Optional[float]:
    """Compute mean AlphaMissense score for all divergent positions in a motif.

    Args:
        human_seq: Human residues in the motif window (may contain gaps '-').
        animal_seq: Animal residues in the motif window (aligned, same length).
        human_protein: Full human protein sequence (gap-stripped, no '*').
        am_index: {position_1indexed: {alt_aa: am_score}} for the human protein.

    Returns mean score, or None if no variants could be looked up.
    """
    # Find the true 0-indexed start position of this motif in the human protein.
    # motif.start_pos is an alignment column coordinate and must NOT be used
    # directly as a residue position — alignments may have many gaps before the
    # motif, making the column number much larger than the residue index.
    h_seq_stripped = human_seq.replace("-", "")
    if not h_seq_stripped:
        return None

    # Search near the rough alignment position first, then fall back to full search.
    protein_start = human_protein.find(h_seq_stripped)
    if protein_start < 0:
        return None  # motif not locatable in the human protein

    scores = []
    human_residue_pos = protein_start  # 0-indexed position of first non-gap residue

    for h_aa, a_aa in zip(human_seq, animal_seq):
        if h_aa == "-":
            continue  # insertion in animal — no human residue consumed
        human_residue_pos_1 = human_residue_pos + 1  # AlphaMissense is 1-indexed
        if h_aa != a_aa and a_aa != "-":
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

    Two-phase approach to prevent deadlocks with concurrently-running steps
    (e.g. step 4c ESM scoring also updates divergent_motif rows):

      Phase 1 — read-only: compute scores in Python from in-memory AM index.
                            No DB writes, no ORM autoflush, no row locks held.
      Phase 2 — single atomic write: COPY scores to a temp table, then issue
                            one bulk UPDATE.  Lock window is milliseconds, not
                            minutes, eliminating circular-wait deadlocks.

    Args:
        am_index: Full in-memory AM index (from build_am_index()).
        gene_ids: Optional filter — only process these gene IDs.

    Returns number of motifs that received a consequence_score.
    """
    import io
    import os

    import psycopg2

    # ── Phase 1: compute all scores in Python (zero DB writes) ──────────────

    scored_pairs: list[tuple[str, float]] = []   # (motif_id, score)
    gene_accession_map = _get_gene_accession_map()

    with get_session() as session:
        q = session.query(Gene)
        if gene_ids:
            q = q.filter(Gene.id.in_(gene_ids))
        genes = q.all()

        log.info("Scoring AlphaMissense consequences for %d genes...", len(genes))
        no_acc = 0
        no_data = 0

        # Pre-fetch ALL human protein sequences in one query — avoids 12,000+
        # per-gene DB round-trips.
        human_prot_q = (
            session.query(Ortholog.gene_id, Ortholog.protein_seq)
            .filter(Ortholog.species_id == "human", Ortholog.protein_seq.isnot(None))
        )
        if gene_ids:
            human_prot_q = human_prot_q.filter(Ortholog.gene_id.in_(gene_ids))
        human_protein_map: dict[str, str] = {
            row.gene_id: row.protein_seq.replace("-", "").replace("*", "")
            for row in human_prot_q.all()
            if row.protein_seq
        }
        log.info("Pre-fetched human proteins for %d genes.", len(human_protein_map))

        for gene in genes:
            accession = gene_accession_map.get(gene.id)
            if not accession:
                no_acc += 1
                continue

            protein_am = am_index.get(accession)
            if not protein_am:
                no_data += 1
                continue

            human_protein = human_protein_map.get(gene.id, "")
            if not human_protein:
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
                    human_protein,
                    protein_am,
                )
                if score is not None:
                    scored_pairs.append((str(motif.id), score))

        # Session closes here with NO commits — only reads happened.

    if no_acc:
        log.info("  %d genes had no resolvable UniProt accession.", no_acc)
    if no_data:
        log.info("  %d genes had no AlphaMissense data (accession not in TSV).", no_data)

    scored = len(scored_pairs)
    if not scored_pairs:
        log.info("AlphaMissense annotation: no motifs scored.")
        return 0

    # ── Phase 2: COPY → temp table → single bulk UPDATE ─────────────────────
    log.info("Writing %d consequence scores via bulk COPY + UPDATE...", scored)

    db_url = os.environ["DATABASE_URL"]
    conn = psycopg2.connect(db_url)
    try:
        with conn:
            with conn.cursor() as cur:
                cur.execute(
                    "CREATE TEMP TABLE _am_scores "
                    "(motif_id UUID, score DOUBLE PRECISION) ON COMMIT DROP"
                )
                buf = io.StringIO()
                for motif_id, score in scored_pairs:
                    buf.write(f"{motif_id}\t{score}\n")
                buf.seek(0)
                cur.copy_from(buf, "_am_scores", columns=("motif_id", "score"))
                cur.execute(
                    "UPDATE divergent_motif dm "
                    "SET consequence_score = s.score "
                    "FROM _am_scores s "
                    "WHERE dm.id = s.motif_id"
                )
                updated = cur.rowcount
    finally:
        conn.close()

    log.info(
        "AlphaMissense annotation complete: %d motifs scored (%d DB rows updated).",
        scored, updated,
    )
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
