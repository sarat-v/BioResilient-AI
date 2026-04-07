"""Step 8 — Functional Evidence Scoring.

Replaces GEO/DESeq2 (old step 8) and Bgee cross-species expression (old step 8b)
with three well-curated, phenotype-configurable evidence sources:

  1. Open Targets Platform  — gene-disease association score [0, 1]
                              Leverages curated genetics, transcriptomics, and
                              proteomics evidence across thousands of diseases.
  2. GTEx tissue expression — median TPM in phenotype-relevant human tissues,
                              normalised to [0, 1] via log2 scaling.
  3. DepMap gene essentiality (optional, cancer/dna_repair phenotypes only) —
                              CRISPR Chronos score converted to [0, 1] selective
                              essentiality signal.

Combined score → candidate_score.expression_score (trait_id="").
This score is used as a 4th rank-product evidence layer in Phase 1 scoring (step 9).

Performance design:
  - ENSG symbol→ID mapping: Ensembl REST bulk lookup (1000 symbols/request → ~13 calls
    for 12k genes, vs 12k individual OT search calls). Cached to disk.
  - GTEx: 20 parallel threads with short per-request delays.
  - OT association batch: 50 ENSG IDs per GraphQL request.
  - DepMap: single bulk CSV download (~150 MB), entirely local thereafter.

Per-phenotype configuration: config/functional_evidence_config.json
Reference for Open Targets: Ochoa et al. (2021) Nature Genetics.
Reference for DepMap:        Tsherniak et al. (2017) Cell.
Reference for GTEx:          GTEx Consortium (2020) Science.
Reference for Ensembl REST:  Yates et al. (2021) Nucleic Acids Research.
"""

import csv
import json
import logging
import math
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock
from typing import Optional

import requests

from db.models import CandidateScore, ExpressionResult, Gene
from db.session import get_session
from pipeline.config import get_local_storage_root

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OPENTARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
GTEX_API = "https://gtexportal.org/api/v2/expression/geneExpression"
ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
DEPMAP_URL = (
    "https://depmap.org/portal/api/download/file"
    "?file_name=CRISPRGeneDependency.csv&release=DepMap+Public+24Q2"
)

OT_BATCH_SIZE = 50        # genes per GraphQL request to Open Targets
ENSEMBL_BATCH_SIZE = 1000 # symbols per Ensembl bulk lookup request
GTEX_WORKERS = 20         # parallel threads for GTEx
REQUEST_TIMEOUT = 30      # seconds for individual HTTP requests
RATE_DELAY_OT = 0.15      # seconds between Open Targets batches
RATE_DELAY_ENSEMBL = 0.3  # seconds between Ensembl bulk requests
RATE_DELAY_GTEX = 0.1     # delay per GTEx worker (gentle per-thread rate limit)

_REPO_ROOT = Path(__file__).resolve().parents[2]
_CONFIG_PATH = _REPO_ROOT / "config" / "functional_evidence_config.json"


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def _load_phenotype_config(phenotype: str) -> dict:
    """Load per-phenotype functional evidence config, falling back to cancer_resistance."""
    try:
        with open(_CONFIG_PATH) as f:
            all_cfg = json.load(f)
    except FileNotFoundError:
        log.warning("Functional evidence config not found at %s — using defaults.", _CONFIG_PATH)
        return {}

    cfg = all_cfg.get(phenotype)
    if cfg is None:
        log.warning(
            "No config for phenotype %r — falling back to cancer_resistance.", phenotype
        )
        cfg = all_cfg.get("cancer_resistance", {})
    return cfg


# ---------------------------------------------------------------------------
# Gene symbol normalisation
# ---------------------------------------------------------------------------

def _hgnc_symbol(gene: Gene) -> Optional[str]:
    """Extract an approximate HGNC gene symbol from a Gene DB record.

    The DB stores UniProt mnemonic entry names like ``TP53_HUMAN``.
    Stripping the ``_HUMAN`` suffix recovers the HGNC symbol for the vast
    majority of well-characterised human genes (TP53, BRCA1, ATM, etc.).
    """
    sym = gene.gene_symbol or ""
    if sym.endswith("_HUMAN"):
        return sym[: -len("_HUMAN")]
    return sym if sym else None


# ---------------------------------------------------------------------------
# ENSG ID mapping via Ensembl bulk REST API
# ---------------------------------------------------------------------------

def _ensembl_bulk_lookup(symbols: list[str]) -> dict[str, str]:
    """Resolve HGNC symbols to Ensembl IDs via the Ensembl REST bulk endpoint.

    POST /lookup/symbol/homo_sapiens accepts up to 1000 symbols at once.
    Returns {symbol: ensembl_id}. Missing/invalid symbols are absent.
    """
    if not symbols:
        return {}

    result: dict[str, str] = {}
    try:
        r = requests.post(
            ENSEMBL_LOOKUP_URL,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            json={"symbols": symbols},
            timeout=60,
        )
        if r.status_code != 200:
            log.warning("Ensembl bulk lookup HTTP %d for %d symbols", r.status_code, len(symbols))
            return {}
        data = r.json()
        for sym, info in data.items():
            if isinstance(info, dict):
                ensg = info.get("id", "")
                if ensg.startswith("ENSG"):
                    result[sym] = ensg
    except Exception as exc:
        log.warning("Ensembl bulk lookup failed: %s", exc)
    return result


def _build_ensg_map(genes: list[Gene]) -> dict[str, str]:
    """Build {hgnc_symbol: ensembl_id} mapping using the Ensembl bulk REST API.

    Uses a persistent disk cache to avoid repeated lookups across pipeline runs.
    On first call for 12k genes: ~13 HTTP requests (1000 symbols each) → ~5 seconds.
    On subsequent calls: instant (cache hit).
    """
    cache_path = Path(get_local_storage_root()) / "ot_ensg_cache.json"
    cache: dict[str, str] = {}
    if cache_path.exists():
        try:
            with open(cache_path) as f:
                cache = json.load(f)
        except Exception:
            cache = {}

    # Collect symbols not yet in cache
    symbols_needed = []
    for gene in genes:
        sym = _hgnc_symbol(gene)
        if sym and sym not in cache:
            symbols_needed.append(sym)

    if symbols_needed:
        log.info(
            "ENSG lookup: %d new symbols via Ensembl bulk API (~%d requests)...",
            len(symbols_needed),
            math.ceil(len(symbols_needed) / ENSEMBL_BATCH_SIZE),
        )
        n_resolved = 0
        for i in range(0, len(symbols_needed), ENSEMBL_BATCH_SIZE):
            batch = symbols_needed[i: i + ENSEMBL_BATCH_SIZE]
            found = _ensembl_bulk_lookup(batch)
            n_resolved += len(found)
            # Cache every symbol (resolved or not) to avoid re-querying unknowns
            for sym in batch:
                cache[sym] = found.get(sym, "")
            time.sleep(RATE_DELAY_ENSEMBL)
            batch_num = i // ENSEMBL_BATCH_SIZE + 1
            total_batches = math.ceil(len(symbols_needed) / ENSEMBL_BATCH_SIZE)
            log.info(
                "  ENSG: batch %d / %d — %d resolved so far",
                batch_num, total_batches, n_resolved,
            )

        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(cache, f)
        log.info("ENSG map: %d / %d symbols resolved.", n_resolved, len(symbols_needed))

    return cache


# ---------------------------------------------------------------------------
# Open Targets
# ---------------------------------------------------------------------------

_OT_QUERY = """
query FunctionalEvidence($ensemblIds: [String!]!) {
  targets(ensemblIds: $ensemblIds) {
    id
    approvedSymbol
    associatedDiseases(size: 200) {
      rows {
        score
        disease { id name }
      }
    }
  }
}
"""


def _fetch_ot_batch(ensg_ids: list[str], disease_set: set[str]) -> dict[str, float]:
    """Fetch Open Targets association scores for a batch of ENSEMBL IDs.

    Returns {ensembl_id: max_score_across_requested_diseases}.
    """
    valid = [e for e in ensg_ids if e and e.startswith("ENSG")]
    if not valid:
        return {}

    try:
        r = requests.post(
            OPENTARGETS_GRAPHQL,
            json={"query": _OT_QUERY, "variables": {"ensemblIds": valid}},
            timeout=REQUEST_TIMEOUT,
        )
        r.raise_for_status()
        data = r.json()
    except Exception as exc:
        log.warning("Open Targets batch failed: %s", exc)
        return {}

    result: dict[str, float] = {}
    for target in (data.get("data") or {}).get("targets") or []:
        eid = target.get("id", "")
        scores = [
            float(row["score"])
            for row in (target.get("associatedDiseases") or {}).get("rows") or []
            if row.get("disease", {}).get("id") in disease_set and row.get("score") is not None
        ]
        if scores:
            result[eid] = max(scores)
    return result


def fetch_open_targets_scores(
    genes: list[Gene],
    disease_ids: list[str],
) -> dict[str, float]:
    """Fetch Open Targets gene-disease association scores for all genes.

    Returns {gene_id (DB primary key): score [0, 1]}.
    """
    if not disease_ids:
        return {}

    disease_set = set(disease_ids)
    log.info(
        "Open Targets: querying %d genes against %d disease IDs...",
        len(genes), len(disease_ids),
    )

    # Build symbol → ENSG map via Ensembl bulk API (cached)
    ensg_map = _build_ensg_map(genes)

    # Map gene DB id → ENSG
    gid_to_ensg: dict[str, str] = {}
    for gene in genes:
        sym = _hgnc_symbol(gene)
        if sym:
            ensg = ensg_map.get(sym, "")
            if ensg:
                gid_to_ensg[gene.id] = ensg

    if not gid_to_ensg:
        log.info("No ENSG IDs resolved — skipping Open Targets scoring.")
        return {}

    log.info("  %d / %d genes have resolved ENSG IDs.", len(gid_to_ensg), len(genes))

    ensg_to_gid: dict[str, str] = {v: k for k, v in gid_to_ensg.items()}
    ensg_list = list(ensg_to_gid.keys())
    n_batches = math.ceil(len(ensg_list) / OT_BATCH_SIZE)

    all_scores: dict[str, float] = {}
    for i in range(0, len(ensg_list), OT_BATCH_SIZE):
        batch = ensg_list[i: i + OT_BATCH_SIZE]
        scores = _fetch_ot_batch(batch, disease_set)
        for ensg, score in scores.items():
            gid = ensg_to_gid.get(ensg)
            if gid:
                all_scores[gid] = score
        time.sleep(RATE_DELAY_OT)
        batch_num = i // OT_BATCH_SIZE + 1
        if batch_num % 10 == 0 or batch_num == n_batches:
            log.info("  OT: %d / %d batches complete", batch_num, n_batches)

    log.info("Open Targets: %d / %d genes scored.", len(all_scores), len(genes))
    return all_scores


# ---------------------------------------------------------------------------
# GTEx  (parallelised with ThreadPoolExecutor)
# ---------------------------------------------------------------------------

def _fetch_gtex_gene(gene_symbol: str, tissue_ids: list[str]) -> Optional[dict[str, float]]:
    """Fetch GTEx v10 median TPM for a gene across the requested tissues.

    Issues one bulk request (all tissues, filter client-side).
    Returns {tissue_id: median_tpm} or None if no data.
    """
    result: dict[str, float] = {}
    tissue_set = set(tissue_ids)

    try:
        r = requests.get(
            GTEX_API,
            params={"geneSymbol": gene_symbol, "datasetId": "gtex_v10"},
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code == 200:
            for row in r.json().get("data", []):
                tissue = row.get("tissueSiteDetailId") or row.get("tissueSiteDetail", "")
                if tissue in tissue_set:
                    tpm = row.get("median")
                    if tpm is not None:
                        try:
                            result[tissue] = float(tpm)
                        except (TypeError, ValueError):
                            pass
    except Exception as exc:
        log.debug("GTEx fetch failed for %s: %s", gene_symbol, exc)

    return result if result else None


def _gtex_normalize(tpm_values: list[float]) -> float:
    """Normalise a list of tissue TPM values to [0, 1].

    Uses log2(1 + mean_tpm) / 10, capped at 1.0.
    Interpretation: mean TPM ≥ 1023 → score 1.0; mean TPM 0 → score 0.0.
    """
    if not tpm_values:
        return 0.0
    mean_tpm = sum(tpm_values) / len(tpm_values)
    return round(min(math.log2(1.0 + mean_tpm) / 10.0, 1.0), 4)


def fetch_gtex_scores(
    genes: list[Gene],
    tissue_ids: list[str],
) -> dict[str, float]:
    """Fetch GTEx expression scores for all genes in phenotype-relevant tissues.

    Uses a thread pool (GTEX_WORKERS=20) so 12k genes complete in ~2 minutes
    instead of ~26 minutes sequentially.

    Returns {gene_id (DB primary key): score [0, 1]}.
    """
    if not tissue_ids:
        return {}

    # Filter to genes with valid symbols
    gene_pairs = [(g, _hgnc_symbol(g)) for g in genes]
    gene_pairs = [(g, sym) for g, sym in gene_pairs if sym]
    log.info("GTEx: querying %d genes × %d tissues with %d workers...",
             len(gene_pairs), len(tissue_ids), GTEX_WORKERS)

    scores: dict[str, float] = {}
    lock = Lock()
    done_count = [0]

    def _worker(g: Gene, sym: str) -> None:
        tpm_map = _fetch_gtex_gene(sym, tissue_ids)
        time.sleep(RATE_DELAY_GTEX)
        if tpm_map:
            score = _gtex_normalize(list(tpm_map.values()))
            if score > 0:
                with lock:
                    scores[g.id] = score
        with lock:
            done_count[0] += 1
            if done_count[0] % 500 == 0:
                log.info("  GTEx: %d / %d genes done...", done_count[0], len(gene_pairs))

    with ThreadPoolExecutor(max_workers=GTEX_WORKERS) as pool:
        futures = [pool.submit(_worker, g, sym) for g, sym in gene_pairs]
        for f in as_completed(futures):
            exc = f.exception()
            if exc:
                log.debug("GTEx worker error: %s", exc)

    log.info("GTEx: %d / %d genes scored.", len(scores), len(gene_pairs))
    return scores


# ---------------------------------------------------------------------------
# DepMap
# ---------------------------------------------------------------------------

def _depmap_cache_path() -> Path:
    root = Path(get_local_storage_root()) / "depmap"
    root.mkdir(parents=True, exist_ok=True)
    return root / "CRISPRGeneDependency.csv"


def _load_depmap_index() -> dict[str, float]:
    """Return {gene_symbol: mean_chronos_score} from the DepMap public CSV.

    Downloads the CSV (~150 MB) on first call and caches it locally.
    Chronos score interpretation: ~-1 = essential, ~0 = non-essential.
    """
    path = _depmap_cache_path()
    if not path.exists():
        log.info("Downloading DepMap CRISPR scores (~150 MB)...")
        try:
            with requests.get(DEPMAP_URL, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        f.write(chunk)
            log.info("DepMap download complete.")
        except Exception as exc:
            log.warning("DepMap download failed: %s — skipping essentiality scoring.", exc)
            if path.exists():
                path.unlink()
            return {}

    log.info("Building DepMap index from %s...", path)
    gene_scores: dict[str, list[float]] = {}

    try:
        with open(path, newline="", errors="replace") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if not header:
                return {}
            symbols = [h.split(" (")[0].strip() for h in header[1:]]
            for sym in symbols:
                gene_scores[sym] = []
            for row in reader:
                for i, val in enumerate(row[1:]):
                    if i >= len(symbols):
                        break
                    try:
                        gene_scores[symbols[i]].append(float(val))
                    except (ValueError, TypeError):
                        pass
    except Exception as exc:
        log.warning("DepMap parse failed: %s", exc)
        return {}

    result = {
        sym: round(sum(sc) / len(sc), 4)
        for sym, sc in gene_scores.items() if sc
    }
    log.info("DepMap index built: %d genes.", len(result))
    return result


def _depmap_to_score(chronos: float) -> float:
    """Convert mean Chronos score to a [0, 1] functional relevance score.

    Genes with strong selective essentiality (chronos << -0.5) score near 1.
    Non-essential genes (chronos ≈ 0) score near 0.
    Uses a sigmoid centred at -0.5: score = σ(-10 × (chronos + 0.5)).
    """
    return round(1.0 / (1.0 + math.exp(10.0 * (chronos + 0.5))), 4)


def fetch_depmap_scores(genes: list[Gene]) -> dict[str, float]:
    """Fetch DepMap essentiality-based scores for all genes.

    Returns {gene_id (DB primary key): score [0, 1]}.
    """
    log.info("DepMap: loading CRISPR essentiality index...")
    index = _load_depmap_index()
    if not index:
        return {}

    scores: dict[str, float] = {}
    for gene in genes:
        sym = _hgnc_symbol(gene)
        if sym and sym in index:
            score = _depmap_to_score(index[sym])
            if score > 0:
                scores[gene.id] = score

    log.info("DepMap: %d / %d genes scored.", len(scores), len(genes))
    return scores


# ---------------------------------------------------------------------------
# Score combination
# ---------------------------------------------------------------------------

def _combine_scores(
    gene_ids: list[str],
    ot_scores: dict[str, float],
    gtex_scores: dict[str, float],
    depmap_scores: dict[str, float],
) -> dict[str, float]:
    """Weighted mean of available sub-scores per gene.

    Default weights: Open Targets=0.5, GTEx=0.3, DepMap=0.2.
    When a source is absent for a gene, its weight is redistributed
    proportionally across the remaining sources.

    Returns {gene_id: combined_score [0, 1]}.
    """
    source_weights = [
        ("ot",     ot_scores,     0.5),
        ("gtex",   gtex_scores,   0.3),
        ("depmap", depmap_scores, 0.2),
    ]

    combined: dict[str, float] = {}
    for gid in gene_ids:
        parts = [(w, s[gid]) for _, s, w in source_weights if gid in s]
        if not parts:
            combined[gid] = 0.0
            continue
        total_w = sum(w for w, _ in parts)
        combined[gid] = round(sum(w * v for w, v in parts) / total_w, 4)

    return combined


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def _persist_results(
    gene_ids: list[str],
    combined_scores: dict[str, float],
    ot_scores: dict[str, float],
    gtex_scores: dict[str, float],
    depmap_scores: dict[str, float],
) -> int:
    """Write evidence sub-scores to expression_result and update candidate_score.

    Writes to trait_id="" (the default used by step 9 scoring) so that the
    functional evidence score feeds directly into the rank-product in step 9.
    Phenotype info is preserved in expression_result.comparison and .geo_accession.

    Returns number of genes with a non-zero combined score.
    """
    from sqlalchemy import delete

    rows: list[ExpressionResult] = []
    for gid in gene_ids:
        if gid in ot_scores:
            rows.append(ExpressionResult(
                gene_id=gid,
                geo_accession="OT:disease_association",
                comparison="OpenTargets gene-disease association score",
                log2fc=round(ot_scores[gid], 4),
                padj=None,
                n_samples=None,
            ))
        if gid in gtex_scores:
            rows.append(ExpressionResult(
                gene_id=gid,
                geo_accession="GTEX:tissue_expression",
                comparison="GTEx median TPM in phenotype-relevant tissues (normalised)",
                log2fc=round(gtex_scores[gid], 4),
                padj=None,
                n_samples=None,
            ))
        if gid in depmap_scores:
            rows.append(ExpressionResult(
                gene_id=gid,
                geo_accession="DEPMAP:essentiality",
                comparison="DepMap CRISPR Chronos selective essentiality",
                log2fc=round(depmap_scores[gid], 4),
                padj=None,
                n_samples=None,
            ))

    nonzero = sum(1 for v in combined_scores.values() if v > 0)

    with get_session() as session:
        # Clear previous expression_result rows
        session.execute(delete(ExpressionResult))
        if rows:
            session.add_all(rows)

        # Write expression_score to CandidateScore rows with trait_id=""
        # (matches what step 9 run_scoring uses as default trait_id)
        cs_map = {
            r.gene_id: r
            for r in session.query(CandidateScore).filter_by(trait_id="").all()
        }

        updated = 0
        for gid, score in combined_scores.items():
            if score <= 0:
                continue
            cs = cs_map.get(gid)
            if cs is None:
                cs = CandidateScore(gene_id=gid, trait_id="")
                session.add(cs)
            cs.expression_score = score
            updated += 1

        session.commit()

    log.info(
        "Persisted %d expression_result rows; updated expression_score for %d genes.",
        len(rows), updated,
    )
    return nonzero


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_functional_evidence(phenotype: str = "cancer_resistance") -> dict:
    """Run the full functional evidence pipeline for all genes in the DB.

    Steps:
      1. Load per-phenotype config (OT disease IDs, GTEx tissues, DepMap flag).
      2. Fetch Open Targets association scores (Ensembl bulk ENSG lookup + batched GraphQL).
      3. Fetch GTEx tissue expression scores (parallelised with 20 threads).
      4. Optionally fetch DepMap essentiality scores.
      5. Combine into a single expression_score per gene.
      6. Persist to expression_result table and candidate_score.expression_score (trait_id="").

    Returns a summary dict with counts for monitoring/reporting.
    """
    cfg = _load_phenotype_config(phenotype)
    if not cfg:
        log.warning("Empty config for phenotype %r — nothing to do.", phenotype)
        return {"phenotype": phenotype, "n_genes": 0, "n_ot": 0, "n_gtex": 0, "n_depmap": 0}

    with get_session() as session:
        genes = session.query(Gene).all()

    if not genes:
        log.warning("No genes in DB — skipping functional evidence step.")
        return {"phenotype": phenotype, "n_genes": 0, "n_ot": 0, "n_gtex": 0, "n_depmap": 0}

    gene_ids = [g.id for g in genes]
    log.info(
        "Functional evidence (step 8): %d genes, phenotype=%r", len(genes), phenotype
    )

    # Open Targets -------------------------------------------------------
    ot_cfg = cfg.get("opentargets", {})
    ot_scores: dict[str, float] = {}
    if ot_cfg.get("enabled", True):
        disease_ids = ot_cfg.get("disease_ids", [])
        if disease_ids:
            ot_scores = fetch_open_targets_scores(genes, disease_ids)
        else:
            log.info("Open Targets: no disease_ids configured — skipping.")
    else:
        log.info("Open Targets: disabled for phenotype %r.", phenotype)

    # GTEx ---------------------------------------------------------------
    gtex_cfg = cfg.get("gtex", {})
    gtex_scores: dict[str, float] = {}
    if gtex_cfg.get("enabled", True):
        tissues = gtex_cfg.get("tissues", [])
        if tissues:
            gtex_scores = fetch_gtex_scores(genes, tissues)
        else:
            log.info("GTEx: no tissues configured — skipping.")
    else:
        log.info("GTEx: disabled for phenotype %r.", phenotype)

    # DepMap -------------------------------------------------------------
    depmap_cfg = cfg.get("depmap", {})
    depmap_scores: dict[str, float] = {}
    if depmap_cfg.get("enabled", False):
        depmap_scores = fetch_depmap_scores(genes)
    else:
        log.info("DepMap: disabled for phenotype %r.", phenotype)

    # Combine and persist ------------------------------------------------
    combined = _combine_scores(gene_ids, ot_scores, gtex_scores, depmap_scores)
    nonzero = _persist_results(gene_ids, combined, ot_scores, gtex_scores, depmap_scores)

    summary = {
        "phenotype": phenotype,
        "n_genes": len(genes),
        "n_ot": len(ot_scores),
        "n_gtex": len(gtex_scores),
        "n_depmap": len(depmap_scores),
        "n_scored": nonzero,
    }
    log.info("Functional evidence complete: %s", summary)
    return summary
