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

Scope: Only human orthologs (gene_symbol ending in _HUMAN) are queried against
OT, GTEx, and DepMap. These are human-centric databases and only human gene symbols
return valid results. Functional evidence scores are propagated for each orthogroup
via the human representative gene's DB entry.

Performance design:
  - ENSG symbol→ID mapping: Ensembl REST bulk lookup (1000 symbols/request → ~6 calls
    for ~5800 human genes). Cached to disk.
  - GTEx gencodeId lookup: GTEx reference/gene endpoint (per-gene, 20 parallel workers,
    cached to disk). Required because GTEx medianGeneExpression needs versioned gencodeId
    (e.g. ENSG00000141510.18 for Gencode v39, not .20 from current Ensembl).
  - GTEx expression: medianGeneExpression endpoint with parallel workers.
  - OT association batch: 50 ENSG IDs per GraphQL request.
  - DepMap: single bulk CSV download (~400 MB, cached locally).

References:
  Open Targets: Ochoa et al. (2021) Nature Genetics.
  DepMap:       Tsherniak et al. (2017) Cell. Data: 24Q4 Public.
  GTEx:         GTEx Consortium (2020) Science. API v2, dataset gtex_v10 (Gencode v39).
  Ensembl REST: Yates et al. (2021) Nucleic Acids Research.
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
GTEX_EXPRESSION_API = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
GTEX_GENE_REF_API = "https://gtexportal.org/api/v2/reference/gene"
ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
# DepMap 24Q4 Public — CRISPRGeneEffect.csv (Chronos scores, ~408 MB) via Figshare.
# NOTE: CRISPRGeneEffect.csv contains raw Chronos scores (negative = essential, ~0 = non-essential).
# Do NOT use CRISPRGeneDependency.csv — those are probability scores (0-1) incompatible
# with the Chronos-calibrated sigmoid in _depmap_to_score().
DEPMAP_URL = "https://ndownloader.figshare.com/files/51064667"
DEPMAP_HEADERS = {"User-Agent": "BioResilient research pipeline (Mozilla/5.0 compatible)"}

OT_BATCH_SIZE = 50        # genes per GraphQL request to Open Targets
ENSEMBL_BATCH_SIZE = 1000 # symbols per Ensembl bulk lookup request
GTEX_WORKERS = 20         # parallel threads for GTEx gencodeId lookup + expression
REQUEST_TIMEOUT = 30      # seconds for individual HTTP requests
RATE_DELAY_OT = 0.15      # seconds between Open Targets batches
RATE_DELAY_ENSEMBL = 0.3  # seconds between Ensembl bulk requests
RATE_DELAY_GTEX = 0.05    # per-worker GTEx delay (gentle rate limit)

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
# Gene symbol helpers
# ---------------------------------------------------------------------------

def _hgnc_symbol(gene: Gene) -> Optional[str]:
    """Extract the HGNC gene symbol from a Gene record.

    Only returns a symbol for human genes (UniProt mnemonic ending in _HUMAN,
    e.g. TP53_HUMAN → TP53). Returns None for all other species since OT,
    GTEx, and DepMap are human-centric databases.
    """
    sym = gene.gene_symbol or ""
    if sym.endswith("_HUMAN"):
        return sym[: -len("_HUMAN")]
    return None


def _human_genes(genes: list[Gene]) -> list[Gene]:
    """Filter to human genes only (those with HGNC-resolvable symbols)."""
    return [g for g in genes if _hgnc_symbol(g) is not None]


# ---------------------------------------------------------------------------
# ENSG ID mapping via Ensembl bulk REST API
# ---------------------------------------------------------------------------

def _ensembl_bulk_lookup(symbols: list[str]) -> dict[str, str]:
    """Resolve HGNC symbols to Ensembl stable IDs via the Ensembl REST bulk endpoint.

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
    """Build {hgnc_symbol: ensembl_stable_id} for human genes.

    Uses Ensembl REST bulk API (1000 symbols/request). Cached to disk.
    Returns bare ENSG IDs without version suffix (e.g. ENSG00000141510).
    """
    cache_path = Path(get_local_storage_root()) / "ot_ensg_cache.json"
    cache: dict[str, str] = {}
    if cache_path.exists():
        try:
            with open(cache_path) as f:
                cache = json.load(f)
        except Exception:
            cache = {}

    symbols_needed = [
        sym for g in genes
        if (sym := _hgnc_symbol(g)) and sym not in cache
    ]

    if symbols_needed:
        log.info(
            "ENSG lookup: %d symbols via Ensembl bulk API (~%d requests)...",
            len(symbols_needed),
            math.ceil(len(symbols_needed) / ENSEMBL_BATCH_SIZE),
        )
        n_resolved = 0
        for i in range(0, len(symbols_needed), ENSEMBL_BATCH_SIZE):
            batch = symbols_needed[i: i + ENSEMBL_BATCH_SIZE]
            found = _ensembl_bulk_lookup(batch)
            n_resolved += len(found)
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

# OT GraphQL query — uses page: {index: 0, size: N} not the deprecated size arg
_OT_QUERY = """
query FunctionalEvidence($ensemblIds: [String!]!) {
  targets(ensemblIds: $ensemblIds) {
    id
    approvedSymbol
    associatedDiseases(page: {index: 0, size: 200}) {
      rows {
        score
        disease { id name }
      }
    }
  }
}
"""


def _fetch_ot_batch(ensg_ids: list[str], disease_set: set[str]) -> dict[str, float]:
    """Fetch Open Targets association scores for a batch of ENSG IDs.

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
    """Fetch Open Targets gene-disease association scores for human genes.

    Returns {gene_id (DB primary key): score [0, 1]}.
    """
    if not disease_ids:
        return {}

    human = _human_genes(genes)
    if not human:
        return {}

    disease_set = set(disease_ids)
    log.info(
        "Open Targets: querying %d human genes against %d disease IDs...",
        len(human), len(disease_ids),
    )

    ensg_map = _build_ensg_map(human)

    gid_to_ensg: dict[str, str] = {}
    for gene in human:
        sym = _hgnc_symbol(gene)
        if sym:
            ensg = ensg_map.get(sym, "")
            if ensg:
                gid_to_ensg[gene.id] = ensg

    if not gid_to_ensg:
        log.info("No ENSG IDs resolved — skipping Open Targets scoring.")
        return {}

    log.info("  %d / %d human genes have ENSG IDs.", len(gid_to_ensg), len(human))

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

    log.info("Open Targets: %d / %d human genes scored.", len(all_scores), len(human))
    return all_scores


# ---------------------------------------------------------------------------
# GTEx  (versioned gencodeId + parallel workers)
# ---------------------------------------------------------------------------

def _build_gtex_gencode_map(symbols: list[str]) -> dict[str, str]:
    """Build {gene_symbol: versioned_gencodeId} for GTEx v10 (Gencode v39).

    GTEx medianGeneExpression requires versioned gencodeId (e.g. ENSG00000141510.18),
    not the bare ENSG ID. Uses GTEx reference/gene endpoint with gencodeVersion=v39.
    Parallelised with 20 threads and cached to disk.
    """
    cache_path = Path(get_local_storage_root()) / "gtex_gencode_cache.json"
    cache: dict[str, str] = {}
    if cache_path.exists():
        try:
            with open(cache_path) as f:
                cache = json.load(f)
        except Exception:
            cache = {}

    needed = [s for s in symbols if s not in cache]
    if not needed:
        return cache

    log.info(
        "GTEx gencodeId lookup: %d symbols with %d workers...", len(needed), GTEX_WORKERS
    )

    lock = Lock()
    resolved = [0]

    def _lookup(sym: str) -> None:
        try:
            r = requests.get(
                GTEX_GENE_REF_API,
                params={
                    "geneId": sym,
                    "gencodeVersion": "v39",
                    "genomeBuild": "GRCh38/hg38",
                },
                timeout=REQUEST_TIMEOUT,
            )
            if r.status_code == 200:
                rows = r.json().get("data", [])
                if rows:
                    gencode_id = rows[0].get("gencodeId", "")
                    with lock:
                        cache[sym] = gencode_id
                        resolved[0] += 1
                    return
        except Exception as exc:
            log.debug("GTEx ref lookup failed for %s: %s", sym, exc)
        with lock:
            cache[sym] = ""  # Mark as attempted even on failure
        time.sleep(RATE_DELAY_GTEX)

    with ThreadPoolExecutor(max_workers=GTEX_WORKERS) as pool:
        futures = [pool.submit(_lookup, s) for s in needed]
        for f in as_completed(futures):
            pass  # Errors handled in _lookup

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(cache, f)
    log.info("GTEx gencodeId map: %d / %d symbols resolved.", resolved[0], len(needed))
    return cache


def _fetch_gtex_expression(gencode_id: str, tissue_set: set[str]) -> Optional[dict[str, float]]:
    """Fetch GTEx v10 median TPM for a gene across the requested tissues.

    Requires versioned gencodeId (e.g. ENSG00000141510.18).
    Returns {tissue_id: median_tpm} or None if no data.
    """
    if not gencode_id:
        return None

    result: dict[str, float] = {}
    try:
        r = requests.get(
            GTEX_EXPRESSION_API,
            params={"gencodeId": gencode_id, "datasetId": "gtex_v10"},
            timeout=REQUEST_TIMEOUT,
        )
        if r.status_code == 200:
            for row in r.json().get("data", []):
                tissue = row.get("tissueSiteDetailId", "")
                if tissue in tissue_set:
                    tpm = row.get("median")
                    if tpm is not None:
                        try:
                            result[tissue] = float(tpm)
                        except (TypeError, ValueError):
                            pass
    except Exception as exc:
        log.debug("GTEx expression fetch failed for %s: %s", gencode_id, exc)

    return result if result else None


def _gtex_normalize(tpm_values: list[float]) -> float:
    """Normalise tissue TPM values to [0, 1] using log2(1 + mean_tpm) / 10, capped at 1."""
    if not tpm_values:
        return 0.0
    mean_tpm = sum(tpm_values) / len(tpm_values)
    return round(min(math.log2(1.0 + mean_tpm) / 10.0, 1.0), 4)


def fetch_gtex_scores(
    genes: list[Gene],
    tissue_ids: list[str],
) -> dict[str, float]:
    """Fetch GTEx expression scores for human genes in phenotype-relevant tissues.

    Process:
      1. Build symbol → versioned gencodeId map (parallel, cached).
      2. Fetch medianGeneExpression for each gene (parallel 20 workers).

    Returns {gene_id (DB primary key): score [0, 1]}.
    """
    if not tissue_ids:
        return {}

    human = _human_genes(genes)
    if not human:
        return {}

    tissue_set = set(tissue_ids)
    symbols = [sym for g in human if (sym := _hgnc_symbol(g))]

    log.info(
        "GTEx: building gencodeId map for %d human gene symbols...", len(symbols)
    )
    gencode_map = _build_gtex_gencode_map(symbols)

    # Build (gene, gencode_id) pairs for expression lookup
    gene_gencode_pairs = []
    for gene in human:
        sym = _hgnc_symbol(gene)
        if sym:
            gid = gencode_map.get(sym, "")
            if gid:
                gene_gencode_pairs.append((gene, gid))

    log.info(
        "GTEx: fetching expression for %d / %d human genes with %d workers...",
        len(gene_gencode_pairs), len(human), GTEX_WORKERS,
    )

    scores: dict[str, float] = {}
    lock = Lock()
    done_count = [0]

    def _worker(g: Gene, gencode_id: str) -> None:
        tpm_map = _fetch_gtex_expression(gencode_id, tissue_set)
        time.sleep(RATE_DELAY_GTEX)
        if tpm_map:
            score = _gtex_normalize(list(tpm_map.values()))
            if score > 0:
                with lock:
                    scores[g.id] = score
        with lock:
            done_count[0] += 1
            n = done_count[0]
            total = len(gene_gencode_pairs)
            if n % 500 == 0 or n == total:
                log.info("  GTEx expression: %d / %d genes done...", n, total)

    with ThreadPoolExecutor(max_workers=GTEX_WORKERS) as pool:
        futures = [pool.submit(_worker, g, gid) for g, gid in gene_gencode_pairs]
        for f in as_completed(futures):
            exc = f.exception()
            if exc:
                log.debug("GTEx worker error: %s", exc)

    log.info("GTEx: %d / %d human genes scored.", len(scores), len(human))
    return scores


# ---------------------------------------------------------------------------
# DepMap
# ---------------------------------------------------------------------------

def _depmap_cache_path() -> Path:
    root = Path(get_local_storage_root()) / "depmap"
    root.mkdir(parents=True, exist_ok=True)
    return root / "CRISPRGeneEffect_24Q4.csv"


def _load_depmap_index() -> dict[str, float]:
    """Return {gene_symbol: mean_chronos_score} from DepMap 24Q4 CRISPRGeneEffect.

    CRISPRGeneEffect.csv stores raw Chronos scores per gene per cell line.
    Chronos score: ~-1 = strongly essential, ~0 = non-essential, rarely positive.
    Mean across all cell lines gives the typical essentiality across cancer lines.
    Downloads (~408 MB) on first call and caches locally.
    """
    path = _depmap_cache_path()
    if not path.exists():
        log.info("Downloading DepMap 24Q4 CRISPRGeneEffect (Chronos scores, ~408 MB)...")
        try:
            with requests.get(
                DEPMAP_URL, headers=DEPMAP_HEADERS, stream=True, timeout=300
            ) as r:
                r.raise_for_status()
                with open(path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        f.write(chunk)
            log.info("DepMap download complete: %s", path)
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
            # Column format: "SYMBOL (ENTREZ_ID)" → extract symbol
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
    """Convert mean Chronos score to [0, 1] selective essentiality.

    Sigmoid centred at -0.5: σ(-10 × (chronos + 0.5)).
    Chronos < -0.5 (strongly essential) → score > 0.5.
    """
    return round(1.0 / (1.0 + math.exp(10.0 * (chronos + 0.5))), 4)


def fetch_depmap_scores(genes: list[Gene]) -> dict[str, float]:
    """Fetch DepMap CRISPR essentiality scores for human genes.

    Returns {gene_id (DB primary key): score [0, 1]}.
    """
    human = _human_genes(genes)
    if not human:
        return {}

    log.info("DepMap: loading CRISPR essentiality index for %d human genes...", len(human))
    index = _load_depmap_index()
    if not index:
        return {}

    scores: dict[str, float] = {}
    for gene in human:
        sym = _hgnc_symbol(gene)
        if sym and sym in index:
            score = _depmap_to_score(index[sym])
            if score > 0:
                scores[gene.id] = score

    log.info("DepMap: %d / %d human genes scored.", len(scores), len(human))
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

    Writes CandidateScore with trait_id="" (the default used by step 9 run_scoring)
    so the functional evidence score feeds directly into the rank-product.

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
                comparison="GTEx v10 median TPM in phenotype-relevant tissues (normalised)",
                log2fc=round(gtex_scores[gid], 4),
                padj=None,
                n_samples=None,
            ))
        if gid in depmap_scores:
            rows.append(ExpressionResult(
                gene_id=gid,
                geo_accession="DEPMAP:essentiality",
                comparison="DepMap 24Q4 CRISPRGeneEffect Chronos selective essentiality",
                log2fc=round(depmap_scores[gid], 4),
                padj=None,
                n_samples=None,
            ))

    nonzero = sum(1 for v in combined_scores.values() if v > 0)

    with get_session() as session:
        session.execute(delete(ExpressionResult))
        if rows:
            session.add_all(rows)

        # Write expression_score to CandidateScore rows (trait_id="" = step 9 default)
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
      2. Filter to human genes (_HUMAN suffix in UniProt mnemonic).
      3. Fetch Open Targets association scores (Ensembl bulk ENSG + batched GraphQL).
      4. Fetch GTEx expression scores (gencodeId lookup + parallelised expression fetch).
      5. Optionally fetch DepMap essentiality scores.
      6. Combine into a single expression_score per gene.
      7. Persist to expression_result and candidate_score.expression_score (trait_id="").

    Returns summary dict for monitoring/reporting.
    """
    cfg = _load_phenotype_config(phenotype)
    if not cfg:
        log.warning("Empty config for phenotype %r — nothing to do.", phenotype)
        return {"phenotype": phenotype, "n_genes": 0, "n_human": 0,
                "n_ot": 0, "n_gtex": 0, "n_depmap": 0, "n_scored": 0}

    with get_session() as session:
        genes = session.query(Gene).all()

    if not genes:
        log.warning("No genes in DB — skipping functional evidence step.")
        return {"phenotype": phenotype, "n_genes": 0, "n_human": 0,
                "n_ot": 0, "n_gtex": 0, "n_depmap": 0, "n_scored": 0}

    human = _human_genes(genes)
    gene_ids = [g.id for g in genes]
    log.info(
        "Functional evidence (step 8): %d total genes, %d human genes, phenotype=%r",
        len(genes), len(human), phenotype,
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
        "n_human": len(human),
        "n_ot": len(ot_scores),
        "n_gtex": len(gtex_scores),
        "n_depmap": len(depmap_scores),
        "n_scored": nonzero,
    }
    log.info("Functional evidence complete: %s", summary)
    return summary
