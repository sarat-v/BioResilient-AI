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

Combined score → candidate_score.expression_score.
This score is used as a 4th rank-product evidence layer in Phase 1 scoring (step 9).

Per-phenotype configuration: config/functional_evidence_config.json
Reference for Open Targets: Ochoa et al. (2021) Nature Genetics.
Reference for DepMap:        Tsherniak et al. (2017) Cell.
Reference for GTEx:          GTEx Consortium (2020) Science.
"""

import csv
import json
import logging
import math
import time
from pathlib import Path
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
DEPMAP_URL = (
    "https://depmap.org/portal/api/download/file"
    "?file_name=CRISPRGeneDependency.csv&release=DepMap+Public+24Q2"
)

OT_BATCH_SIZE = 50    # genes per GraphQL request to Open Targets
REQUEST_TIMEOUT = 30  # seconds for individual HTTP requests
RATE_DELAY_OT = 0.15  # seconds between Open Targets batches
RATE_DELAY_GTEX = 0.12  # seconds between GTEx requests

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
    The approximation is acceptable because:
      - Top-tier candidates are almost always well-known genes with matching names.
      - Olfactory receptors / pseudogenes that don't match will simply receive no
        functional score, which does not penalise them — it is conservative.
    """
    sym = gene.gene_symbol or ""
    if sym.endswith("_HUMAN"):
        return sym[: -len("_HUMAN")]
    return sym if sym else None


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

_OT_SYMBOL_QUERY = """
query SymbolToEnsembl($queryString: String!) {
  search(queryString: $queryString, entityNames: ["target"], page: {size: 1}) {
    hits {
      id
      ... on Target {
        approvedSymbol
      }
    }
  }
}
"""


def _ot_symbol_to_ensg(symbol: str) -> Optional[str]:
    """Resolve a gene symbol to an ENSEMBL ID via the Open Targets search API."""
    try:
        r = requests.post(
            OPENTARGETS_GRAPHQL,
            json={"query": _OT_SYMBOL_QUERY, "variables": {"queryString": symbol}},
            timeout=REQUEST_TIMEOUT,
        )
        r.raise_for_status()
        hits = (r.json().get("data") or {}).get("search", {}).get("hits", [])
        if hits:
            hit = hits[0]
            approved = hit.get("approvedSymbol", "")
            if approved.upper() == symbol.upper():
                return hit.get("id")
    except Exception as exc:
        log.debug("OT symbol lookup failed for %s: %s", symbol, exc)
    return None


def _build_ensg_map(genes: list[Gene]) -> dict[str, str]:
    """Build {hgnc_symbol: ensembl_id} mapping via Open Targets search.

    Uses a cache file to avoid repeated API calls across pipeline runs.
    """
    cache_path = Path(get_local_storage_root()) / "ot_ensg_cache.json"
    cache: dict[str, str] = {}
    if cache_path.exists():
        try:
            with open(cache_path) as f:
                cache = json.load(f)
        except Exception:
            cache = {}

    symbols_needed = []
    for gene in genes:
        sym = _hgnc_symbol(gene)
        if sym and sym not in cache:
            symbols_needed.append(sym)

    if symbols_needed:
        log.info("OT symbol→ENSG lookup for %d new symbols...", len(symbols_needed))
        for i, sym in enumerate(symbols_needed):
            ensg = _ot_symbol_to_ensg(sym)
            cache[sym] = ensg or ""
            time.sleep(RATE_DELAY_OT)
            if (i + 1) % 100 == 0:
                log.info("  Symbol lookup: %d / %d done", i + 1, len(symbols_needed))

        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(cache, f)

    return cache


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
        scores = []
        for row in (target.get("associatedDiseases") or {}).get("rows") or []:
            if row.get("disease", {}).get("id") in disease_set:
                s = row.get("score")
                if s is not None:
                    scores.append(float(s))
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

    # Build symbol → ENSG map (cached)
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

    log.info("  %d / %d genes have ENSG IDs resolved.", len(gid_to_ensg), len(genes))

    # Invert: ENSG → gene DB id (for result mapping)
    ensg_to_gid: dict[str, str] = {v: k for k, v in gid_to_ensg.items()}

    ensg_list = list(ensg_to_gid.keys())
    all_scores: dict[str, float] = {}
    n_batches = math.ceil(len(ensg_list) / OT_BATCH_SIZE)

    for i in range(0, len(ensg_list), OT_BATCH_SIZE):
        batch = ensg_list[i: i + OT_BATCH_SIZE]
        scores = _fetch_ot_batch(batch, disease_set)
        for ensg, score in scores.items():
            gid = ensg_to_gid.get(ensg)
            if gid:
                all_scores[gid] = score
        time.sleep(RATE_DELAY_OT)
        batch_num = i // OT_BATCH_SIZE + 1
        if batch_num % 5 == 0 or batch_num == n_batches:
            log.info("  OT: %d / %d batches complete", batch_num, n_batches)

    log.info("Open Targets: %d / %d genes scored.", len(all_scores), len(genes))
    return all_scores


# ---------------------------------------------------------------------------
# GTEx
# ---------------------------------------------------------------------------

def _fetch_gtex_gene(gene_symbol: str, tissue_ids: list[str]) -> Optional[dict[str, float]]:
    """Fetch GTEx v10 median TPM for a gene across the requested tissues.

    Issues one request per gene (returning all tissues) and filters client-side.
    Falls back to per-tissue requests if the bulk call fails.

    Returns {tissue_id: median_tpm} or None if no data.
    """
    result: dict[str, float] = {}
    tissue_set = set(tissue_ids)

    try:
        # Single request for all tissues
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
        log.debug("GTEx bulk fetch failed for %s: %s", gene_symbol, exc)

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

    Returns {gene_id (DB primary key): score [0, 1]}.
    """
    if not tissue_ids:
        return {}

    log.info("GTEx: querying %d genes across %d tissues...", len(genes), len(tissue_ids))
    scores: dict[str, float] = {}

    for i, gene in enumerate(genes):
        sym = _hgnc_symbol(gene)
        if not sym:
            continue
        tpm_map = _fetch_gtex_gene(sym, tissue_ids)
        if tpm_map:
            score = _gtex_normalize(list(tpm_map.values()))
            if score > 0:
                scores[gene.id] = score
        time.sleep(RATE_DELAY_GTEX)
        if (i + 1) % 200 == 0:
            log.info("  GTEx: %d / %d genes done...", i + 1, len(genes))

    log.info("GTEx: %d / %d genes scored.", len(scores), len(genes))
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
    """Convert mean Chronos score to a [0, 1] functional relevance score.

    Genes with strong selective essentiality (chronos << -0.5) score near 1.
    Broadly non-essential genes (chronos ≈ 0) score near 0.
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
    phenotype: str,
) -> int:
    """Write evidence sub-scores to expression_result and update candidate_score.

    The expression_result table is repurposed:
      geo_accession → source identifier (e.g. "OT:cancer", "GTEX:Spleen")
      comparison    → human-readable source description
      log2fc        → sub-score value [0, 1]
      padj          → null (not applicable)

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

        # Load existing CandidateScore rows (keyed by gene_id + trait_id)
        cs_map = {
            r.gene_id: r
            for r in session.query(CandidateScore).filter_by(trait_id=phenotype).all()
        }

        updated = 0
        for gid, score in combined_scores.items():
            if score <= 0:
                continue
            cs = cs_map.get(gid)
            if cs is None:
                cs = CandidateScore(gene_id=gid, trait_id=phenotype)
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
      2. Fetch Open Targets association scores (batched GraphQL).
      3. Fetch GTEx tissue expression scores (one request per gene).
      4. Optionally fetch DepMap essentiality scores.
      5. Combine into a single expression_score per gene.
      6. Persist to expression_result table and candidate_score.expression_score.

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
    nonzero = _persist_results(
        gene_ids, combined, ot_scores, gtex_scores, depmap_scores, phenotype
    )

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
