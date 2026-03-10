"""Step 8 — Expression annotation via GEO + DESeq2.

For species with available GEO datasets (ground squirrel, naked mole rat, axolotl priority):
  1. Search GEO for datasets matching the species + protective condition.
  2. Download raw count matrix via GEOparse.
  3. Run DESeq2 differential expression analysis via rpy2.
  4. Flag candidate genes as upregulated (log2FC > 1, padj < 0.05).
  5. Store expression_score in CandidateScore.
"""

import logging
import tempfile
from pathlib import Path
from typing import Optional

import pandas as pd

from db.models import CandidateScore, ExpressionResult, Gene, Ortholog
from db.session import get_session
from pipeline.config import get_ncbi_api_key, get_ncbi_email, get_storage_root, get_thresholds

log = logging.getLogger(__name__)

# Species known to have useful GEO datasets in Phase 1
GEO_PRIORITY_SPECIES = {"ground_squirrel", "naked_mole_rat", "axolotl", "spiny_mouse", "mouse_lemur"}


def search_geo(species_id: str, search_terms: list[str]) -> list[str]:
    """Search GEO for dataset accessions matching the species + conditions.

    Returns a list of GSE accessions (e.g. ["GSE12345", ...]).
    """
    try:
        import requests

        api_key = get_ncbi_api_key()
        email = get_ncbi_email()

        query = " OR ".join(f'"{t}"' for t in search_terms)
        query += ' AND "expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type]'

        params = {
            "db": "gds",
            "term": query,
            "retmax": 10,
            "retmode": "json",
            "usehistory": "n",
        }
        if api_key:
            params["api_key"] = api_key
        if email:
            params["email"] = email

        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params=params,
            timeout=30,
        )
        r.raise_for_status()
        data = r.json()
        ids = data.get("esearchresult", {}).get("idlist", [])

        if not ids:
            log.info("  GEO: no datasets found for %s", species_id)
            return []

        # Convert GDS IDs to GSE accessions
        gse_accessions = _gds_ids_to_gse(ids)
        log.info("  GEO: found %d datasets for %s: %s", len(gse_accessions), species_id, gse_accessions)
        return gse_accessions

    except Exception as exc:
        log.warning("  GEO search failed for %s: %s", species_id, exc)
        return []


def _gds_ids_to_gse(gds_ids: list[str]) -> list[str]:
    """Convert GDS numeric IDs to GSE accession strings via Entrez summary."""
    try:
        import requests

        params = {
            "db": "gds",
            "id": ",".join(gds_ids),
            "retmode": "json",
        }
        r = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params=params,
            timeout=30,
        )
        r.raise_for_status()
        data = r.json()
        accessions = []
        for uid, summary in data.get("result", {}).items():
            if uid == "uids":
                continue
            acc = summary.get("accession", "")
            if acc.startswith("GSE"):
                accessions.append(acc)
        return accessions
    except Exception:
        return []


def download_geo_dataset(gse_accession: str, out_dir: Path) -> Optional[Path]:
    """Download a GEO dataset using GEOparse.

    Returns path to directory containing the soft file, or None on failure.
    """
    try:
        import GEOparse

        geo_dir = out_dir / gse_accession
        geo_dir.mkdir(parents=True, exist_ok=True)

        gse = GEOparse.get_GEO(geo=gse_accession, destdir=str(geo_dir), silent=True)
        return geo_dir, gse

    except Exception as exc:
        log.warning("  GEOparse failed for %s: %s", gse_accession, exc)
        return None, None


def extract_count_matrix(gse) -> Optional[pd.DataFrame]:
    """Extract a gene × sample count matrix from a GEOparse GSE object.

    Handles both microarray (intensity) and RNA-seq (counts) data.
    Returns None if no suitable data is found.
    """
    try:
        tables = []
        for gsm_name, gsm in gse.gsms.items():
            if gsm.table is not None and not gsm.table.empty:
                col = gsm.table.set_index(gsm.table.columns[0])[gsm.table.columns[-1]]
                col.name = gsm_name
                tables.append(col)

        if not tables:
            return None

        matrix = pd.concat(tables, axis=1).fillna(0)
        # Keep only integer-like columns (RNA-seq counts)
        matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0)
        return matrix

    except Exception as exc:
        log.warning("  Count matrix extraction failed: %s", exc)
        return None


def run_deseq2(count_matrix: pd.DataFrame, condition_map: dict[str, str]) -> Optional[pd.DataFrame]:
    """Run DESeq2 differential expression analysis via rpy2.

    Args:
        count_matrix: Gene × sample DataFrame with integer counts.
        condition_map: {sample_name: "treatment" | "control"}

    Returns:
        DataFrame with columns [gene_id, log2FoldChange, pvalue, padj] or None.
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr

        pandas2ri.activate()
        deseq2 = importr("DESeq2")
        base = importr("base")

        # Filter matrix to samples with condition labels
        samples = [s for s in count_matrix.columns if s in condition_map]
        matrix = count_matrix[samples].astype(int)

        conditions = [condition_map[s] for s in samples]
        col_data = pd.DataFrame({"condition": conditions}, index=samples)

        r_matrix = pandas2ri.py2rpy(matrix)
        r_col_data = pandas2ri.py2rpy(col_data)

        dds = deseq2.DESeqDataSetFromMatrix(
            countData=r_matrix,
            colData=r_col_data,
            design=ro.Formula("~ condition"),
        )
        dds = deseq2.DESeq(dds)
        res = deseq2.results(dds, contrast=ro.StrVector(["condition", "treatment", "control"]))
        res_df = pandas2ri.rpy2py(base.as_data_frame(res))
        res_df.index.name = "gene_id"
        res_df = res_df.reset_index()
        return res_df

    except ImportError:
        log.warning("rpy2 or DESeq2 not available — skipping DESeq2 analysis.")
        return None
    except Exception as exc:
        log.warning("DESeq2 analysis failed: %s", exc)
        return None


def compute_expression_score(
    deseq2_results: pd.DataFrame,
    candidate_gene_symbols: set[str],
) -> dict[str, float]:
    """Convert DESeq2 results to expression scores for candidate genes.

    A gene is upregulated if log2FC > threshold AND padj < threshold.
    Score = 1.0 if upregulated, 0.5 if detected but not DE, 0.0 if not detected.

    Returns {gene_symbol: expression_score}
    """
    thresholds = get_thresholds()
    log2fc_min = thresholds.get("expression_log2fc_min", 1.0)
    padj_max = thresholds.get("expression_padj_max", 0.05)

    scores = {}
    for _, row in deseq2_results.iterrows():
        gid = str(row.get("gene_id", ""))
        if gid not in candidate_gene_symbols:
            continue

        log2fc = row.get("log2FoldChange", 0)
        padj = row.get("padj", 1.0)

        if pd.isna(log2fc) or pd.isna(padj):
            scores[gid] = 0.5
        elif log2fc >= log2fc_min and padj <= padj_max:
            scores[gid] = 1.0
        else:
            scores[gid] = 0.5

    return scores


def save_expression_evidence(
    results: pd.DataFrame,
    geo_accession: str,
    comparison: str,
    gene_symbol_to_id: dict[str, str],
) -> int:
    """Write per-gene expression evidence to ExpressionResult for traceability."""
    if results is None or results.empty:
        return 0
    n = 0
    with get_session() as session:
        for _, row in results.iterrows():
            symbol = str(row.get("gene_id", "")).strip()
            gene_id = gene_symbol_to_id.get(symbol)
            if not gene_id:
                continue
            log2fc = row.get("log2FoldChange")
            padj = row.get("padj")
            if pd.isna(log2fc):
                log2fc = None
            if pd.isna(padj):
                padj = None
            session.add(ExpressionResult(
                gene_id=gene_id,
                geo_accession=geo_accession,
                comparison=comparison,
                log2fc=float(log2fc) if log2fc is not None else None,
                padj=float(padj) if padj is not None else None,
                n_samples=int(results.shape[1]) - 2 if results.shape[1] > 2 else None,
            ))
            n += 1
        session.commit()
    return n


def save_expression_scores(scores_by_species: dict[str, dict[str, float]]) -> None:
    """Aggregate expression scores across species and update CandidateScore."""
    with get_session() as session:
        genes = session.query(Gene).all()
        gene_map = {g.gene_symbol: g.id for g in genes}

        for gene_symbol, gene_id in gene_map.items():
            species_scores = []
            for species_id, scores in scores_by_species.items():
                if gene_symbol in scores:
                    species_scores.append(scores[gene_symbol])

            if not species_scores:
                continue

            # Average expression score across all species with data
            expr_score = sum(species_scores) / len(species_scores)

            cs = session.query(CandidateScore).filter_by(gene_id=gene_id, trait_id="").first()
            if cs is None:
                cs = CandidateScore(gene_id=gene_id, trait_id="")
                session.add(cs)
            cs.expression_score = round(expr_score, 4)

    log.info("Expression scores saved for %d genes.", len(gene_map))


def run_expression_pipeline(species_list: list[dict]) -> dict[str, dict[str, float]]:
    """Run the full expression pipeline for all priority species.

    Returns {species_id: {gene_symbol: expression_score}}
    """
    root = Path(get_storage_root())
    geo_dir = root / "geo_data"
    geo_dir.mkdir(parents=True, exist_ok=True)

    scores_by_species: dict[str, dict[str, float]] = {}

    for entry in species_list:
        sid = entry["id"]
        if sid == "human" or sid not in GEO_PRIORITY_SPECIES:
            continue
        if not entry.get("geo_search_terms"):
            continue

        log.info("Processing GEO expression for %s...", sid)
        accessions = search_geo(sid, entry["geo_search_terms"])

        for gse_acc in accessions[:2]:   # Limit to 2 datasets per species in Phase 1
            geo_path, gse = download_geo_dataset(gse_acc, geo_dir)
            if gse is None:
                continue

            matrix = extract_count_matrix(gse)
            if matrix is None or matrix.empty:
                continue

            # Infer condition labels from GSM titles (heuristic)
            condition_map = _infer_conditions(gse)
            if not condition_map or len(set(condition_map.values())) < 2:
                log.info("  %s/%s: could not infer treatment/control groups", sid, gse_acc)
                continue

            results = run_deseq2(matrix, condition_map)
            if results is None:
                continue

            with get_session() as session:
                genes = session.query(Gene).all()
                gene_symbols = {g.gene_symbol for g in genes}
                gene_symbol_to_id = {g.gene_symbol: g.id for g in genes}

            species_scores = compute_expression_score(results, gene_symbols)
            scores_by_species[sid] = species_scores
            n_ev = save_expression_evidence(results, gse_acc, sid, gene_symbol_to_id)
            log.info("  %s/%s: scored %d genes, saved %d evidence rows", sid, gse_acc, len(species_scores), n_ev)

    return scores_by_species


def _infer_conditions(gse) -> dict[str, str]:
    """Heuristic: assign treatment/control based on GSM title keywords."""
    condition_map = {}
    treatment_keywords = {
        "hibernat", "torpor", "cancer", "tumor", "wound", "injury",
        "regenerat", "cold", "stress", "treatment", "treated",
    }
    control_keywords = {
        "control", "normal", "baseline", "untreated", "active", "summer",
    }

    for gsm_name, gsm in gse.gsms.items():
        title = (gsm.metadata.get("title", [""])[0] or "").lower()
        if any(k in title for k in treatment_keywords):
            condition_map[gsm_name] = "treatment"
        elif any(k in title for k in control_keywords):
            condition_map[gsm_name] = "control"

    return condition_map
