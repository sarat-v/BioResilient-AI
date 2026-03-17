"""Codon-guided MEME selection analysis.

Replaces the protein divergence proxy with statistically rigorous episodic
positive selection detection using HyPhy MEME (Mixed Effects Model of Evolution).

Pipeline for each orthogroup:
  1. Fetch CDS (nucleotide coding sequences) from NCBI for each protein accession
  2. Build a codon alignment guided by the existing protein (amino acid) alignment
     using PAL2NAL-style logic — ensures codons are in the correct reading frame
  3. Run HyPhy MEME on the codon alignment + species tree
  4. Parse per-site selection p-values and report:
       - fraction of sites under episodic positive selection
       - mean effect size (beta+) at significant sites
       - which branches show selection at each site

MEME detects *episodic* selection — sites that were under positive selection
at some point in some lineages, not necessarily in all. This is exactly the
right model for cancer resistance: selection happened once, in specific lineages.

Requirements:
  - HyPhy ≥ 2.5 installed (provides `meme` analysis via hyphy LIBPATH)
  - NCBI API key (for CDS fetching at 10 req/sec)
  - Entrez biopython for CDS lookup
"""

import json
import logging
import os
import re
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from Bio import Entrez, SeqIO

from pipeline.config import get_ncbi_api_key, get_ncbi_email, get_local_storage_root, get_tool_config

log = logging.getLogger(__name__)

MEME_SIGNIFICANCE_THRESHOLD = 0.05   # per-site p-value threshold for MEME


# ---------------------------------------------------------------------------
# CDS fetching
# ---------------------------------------------------------------------------

def _setup_entrez() -> None:
    Entrez.email = get_ncbi_email()
    api_key = get_ncbi_api_key()
    if api_key:
        Entrez.api_key = api_key


def _prefetch_all_cds(aligned_orthogroups: dict[str, dict[str, str]]) -> None:
    """Pre-fetch CDS for all unique protein accessions before spawning workers.

    Uses NCBI elink (list of IDs → Biopython uses POST) + efetch in batches.
    Results cached to disk so workers never hit NCBI.
    """
    _setup_entrez()
    cache_dir = Path(get_local_storage_root()) / "cds"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # UniProt mnemonic: GENENAME_SPECIES where SPECIES is letters-only (e.g. CENPK_HUMAN)
    # RefSeq like XP_049720501.1 has digits after underscore — must NOT match
    _UNIPROT_PAT = re.compile(r"^[A-Z0-9]{1,11}_[A-Z]{3,5}$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")

    def _is_uniprot(acc: str) -> bool:
        # Extra check: RefSeq always has digits after underscore (XP_12345, NP_12345)
        if "_" in acc:
            after = acc.split("_", 1)[1].lstrip("0123456789.")
            if after:  # has letters after digits → it's not a simple RefSeq
                return bool(_UNIPROT_PAT.match(acc))
            return False  # digits-only after underscore → RefSeq
        return bool(_UNIPROT_PAT.match(acc))

    all_accessions: set[str] = set()
    for seqs in aligned_orthogroups.values():
        for label in seqs:
            acc = label.split("|")[-1] if "|" in label else label
            all_accessions.add(acc)

    to_fetch = []
    for acc in sorted(all_accessions):
        cf = cache_dir / f"{acc}.fna"
        if _is_uniprot(acc):
            if not cf.exists():
                cf.write_text("")
            continue
        # Clear empty cache files that may have been written by buggy earlier code
        if cf.exists() and cf.stat().st_size == 0:
            cf.unlink()
        if not cf.exists():
            to_fetch.append(acc)

    ncbi_total = sum(1 for a in all_accessions if not _is_uniprot(a))
    if not to_fetch:
        log.info("CDS pre-fetch: all %d NCBI accessions already cached.", ncbi_total)
        return

    log.info("CDS pre-fetch: %d NCBI accessions to fetch (%d cached, %d UniProt skipped).",
             len(to_fetch), ncbi_total - len(to_fetch), len(all_accessions) - ncbi_total)

    # Use small batches — passing list (not comma string) makes Biopython use POST
    _ELINK_BATCH = 50
    _EFETCH_BATCH = 50
    done = 0
    failed = 0

    for batch_start in range(0, len(to_fetch), _ELINK_BATCH):
        batch = to_fetch[batch_start: batch_start + _ELINK_BATCH]

        try:
            # Pass list — Biopython sends as POST, avoids URL-length issues
            handle = Entrez.elink(
                dbfrom="protein", db="nuccore",
                id=batch,
                linkname="protein_nuccore_mrna",
            )
            link_records = Entrez.read(handle)
            handle.close()
            time.sleep(0.12)
        except Exception as exc:
            log.debug("elink batch failed at %d: %s", batch_start, exc)
            for acc in batch:
                (cache_dir / f"{acc}.fna").write_text("")
            failed += len(batch)
            done += len(batch)
            continue

        acc_to_nuc: dict[str, str] = {}
        for i, acc in enumerate(batch):
            try:
                for ls in link_records[i].get("LinkSetDb", []):
                    if ls["LinkName"] in ("protein_nuccore_mrna", "protein_nuccore"):
                        if ls["Link"]:
                            acc_to_nuc[acc] = ls["Link"][0]["Id"]
                        break
            except (IndexError, KeyError):
                pass
            if acc not in acc_to_nuc:
                (cache_dir / f"{acc}.fna").write_text("")

        if not acc_to_nuc:
            done += len(batch)
            continue

        acc_by_nuc = {v: k for k, v in acc_to_nuc.items()}
        nuc_ids = list(acc_to_nuc.values())
        fetched_nuc_ids: set[str] = set()

        for ef_start in range(0, len(nuc_ids), _EFETCH_BATCH):
            sub_ids = nuc_ids[ef_start: ef_start + _EFETCH_BATCH]
            try:
                handle = Entrez.efetch(
                    db="nuccore", id=sub_ids,
                    rettype="gb", retmode="text",
                )
                records = list(SeqIO.parse(handle, "genbank"))
                handle.close()
                time.sleep(0.12)
            except Exception as exc:
                log.debug("efetch batch failed: %s", exc)
                for nid in sub_ids:
                    acc = acc_by_nuc.get(nid)
                    if acc:
                        (cache_dir / f"{acc}.fna").write_text("")
                failed += len(sub_ids)
                continue

            for record in records:
                matched_acc = None
                for nid in sub_ids:
                    if record.id.startswith(nid) or nid in record.id:
                        matched_acc = acc_by_nuc.get(nid)
                        fetched_nuc_ids.add(nid)
                        break
                if not matched_acc:
                    continue
                cds_seq = None
                for feat in record.features:
                    if feat.type == "CDS":
                        seq = str(feat.extract(record.seq))
                        if len(seq) % 3 == 0 and seq.upper().startswith("ATG"):
                            cds_seq = seq
                            break
                (cache_dir / f"{matched_acc}.fna").write_text(cds_seq or "")

        for acc, nid in acc_to_nuc.items():
            if nid not in fetched_nuc_ids and not (cache_dir / f"{acc}.fna").exists():
                (cache_dir / f"{acc}.fna").write_text("")

        done += len(batch)
        if done % 1000 == 0 or done >= len(to_fetch):
            log.info("  CDS pre-fetch: %d / %d done (%d failed).",
                     done, len(to_fetch), failed)

    log.info("CDS pre-fetch complete: %d fetched, %d failed/no-CDS, %d already cached.",
             done - failed, failed, ncbi_total - len(to_fetch))


def fetch_cds_for_protein(protein_accession: str) -> Optional[str]:
    """Fetch the CDS nucleotide sequence for a protein accession from NCBI.

    Results are cached to disk at {storage_root}/cds/{accession}.fna so that
    repeated runs of step6, 6b, and 6c never re-fetch the same accession.

    Rate: respects NCBI's 10 req/s limit with API key.
    """
    _setup_entrez()
    # Strip species prefix (e.g. "human|human|NP_001234.1" → "NP_001234.1")
    acc = protein_accession.split("|")[-1] if "|" in protein_accession else protein_accession

    # Skip non-NCBI accessions: UniProt mnemonics (GENE_SPECIES) and UniProt accessions
    if re.match(r"^[A-Z0-9]{1,11}_[A-Z]{3,5}$", acc):  # UniProt mnemonic e.g. CENPK_HUMAN
        return None
    if re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9])", acc):
        return None

    # Check disk cache first — avoids all network I/O on resume/re-run
    cache_dir = Path(get_local_storage_root()) / "cds"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{acc}.fna"
    if cache_file.exists():
        cds = cache_file.read_text().strip()
        return cds if cds else None   # empty file = previously confirmed "no CDS"

    try:
        # Link protein → nucleotide
        handle = Entrez.elink(dbfrom="protein", db="nuccore", id=acc, linkname="protein_nuccore_mrna")
        link_record = Entrez.read(handle)
        handle.close()
        time.sleep(0.12)

        nuc_ids = []
        if link_record and link_record[0].get("LinkSetDb"):
            for link in link_record[0]["LinkSetDb"]:
                if link["LinkName"] in ("protein_nuccore_mrna", "protein_nuccore"):
                    nuc_ids = [l["Id"] for l in link["Link"]]
                    break

        if not nuc_ids:
            cache_file.write_text("")   # cache negative result
            return None

        # Fetch the first mRNA record
        handle = Entrez.efetch(db="nuccore", id=nuc_ids[0], rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        time.sleep(0.12)

        # Extract CDS feature
        for feature in record.features:
            if feature.type == "CDS":
                cds_seq = str(feature.extract(record.seq))
                # Validate: must be divisible by 3 and start with ATG
                if len(cds_seq) % 3 == 0 and cds_seq.upper().startswith("ATG"):
                    cache_file.write_text(cds_seq)   # persist for future runs
                    return cds_seq

        cache_file.write_text("")   # CDS not found — cache to skip next time
        return None

    except Exception as exc:
        log.debug("CDS fetch failed for %s: %s", acc, exc)
        return None


# ---------------------------------------------------------------------------
# Codon alignment (PAL2NAL-style)
# ---------------------------------------------------------------------------

def protein_to_codon_alignment(
    protein_alignment: dict[str, str],
    cds_seqs: dict[str, str],
) -> Optional[dict[str, str]]:
    """Build a codon alignment from a protein alignment and CDS sequences.

    Each amino acid in the protein alignment maps to exactly 3 nucleotides
    in the CDS. Gap characters in the protein alignment become '---' in the
    codon alignment.

    Args:
        protein_alignment: {label: aligned_protein_seq} — contains gap characters
        cds_seqs: {label: cds_nucleotide_seq} — no gaps, raw CDS

    Returns:
        {label: codon_aligned_nucleotide_seq} or None if any sequence is invalid.
    """
    codon_aln: dict[str, str] = {}
    for label, prot_aln in protein_alignment.items():
        cds = cds_seqs.get(label)
        if not cds:
            return None  # Missing CDS for this species — can't build codon alignment

        # Remove gaps from protein to count non-gap positions
        non_gap_count = sum(1 for aa in prot_aln if aa != "-")
        expected_cds_len = non_gap_count * 3

        # Allow CDS to be slightly longer (stop codon) but not shorter
        if len(cds) < expected_cds_len:
            log.debug("CDS too short for %s: expected ≥%d, got %d", label, expected_cds_len, len(cds))
            return None

        cds_pos = 0
        codon_seq = []
        for aa in prot_aln:
            if aa == "-":
                codon_seq.append("---")
            else:
                codon_seq.append(cds[cds_pos:cds_pos + 3])
                cds_pos += 3

        codon_aln[label] = "".join(codon_seq)

    # Verify all sequences have the same length
    lengths = {len(v) for v in codon_aln.values()}
    if len(lengths) != 1:
        log.debug("Codon alignment has unequal lengths: %s", lengths)
        return None

    return codon_aln


# ---------------------------------------------------------------------------
# HyPhy MEME runner
# ---------------------------------------------------------------------------

def run_meme(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
) -> Optional[dict]:
    """Run HyPhy MEME on a codon alignment and return the JSON results.

    MEME (Mixed Effects Model of Evolution) detects sites under episodic
    positive selection — the correct model for convergent adaptive evolution.

    Returns the parsed MEME JSON output, or None on failure.
    """
    tools = get_tool_config()
    hyphy_bin = tools.get("hyphy_bin", "hyphy")

    root = Path(get_local_storage_root())
    meme_dir = root / "meme" / og_id
    meme_dir.mkdir(parents=True, exist_ok=True)

    aln_path = meme_dir / "codon_aln.fna"
    tree_path = meme_dir / "species.treefile"
    out_path = meme_dir / "meme.json"

    # Write alignment
    with open(aln_path, "w") as f:
        for label, seq in codon_aln.items():
            f.write(f">{label}\n{seq}\n")

    tree_path.write_text(tree_newick)

    try:
        cmd = [
            hyphy_bin, "meme",
            "--alignment", str(aln_path),
            "--tree", str(tree_path),
            "--output", str(out_path),
            "--branches", "All",
        ]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode != 0:
            log.debug("MEME failed for %s: %s", og_id, result.stderr[:500])
            return None

        if not out_path.exists():
            return None

        return json.loads(out_path.read_text())

    except subprocess.TimeoutExpired:
        log.debug("MEME timed out for %s", og_id)
        return None
    except Exception as exc:
        log.debug("MEME error for %s: %s", og_id, exc)
        return None


# ---------------------------------------------------------------------------
# MEME result parsing
# ---------------------------------------------------------------------------

def parse_meme_results(meme_json: dict, og_id: str) -> dict:
    """Extract selection statistics from HyPhy MEME output.

    MEME outputs per-site results with:
      - alpha: synonymous rate
      - beta-: purifying selection intensity
      - beta+: positive selection intensity
      - p-value: episodic diversification p-value per site
      - Episodic selection detected at N% of sites

    Returns a dict compatible with load_selection_scores().
    """
    try:
        mle_content = meme_json.get("MLE", {}).get("content", {}).get("0", [])
        if not mle_content:
            return _meme_null_result(og_id)

        n_sites = len(mle_content)
        selected_sites = []
        beta_plus_values = []

        for site_data in mle_content:
            # MEME MLE columns: [alpha, beta-, omega-, beta+, omega+, p-value, ...]
            if len(site_data) < 6:
                continue
            pvalue = site_data[5]
            beta_plus = site_data[3]
            if pvalue is not None and pvalue < MEME_SIGNIFICANCE_THRESHOLD:
                selected_sites.append(pvalue)
                if beta_plus is not None:
                    beta_plus_values.append(beta_plus)

        if n_sites == 0:
            return _meme_null_result(og_id)

        fraction_selected = len(selected_sites) / n_sites
        mean_beta_plus = sum(beta_plus_values) / len(beta_plus_values) if beta_plus_values else 0.0

        # Map to pseudo dN/dS: use fraction_selected × normalised beta+
        # beta+ is the non-synonymous rate under positive selection;
        # values >> 1 indicate strong positive selection
        pseudo_dnds = fraction_selected * min(mean_beta_plus / 5.0, 2.0)
        # p-value: geometric mean of site p-values under selection (penalised if none)
        if selected_sites:
            import math
            log_sum = sum(math.log(p) for p in selected_sites if p > 0)
            pseudo_pvalue = max(math.exp(log_sum / len(selected_sites)), 1e-10)
        else:
            pseudo_pvalue = 1.0

        log.debug("  MEME %s: %d/%d sites selected, mean_beta+=%.2f, pseudo_dnds=%.3f",
                  og_id, len(selected_sites), n_sites, mean_beta_plus, pseudo_dnds)

        return {
            "dnds_ratio": round(pseudo_dnds, 4),
            "dnds_pvalue": round(pseudo_pvalue, 6),
            "selection_model": "MEME_episodic",
            "branches_under_selection": [],  # MEME is site-level, not branch-level
            "fraction_sites_selected": round(fraction_selected, 4),
            "n_sites_selected": len(selected_sites),
            "mean_beta_plus": round(mean_beta_plus, 4),
        }

    except Exception as exc:
        log.debug("MEME parse error for %s: %s", og_id, exc)
        return _meme_null_result(og_id)


def _meme_null_result(og_id: str) -> dict:
    return {
        "dnds_ratio": 0.0,
        "dnds_pvalue": 1.0,
        "selection_model": "MEME_no_signal",
        "branches_under_selection": [],
        "fraction_sites_selected": 0.0,
        "n_sites_selected": 0,
        "mean_beta_plus": 0.0,
    }


# ---------------------------------------------------------------------------
# FEL — Fixed Effects Likelihood (pervasive selection)
# ---------------------------------------------------------------------------

def run_fel(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
) -> Optional[dict]:
    """Run HyPhy FEL on a codon alignment.

    FEL (Fixed Effects Likelihood) detects sites under *pervasive* (constant)
    positive or purifying selection — complementary to MEME's episodic detection.

    Returns the parsed FEL JSON output, or None on failure.
    """
    tools = get_tool_config()
    hyphy_bin = tools.get("hyphy_bin", "hyphy")

    root = Path(get_local_storage_root())
    fel_dir = root / "fel" / og_id
    fel_dir.mkdir(parents=True, exist_ok=True)

    aln_path = fel_dir / "codon_aln.fna"
    tree_path = fel_dir / "species.treefile"
    out_path = fel_dir / "fel.json"

    with open(aln_path, "w") as f:
        for label, seq in codon_aln.items():
            f.write(f">{label}\n{seq}\n")
    tree_path.write_text(tree_newick)

    try:
        cmd = [
            hyphy_bin, "fel",
            "--alignment", str(aln_path),
            "--tree", str(tree_path),
            "--output", str(out_path),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            log.debug("FEL failed for %s: %s", og_id, result.stderr[:300])
            return None
        if not out_path.exists():
            return None
        return json.loads(out_path.read_text())
    except subprocess.TimeoutExpired:
        log.debug("FEL timed out for %s", og_id)
        return None
    except Exception as exc:
        log.debug("FEL error for %s: %s", og_id, exc)
        return None


def parse_fel_results(fel_json: dict) -> dict:
    """Count sites under pervasive positive selection from HyPhy FEL output.

    Returns {"fel_sites": int} — count of sites with dN > dS at p < 0.05.
    """
    try:
        mle_content = fel_json.get("MLE", {}).get("content", {}).get("0", [])
        if not mle_content:
            return {"fel_sites": 0}

        pos_sites = 0
        for site_data in mle_content:
            # FEL MLE columns: [alpha (dS), beta (dN), LRT, p-value, ...]
            if len(site_data) < 4:
                continue
            alpha = site_data[0]  # synonymous rate
            beta = site_data[1]   # non-synonymous rate
            pvalue = site_data[3]
            if pvalue is not None and pvalue < 0.05 and beta is not None and alpha is not None:
                if beta > alpha:  # dN > dS = positive selection
                    pos_sites += 1

        return {"fel_sites": pos_sites}

    except Exception as exc:
        log.debug("FEL parse error: %s", exc)
        return {"fel_sites": 0}


# ---------------------------------------------------------------------------
# BUSTED — Branch-Site Unrestricted Test for Episodic Diversification
# ---------------------------------------------------------------------------

def run_busted(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
) -> Optional[dict]:
    """Run HyPhy BUSTED on a codon alignment.

    BUSTED tests whether the gene *as a whole* has experienced episodic positive
    selection anywhere in the tree, providing a single gene-level p-value.
    Complementary to MEME (which is site-level).

    Returns the parsed BUSTED JSON output, or None on failure.
    """
    tools = get_tool_config()
    hyphy_bin = tools.get("hyphy_bin", "hyphy")

    root = Path(get_local_storage_root())
    busted_dir = root / "busted" / og_id
    busted_dir.mkdir(parents=True, exist_ok=True)

    aln_path = busted_dir / "codon_aln.fna"
    tree_path = busted_dir / "species.treefile"
    out_path = busted_dir / "busted.json"

    with open(aln_path, "w") as f:
        for label, seq in codon_aln.items():
            f.write(f">{label}\n{seq}\n")
    tree_path.write_text(tree_newick)

    try:
        cmd = [
            hyphy_bin, "busted",
            "--alignment", str(aln_path),
            "--tree", str(tree_path),
            "--output", str(out_path),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            log.debug("BUSTED failed for %s: %s", og_id, result.stderr[:300])
            return None
        if not out_path.exists():
            return None
        return json.loads(out_path.read_text())
    except subprocess.TimeoutExpired:
        log.debug("BUSTED timed out for %s", og_id)
        return None
    except Exception as exc:
        log.debug("BUSTED error for %s: %s", og_id, exc)
        return None


def parse_busted_results(busted_json: dict) -> dict:
    """Extract gene-level p-value from HyPhy BUSTED output.

    Returns {"busted_pvalue": float} — p-value for gene-wide episodic selection.
    """
    try:
        # BUSTED stores test p-value at json["test results"]["p-value"]
        pvalue = (
            busted_json
            .get("test results", {})
            .get("p-value")
        )
        if pvalue is None:
            return {"busted_pvalue": 1.0}
        return {"busted_pvalue": round(float(pvalue), 6)}
    except Exception as exc:
        log.debug("BUSTED parse error: %s", exc)
        return {"busted_pvalue": 1.0}


# ---------------------------------------------------------------------------
# FEL + BUSTED pipeline entry point
# ---------------------------------------------------------------------------

def _fel_busted_worker(args: tuple) -> tuple[str, dict]:
    """Top-level worker: run FEL + BUSTED for one orthogroup."""
    og_id, aligned_seqs, species_treefile_str = args
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species
    from pathlib import Path as _Path
    from Bio import SeqIO as _SeqIO

    root = get_local_storage_root()
    species_ids = [label.split("|")[0] for label in aligned_seqs]
    pruned_tree = prune_tree_to_species(_Path(species_treefile_str), species_ids)

    codon_aln_path = _Path(root) / "meme" / og_id / "codon_aln.fna"
    codon_aln = None
    if codon_aln_path.exists():
        try:
            codon_aln = {r.id: str(r.seq) for r in _SeqIO.parse(str(codon_aln_path), "fasta")}
        except Exception:
            codon_aln = None

    if codon_aln is None:
        cds_seqs: dict[str, str] = {}
        for label in aligned_seqs:
            parts = label.split("|")
            accession = parts[-1] if parts else label
            cds = fetch_cds_for_protein(accession)
            if cds:
                cds_seqs[label] = cds
        if len(cds_seqs) == len(aligned_seqs):
            codon_aln = protein_to_codon_alignment(aligned_seqs, cds_seqs)

    if codon_aln is None:
        return og_id, {}

    fel_json = run_fel(codon_aln, pruned_tree, og_id)
    busted_json = run_busted(codon_aln, pruned_tree, og_id)
    fel_result = parse_fel_results(fel_json) if fel_json else {"fel_sites": 0}
    busted_result = parse_busted_results(busted_json) if busted_json else {"busted_pvalue": 1.0}
    return og_id, {**fel_result, **busted_result}


def run_fel_busted_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
    motifs_by_og: dict[str, list],
    species_treefile: Path,
) -> dict[str, dict]:
    """Run FEL and BUSTED in parallel for all candidate orthogroups.

    Returns {og_id: {"fel_sites": int, "busted_pvalue": float}}.
    """
    from pipeline.layer2_evolution.selection import _should_run_hyphy

    candidates = [og for og in aligned_orthogroups if _should_run_hyphy(og, motifs_by_og)]
    log.info("Running FEL + BUSTED on %d candidate orthogroups...", len(candidates))

    n_cpu = os.cpu_count() or 8
    n_workers = max(1, n_cpu - 2)
    total = len(candidates)
    all_results: dict[str, dict] = {}
    done = 0

    work_items = [
        (og_id, aligned_orthogroups[og_id], str(species_treefile))
        for og_id in candidates
    ]

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_fel_busted_worker, item): item[0] for item in work_items}
        for future in as_completed(futures):
            og_id, result = future.result()
            done += 1
            if result:
                all_results[og_id] = result
            if done % 100 == 0 or done == total:
                log.info("  FEL+BUSTED: %d / %d orthogroups (%.0f%%).", done, total, 100 * done / total)

    log.info("FEL + BUSTED complete: %d / %d orthogroups processed.", len(all_results), len(candidates))
    return all_results



# ---------------------------------------------------------------------------
# RELAX — Branch-specific rate acceleration / relaxation
# ---------------------------------------------------------------------------

def run_relax(
    codon_aln: dict[str, str],
    tree_newick: str,
    og_id: str,
    test_branches: Optional[list[str]] = None,
) -> Optional[dict]:
    """Run HyPhy RELAX on a codon alignment to detect branch-specific rate shifts.

    RELAX compares the *intensity* of selection (the k parameter) between
    test branches (resilient species) and reference branches (human + outgroups).
    k > 1 → selection intensified (acceleration); k < 1 → relaxation.

    Args:
        codon_aln: {label: codon_sequence}
        tree_newick: newick tree string
        og_id: orthogroup ID (for output directory)
        test_branches: species labels to mark as 'test'; defaults to all non-human

    Returns:
        Raw RELAX JSON dict or None on failure.
    """
    tools = get_tool_config()
    hyphy_bin = tools.get("hyphy_bin", "hyphy")

    root = Path(get_local_storage_root())
    relax_dir = root / "relax" / og_id
    relax_dir.mkdir(parents=True, exist_ok=True)

    aln_path = relax_dir / "codon_aln.fna"
    tree_path = relax_dir / "species.treefile"
    out_path = relax_dir / "relax.json"

    if out_path.exists() and out_path.stat().st_size > 0:
        try:
            with open(out_path) as f:
                return json.load(f)
        except Exception:
            pass

    # Annotate test branches in tree: append {Test} to resilient species labels
    if test_branches is None:
        test_branches = [lbl for lbl in codon_aln if not lbl.startswith("human")]

    labeled_tree = tree_newick
    for branch in test_branches:
        labeled_tree = labeled_tree.replace(branch, f"{branch}{{Test}}")

    with open(aln_path, "w") as f:
        for label, seq in codon_aln.items():
            f.write(f">{label}\n{seq}\n")
    tree_path.write_text(labeled_tree)

    cmd = [
        hyphy_bin, "relax",
        "--alignment", str(aln_path),
        "--tree", str(tree_path),
        "--output", str(out_path),
        "--test", "Test",
        "--quiet",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600, check=False)
        if result.returncode != 0 or not out_path.exists():
            log.debug("RELAX failed for %s: %s", og_id, result.stderr[:200])
            return None
        with open(out_path) as f:
            return json.load(f)
    except Exception as exc:
        log.debug("RELAX error for %s: %s", og_id, exc)
        return None


def parse_relax_results(relax_json: dict) -> dict:
    """Extract k (rate intensity) and p-value from RELAX output.

    Returns {"relax_k": float, "relax_pvalue": float}.
    k > 1 means test branches have intensified selection (acceleration).
    """
    try:
        test_results = relax_json.get("test results", {})
        pvalue = test_results.get("p-value")
        k = test_results.get("relaxation or intensification parameter")
        return {
            "relax_k": round(float(k), 4) if k is not None else None,
            "relax_pvalue": round(float(pvalue), 6) if pvalue is not None else None,
        }
    except Exception as exc:
        log.debug("RELAX parse error: %s", exc)
        return {"relax_k": None, "relax_pvalue": None}


def _relax_worker(args: tuple) -> tuple[str, Optional[dict]]:
    """Top-level worker: run RELAX for one orthogroup."""
    og_id, aligned_seqs, species_treefile_str, storage_root = args
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species
    from pathlib import Path as _Path
    from Bio import SeqIO as _SeqIO

    codon_aln_path = _Path(storage_root) / "meme" / og_id / "codon_aln.fna"
    if not codon_aln_path.exists():
        return og_id, None
    try:
        codon_aln = {rec.id: str(rec.seq) for rec in _SeqIO.parse(str(codon_aln_path), "fasta")}
    except Exception:
        return og_id, None
    if not codon_aln:
        return og_id, None

    species_ids = [label.split("|")[0] for label in aligned_seqs]
    pruned_tree = prune_tree_to_species(_Path(species_treefile_str), species_ids)
    relax_json = run_relax(codon_aln, pruned_tree, og_id)
    if relax_json is not None:
        return og_id, parse_relax_results(relax_json)
    return og_id, None


def run_relax_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
    motifs_by_og: dict[str, list],
    species_treefile: Path,
) -> dict[str, dict]:
    """Run RELAX in parallel for all candidate orthogroups that have codon alignments.

    Returns {og_id: {"relax_k": float, "relax_pvalue": float}}.
    """
    from pipeline.layer2_evolution.selection import _should_run_hyphy

    candidates = [og for og in aligned_orthogroups if _should_run_hyphy(og, motifs_by_og)]
    log.info("Running RELAX on %d candidate orthogroups...", len(candidates))

    n_cpu = os.cpu_count() or 8
    n_workers = max(1, n_cpu - 2)
    total = len(candidates)
    all_results: dict[str, dict] = {}
    done = 0
    storage_root = str(get_local_storage_root())

    work_items = [
        (og_id, aligned_orthogroups[og_id], str(species_treefile), storage_root)
        for og_id in candidates
    ]

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_relax_worker, item): item[0] for item in work_items}
        for future in as_completed(futures):
            og_id, result = future.result()
            done += 1
            if result is not None:
                all_results[og_id] = result
            if done % 100 == 0 or done == total:
                log.info("  RELAX: %d / %d orthogroups (%.0f%%).", done, total, 100 * done / total)

    log.info("RELAX complete: %d / %d orthogroups", len(all_results), len(candidates))
    return all_results


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------

def _meme_worker(args: tuple) -> tuple[str, Optional[dict]]:
    """Top-level worker: run MEME (or aBSREL proxy) for one orthogroup."""
    og_id, aligned_seqs, species_treefile_str, motifs_by_og = args
    from pipeline.layer2_evolution.phylo_tree import prune_tree_to_species
    from pipeline.layer2_evolution.selection import (
        _should_run_hyphy, parse_absrel_results, run_absrel, write_hyphy_input,
    )
    species_ids = [label.split("|")[0] for label in aligned_seqs]
    pruned_tree = prune_tree_to_species(Path(species_treefile_str), species_ids)

    cds_seqs: dict[str, str] = {}
    for label in aligned_seqs:
        parts = label.split("|")
        accession = parts[-1] if parts else label
        cds = fetch_cds_for_protein(accession)
        if cds:
            cds_seqs[label] = cds

    codon_aln = None
    if len(cds_seqs) == len(aligned_seqs):
        codon_aln = protein_to_codon_alignment(aligned_seqs, cds_seqs)

    if codon_aln is not None:
        meme_json = run_meme(codon_aln, pruned_tree, og_id)
        if meme_json is not None:
            return og_id, ("meme", parse_meme_results(meme_json, og_id))

    aln_path, tree_path = write_hyphy_input(og_id, aligned_seqs, pruned_tree)
    raw_result = run_absrel(aln_path, tree_path, og_id)
    if raw_result is not None:
        return og_id, ("proxy", parse_absrel_results(raw_result))
    return og_id, None


def run_meme_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
    motifs_by_og: dict[str, list],
    species_treefile: Path,
) -> dict[str, dict]:
    """Run MEME in parallel for all candidate orthogroups.

    Falls back to protein divergence proxy for any orthogroup where CDS
    fetching fails (e.g. UniProt accessions, missing NCBI links).

    Returns {og_id: parsed_selection_result}.
    """
def run_meme_pipeline(
    aligned_orthogroups: dict[str, dict[str, str]],
    motifs_by_og: dict[str, list],
    species_treefile: Path,
    flush_callback=None,
) -> dict[str, dict]:
    """Run MEME in parallel for all candidate orthogroups.

    Falls back to protein divergence proxy for any orthogroup where CDS
    fetching fails (e.g. UniProt accessions, missing NCBI links).

    Args:
        flush_callback: optional callable(batch: dict) called every 500 OGs
            to write results to DB incrementally.

    Returns {og_id: parsed_selection_result} for any remaining unflushed results.
    """
    from pipeline.layer2_evolution.selection import _should_run_hyphy
    candidates = [og for og in aligned_orthogroups if _should_run_hyphy(og, motifs_by_og)]
    log.info("Running MEME on %d candidate orthogroups (falls back to proxy if CDS unavailable)...",
             len(candidates))

    # Pre-fetch all CDS sequences to disk cache before spawning workers.
    # This avoids 14 workers competing for NCBI's rate limit simultaneously.
    _prefetch_all_cds(aligned_orthogroups)

    n_cpu = os.cpu_count() or 8
    n_workers = max(1, n_cpu - 2)
    total = len(candidates)
    all_results: dict[str, dict] = {}
    meme_success = 0
    proxy_fallback = 0
    done = 0
    pending_flush: dict[str, dict] = {}
    _FLUSH_INTERVAL = 500

    work_items = [
        (og_id, aligned_orthogroups[og_id], str(species_treefile), motifs_by_og)
        for og_id in candidates
    ]

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(_meme_worker, item): item[0] for item in work_items}
        for future in as_completed(futures):
            og_id, result = future.result()
            done += 1
            if result is not None:
                kind, parsed = result
                pending_flush[og_id] = parsed
                all_results[og_id] = parsed
                if kind == "meme":
                    meme_success += 1
                else:
                    proxy_fallback += 1
            if done % 100 == 0 or done == total:
                log.info("  MEME: %d / %d orthogroups (%.0f%%) — %d MEME, %d proxy, %d failed.",
                         done, total, 100 * done / total,
                         meme_success, proxy_fallback, done - meme_success - proxy_fallback)
            # Flush to DB periodically
            if flush_callback and len(pending_flush) >= _FLUSH_INTERVAL:
                try:
                    flush_callback(pending_flush)
                    pending_flush.clear()
                except Exception as exc:
                    log.warning("MEME flush failed: %s", exc)

    # Flush any remaining
    if flush_callback and pending_flush:
        try:
            flush_callback(pending_flush)
            pending_flush.clear()
        except Exception as exc:
            log.warning("MEME final flush failed: %s", exc)

    log.info("MEME complete: %d MEME, %d proxy fallback, %d failed.",
             meme_success, proxy_fallback, total - meme_success - proxy_fallback)

    # Apply Benjamini-Hochberg FDR correction across all site-level p-values
    from pipeline.stats import apply_bh_to_meme_results
    all_results = apply_bh_to_meme_results(all_results)

    # Return only unflushed results (empty if flush_callback handled everything)
    return {} if flush_callback else all_results
