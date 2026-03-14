"""Step 3c (part 1) — Nucleotide region extraction.

For each gene × species pair present in the Ortholog table, this module:
  1. Downloads the genome annotation (GFF3) and genome sequence (FASTA) from NCBI
     using the assembly accession recorded in the Species table.
  2. Locates the gene's coordinates in the GFF3 by matching against the ortholog
     protein_id (mapped to a gene feature via NCBI Gene ID lookup).
  3. Extracts three sequence windows:
       CDS        — the annotated coding sequence
       promoter   — 5,000 bp upstream of the transcription start site
       downstream — 2,000 bp downstream of the gene end
  4. Upserts NucleotideRegion rows into the database.

Downloaded genome files are cached to {storage_root}/genomes/{species_id}/ and
synced to S3 so they survive EC2 Spot interruptions.

Entry point: run_nucleotide_scan(gene_ids=None) -> int
"""

import gzip
import io
import json
import logging
import re
import time
from pathlib import Path
from typing import Optional
from urllib.parse import urlencode

import requests
from Bio import SeqIO
from Bio.Seq import Seq

from db.models import Gene, NucleotideRegion, Ortholog, Species
from db.session import get_session
from pipeline.config import (
    get_ncbi_api_key,
    get_ncbi_email,
    get_storage_root,
    sync_to_s3,
    sync_from_s3,
)

log = logging.getLogger(__name__)

_PROMOTER_UPSTREAM_BP   = 5000
_DOWNSTREAM_FLANKING_BP = 2000
_NCBI_EFETCH      = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_NCBI_ESEARCH     = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_NCBI_DATASETS    = "https://api.ncbi.nlm.nih.gov/datasets/v2"
_RATE_SLEEP       = 0.11   # 10 req/s with API key


# ─────────────────────────────────────────────────────────────────────────────
# Storage helpers
# ─────────────────────────────────────────────────────────────────────────────

def _genomes_dir(species_id: str) -> Path:
    root = get_storage_root()
    if root.startswith("s3://"):
        base = Path("/tmp/bioresilient/genomes")
    else:
        base = Path(root) / "genomes"
    d = base / species_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def _ncbi_headers() -> dict:
    key = get_ncbi_api_key()
    return {"api-key": key} if key else {}


def _ncbi_params(**kwargs) -> dict:
    params = {"api_key": get_ncbi_api_key(), "email": get_ncbi_email()}
    params.update(kwargs)
    return {k: v for k, v in params.items() if v}


# ─────────────────────────────────────────────────────────────────────────────
# GFF3 download
# ─────────────────────────────────────────────────────────────────────────────

def _download_gff3(species_id: str, assembly: str) -> Optional[Path]:
    """Download and cache the GFF3 annotation for a genome assembly.

    Cache key includes assembly accession so stale files from old assemblies
    are not reused when the assembly is updated.
    """
    gdir = _genomes_dir(species_id)
    gff_path = gdir / f"{species_id}.gff3.gz"
    s3_key   = f"genomes/{species_id}/{assembly}/{species_id}.gff3.gz"

    # Try restoring from S3 first (assembly-versioned key)
    if not gff_path.exists():
        sync_from_s3(s3_key, gff_path)

    if gff_path.exists() and gff_path.stat().st_size > 1000:
        log.info("  GFF3 cache hit: %s", gff_path)
        return gff_path

    log.info("  Downloading GFF3 for %s (assembly=%s)...", species_id, assembly)
    url = f"{_NCBI_DATASETS}/genome/accession/{assembly}/download"
    params = {"include_annotation_type": "GENOME_GFF", "api-key": get_ncbi_api_key()}
    try:
        r = requests.get(url, params={k: v for k, v in params.items() if v},
                         headers={"Accept": "application/zip"}, stream=True, timeout=300)
        if r.status_code != 200:
            log.warning("  GFF3 download HTTP %s for %s", r.status_code, species_id)
            return None
        import zipfile
        zdata = io.BytesIO(r.content)
        with zipfile.ZipFile(zdata) as zf:
            gff_entries = [n for n in zf.namelist() if n.endswith(".gff") or n.endswith(".gff.gz") or n.endswith(".gff3") or n.endswith(".gff3.gz")]
            if not gff_entries:
                log.warning("  No GFF3 file found in Datasets package for %s", species_id)
                return None
            with zf.open(gff_entries[0]) as src:
                raw = src.read()
            # Normalise to gzipped
            if not gff_entries[0].endswith(".gz"):
                import gzip as _gz
                raw = _gz.compress(raw)
            gff_path.write_bytes(raw)
        sync_to_s3(gff_path, s3_key)
        log.info("  ✓ GFF3 downloaded for %s (%d bytes)", species_id, gff_path.stat().st_size)
        return gff_path
    except Exception as exc:
        log.warning("  GFF3 download failed for %s: %s", species_id, exc)
        return None


def _download_genome_fasta(species_id: str, assembly: str) -> Optional[Path]:
    """Download and cache the genomic FASTA for sequence extraction.

    Cache key includes assembly accession to avoid stale files from old assemblies.
    """
    gdir = _genomes_dir(species_id)
    fa_path = gdir / f"{species_id}.genome.fa.gz"
    s3_key  = f"genomes/{species_id}/{assembly}/{species_id}.genome.fa.gz"

    if not fa_path.exists():
        sync_from_s3(s3_key, fa_path)
    if fa_path.exists() and fa_path.stat().st_size > 100_000:
        return fa_path

    log.info("  Downloading genomic FASTA for %s...", species_id)
    url = f"{_NCBI_DATASETS}/genome/accession/{assembly}/download"
    params = {"include_annotation_type": "GENOME_FASTA", "api-key": get_ncbi_api_key()}
    try:
        r = requests.get(url, params={k: v for k, v in params.items() if v},
                         headers={"Accept": "application/zip"}, stream=True, timeout=600)
        if r.status_code != 200:
            return None
        import zipfile
        zdata = io.BytesIO(r.content)
        with zipfile.ZipFile(zdata) as zf:
            fa_entries = [n for n in zf.namelist() if "genomic" in n and (n.endswith(".fna") or n.endswith(".fna.gz") or n.endswith(".fa") or n.endswith(".fa.gz"))]
            if not fa_entries:
                return None
            with zf.open(fa_entries[0]) as src:
                raw = src.read()
            if not fa_entries[0].endswith(".gz"):
                import gzip as _gz
                raw = _gz.compress(raw)
            fa_path.write_bytes(raw)
        sync_to_s3(fa_path, s3_key)
        log.info("  ✓ Genome FASTA for %s (%d bytes)", species_id, fa_path.stat().st_size)
        return fa_path
    except Exception as exc:
        log.warning("  Genomic FASTA download failed for %s: %s", species_id, exc)
        return None


# ─────────────────────────────────────────────────────────────────────────────
# GFF3 parsing
# ─────────────────────────────────────────────────────────────────────────────

def _parse_gff3_genes(gff_path: Path) -> dict[str, dict]:
    """Parse GFF3 and return a dict mapping identifiers → {chrom, start, end, strand, tss}.

    Indexes the following feature types:
    - gene / mRNA / transcript  → keyed by Name, gene, Dbxref, ID attributes
    - CDS                       → keyed by protein_id attribute (matches Ortholog.protein_id)

    For CDS features the coordinates are expanded to the parent gene span so
    promoter / downstream windows are computed correctly. We store the widest
    CDS span seen per protein_id.

    Returns 1-indexed inclusive coordinates (GFF3 standard).
    """
    opener = gzip.open if gff_path.suffix == ".gz" else open
    genes: dict[str, dict] = {}
    # Track CDS spans separately: protein_id → {chrom, min_start, max_end, strand}
    cds_spans: dict[str, dict] = {}

    try:
        with opener(gff_path, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                seqid, _, feat_type, start, end, _, strand, _, attrs = parts
                try:
                    s, e = int(start), int(end)
                except ValueError:
                    continue

                # Parse attributes into dict
                attr_dict: dict[str, str] = {}
                for part in attrs.split(";"):
                    part = part.strip()
                    if "=" in part:
                        k, _, v = part.partition("=")
                        attr_dict[k.strip()] = v.strip()

                coord = {
                    "chrom":  seqid,
                    "start":  s,
                    "end":    e,
                    "strand": strand,
                    "tss":    s if strand == "+" else e,
                }

                if feat_type in ("gene", "mRNA", "transcript"):
                    # Index by every useful attribute
                    for key in ("gene", "Name", "gene_name", "ID"):
                        val = attr_dict.get(key, "")
                        if val:
                            genes[val.upper()] = coord
                    # Dbxref can be comma-separated: "GeneID:12345,Genbank:NP_xxx"
                    for item in attr_dict.get("Dbxref", "").split(","):
                        item = item.strip()
                        if item:
                            genes[item.upper()] = coord

                elif feat_type == "CDS":
                    # protein_id attribute holds the RefSeq protein accession (e.g. XP_004853481.1)
                    prot_id = attr_dict.get("protein_id", "").strip()
                    if prot_id:
                        key = prot_id.upper()
                        if key not in cds_spans:
                            cds_spans[key] = {"chrom": seqid, "min_s": s, "max_e": e, "strand": strand}
                        else:
                            cds_spans[key]["min_s"] = min(cds_spans[key]["min_s"], s)
                            cds_spans[key]["max_e"] = max(cds_spans[key]["max_e"], e)
                        # Also index version-stripped accession (NP_001234.1 → NP_001234)
                        base = prot_id.split(".")[0].upper()
                        if base != key:
                            cds_spans.setdefault(base, cds_spans[key])

    except Exception as exc:
        log.warning("  GFF3 parse error for %s: %s", gff_path, exc)

    # Merge CDS spans into genes index
    for prot_key, span in cds_spans.items():
        gs = span["min_s"]
        ge = span["max_e"]
        st = span["strand"]
        genes[prot_key] = {
            "chrom":  span["chrom"],
            "start":  gs,
            "end":    ge,
            "strand": st,
            "tss":    gs if st == "+" else ge,
        }

    return genes


def _find_gene_coords(gene_symbol: str, protein_id: str, gff_index: dict[str, dict]) -> Optional[dict]:
    """Look up gene coordinates from the GFF3 index.

    The Ortholog table stores protein_id in the format "{species_id}|{accession}"
    (e.g. "naked_mole_rat|XP_004853481.1") because reheader_fasta() prefixes
    headers. We strip the species prefix before lookup.

    For human, gene_symbol is stored in UniProt format "GENE_HUMAN" (e.g. "TP53_HUMAN").
    The human GFF3 uses bare gene symbols (e.g. "TP53") in the gene= attribute.

    Tries multiple matching strategies in priority order:
    1. Raw accession after stripping species prefix — most reliable for non-human
    2. Version-stripped accession
    3. Bare gene symbol — strips UniProt species suffix (TP53_HUMAN → TP53)
    4. Full gene symbol as-is
    5. GeneID: prefixed Dbxref lookup
    """
    # Strip species prefix if present (format: "species_id|accession")
    raw_accession = protein_id.split("|")[-1] if protein_id and "|" in protein_id else protein_id

    # Strip UniProt species suffix from gene symbol (TP53_HUMAN → TP53)
    bare_symbol = gene_symbol
    if gene_symbol and "_" in gene_symbol:
        parts = gene_symbol.rsplit("_", 1)
        # Only strip if suffix looks like a species tag (all caps, 2-6 chars)
        if len(parts[1]) <= 6 and parts[1].isupper():
            bare_symbol = parts[0]

    candidates = [
        raw_accession.upper() if raw_accession else "",
        raw_accession.split(".")[0].upper() if raw_accession else "",
        bare_symbol.upper() if bare_symbol else "",
        gene_symbol.upper() if gene_symbol else "",
        f"GENEID:{raw_accession}".upper() if raw_accession else "",
        f"GENE:{bare_symbol}".upper() if bare_symbol else "",
    ]
    for key in candidates:
        if key and key in gff_index:
            return gff_index[key]
    return None


# ─────────────────────────────────────────────────────────────────────────────
# Sequence extraction
# ─────────────────────────────────────────────────────────────────────────────

def _load_genome_index(fa_path: Path) -> dict[str, Seq]:
    """Load a gzipped genomic FASTA into a dict: chrom → Seq.

    Large genomes are memory-expensive; we stream-parse and keep only the
    sequences needed. For the first call we build and pickle an index of
    chrom → (offset, length) for selective extraction, but for simplicity
    we load the full genome here (manageable at ≤3 GB for most assemblies).
    """
    log.info("    Loading genome index from %s...", fa_path.name)
    opener = gzip.open if fa_path.suffix == ".gz" else open
    seqs: dict[str, Seq] = {}
    with opener(fa_path, "rt") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            seqs[rec.id] = rec.seq
    log.info("    Loaded %d chromosomes/scaffolds", len(seqs))
    return seqs


def _extract_window(genome: dict[str, Seq], chrom: str, start: int, end: int) -> Optional[str]:
    """Extract a sequence window (1-indexed inclusive) from the genome dict."""
    seq = genome.get(chrom)
    if seq is None:
        return None
    # Convert to 0-indexed Python slice
    s = max(0, start - 1)
    e = min(len(seq), end)
    if s >= e:
        return None
    return str(seq[s:e])


# ─────────────────────────────────────────────────────────────────────────────
# Region extraction per gene × species
# ─────────────────────────────────────────────────────────────────────────────

def extract_regions_for_gene(
    gene: Gene,
    species: Species,
    orthologs_by_species: dict[str, Ortholog],
    gff_index: dict[str, dict],
    genome: dict[str, Seq],
) -> list[dict]:
    """Extract CDS / promoter / downstream nucleotide windows for one gene × species.

    Returns a list of dicts with keys matching NucleotideRegion fields.
    Empty list if the gene cannot be located in the GFF3.
    """
    orth = orthologs_by_species.get(species.id)
    if orth is None:
        return []

    coords = _find_gene_coords(gene.gene_symbol, orth.protein_id or "", gff_index)
    if coords is None:
        log.debug("    No GFF3 coords found for %s in %s", gene.gene_symbol, species.id)
        return []

    chrom   = coords["chrom"]
    gstart  = coords["start"]
    gend    = coords["end"]
    strand  = coords["strand"]
    tss     = coords["tss"]

    regions = []
    for region_type, (win_start, win_end) in [
        ("cds",        (gstart, gend)),
        ("promoter",   (tss - _PROMOTER_UPSTREAM_BP,  tss - 1) if strand == "+" else (tss + 1, tss + _PROMOTER_UPSTREAM_BP)),
        ("downstream", (gend + 1, gend + _DOWNSTREAM_FLANKING_BP) if strand == "+" else (gstart - _DOWNSTREAM_FLANKING_BP, gstart - 1)),
    ]:
        seq = _extract_window(genome, chrom, max(1, win_start), win_end)
        if seq is None:
            continue
        regions.append({
            "gene_id":     gene.id,
            "species_id":  species.id,
            "region_type": region_type,
            "sequence":    seq,
            "chrom":       chrom,
            "start":       win_start,
            "end":         win_end,
        })
    return regions


# ─────────────────────────────────────────────────────────────────────────────
# DB upsert
# ─────────────────────────────────────────────────────────────────────────────

def _upsert_regions(session, region_dicts: list[dict]) -> int:
    """Upsert NucleotideRegion rows. Returns count of rows written."""
    import uuid as _uuid_mod
    written = 0
    for rd in region_dicts:
        existing = (
            session.query(NucleotideRegion)
            .filter_by(gene_id=rd["gene_id"], species_id=rd["species_id"], region_type=rd["region_type"])
            .first()
        )
        if existing is None:
            existing = NucleotideRegion(id=str(_uuid_mod.uuid4()), **rd)
            session.add(existing)
        else:
            for k, v in rd.items():
                setattr(existing, k, v)
        written += 1
    return written


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def run_nucleotide_scan(gene_ids: list[str] | None = None) -> int:
    """Extract nucleotide regions for all (or specified) genes across all species.

    For each species:
      - Downloads GFF3 annotation (cached)
      - Downloads genomic FASTA (cached)
      - For each gene with an ortholog in this species, extracts CDS/promoter/downstream
      - Upserts NucleotideRegion rows

    Returns total number of regions written to the database.
    """
    total_written = 0

    with get_session() as session:
        species_list = session.query(Species).all()
        q = session.query(Gene)
        if gene_ids:
            q = q.filter(Gene.id.in_(gene_ids))
        genes = q.all()
        log.info("Nucleotide scan: %d genes × %d species", len(genes), len(species_list))

        # Pre-index orthologs: {gene_id: {species_id: Ortholog}}
        all_orthologs = session.query(Ortholog).all()
        og_index: dict[str, dict[str, Ortholog]] = {}
        for orth in all_orthologs:
            og_index.setdefault(orth.gene_id, {})[orth.species_id] = orth

        for sp in species_list:
            if not sp.genome_assembly:
                log.info("  Skipping %s — no genome assembly recorded", sp.id)
                continue

            log.info("  Processing species: %s (%s)", sp.id, sp.genome_assembly)

            gff_path = _download_gff3(sp.id, sp.genome_assembly)
            if gff_path is None:
                log.warning("  No GFF3 for %s — skipping", sp.id)
                continue

            fa_path = _download_genome_fasta(sp.id, sp.genome_assembly)
            if fa_path is None:
                log.warning("  No genomic FASTA for %s — skipping sequence extraction", sp.id)
                continue

            gff_index = _parse_gff3_genes(gff_path)
            log.info("  GFF3 index: %d entries for %s", len(gff_index), sp.id)

            genome = _load_genome_index(fa_path)
            if not genome:
                log.warning("  Empty genome for %s — skipping", sp.id)
                continue

            species_written = 0
            genes_matched = 0
            for gene in genes:
                species_orthologs = og_index.get(gene.id, {})
                region_dicts = extract_regions_for_gene(
                    gene, sp, species_orthologs, gff_index, genome
                )
                if region_dicts:
                    species_written += _upsert_regions(session, region_dicts)
                    genes_matched += 1

            session.commit()
            log.info("  %s: %d regions written (%d/%d genes matched in GFF3)",
                     sp.id, species_written, genes_matched, len(genes))
            total_written += species_written

            # Free genome memory before next species
            del genome

            time.sleep(_RATE_SLEEP)

    log.info("Nucleotide scan complete: %d total regions written", total_written)
    return total_written
