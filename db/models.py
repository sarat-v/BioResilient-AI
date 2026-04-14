"""SQLAlchemy ORM models for the BioResilient AI database.

All six layers are defined here. Phase 1 populates:
  Species, Gene, Ortholog, DivergentMotif, EvolutionScore, CandidateScore

Phase 2 adds data to:
  DiseaseAnnotation, DrugTarget, GeneTherapyScore, SafetyFlag

Step 9b adds:
  ConvergentPositionAnnotation — per-convergent-motif structural context
"""

import uuid
from datetime import datetime
from typing import Optional

from sqlalchemy import (
    Boolean,
    Column,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.dialects.postgresql import ARRAY, JSON
from sqlalchemy.orm import DeclarativeBase, relationship


class Base(DeclarativeBase):
    pass


def _uuid() -> str:
    return str(uuid.uuid4())


# ---------------------------------------------------------------------------
# Layer 1: Species + Sequence
# ---------------------------------------------------------------------------


class Species(Base):
    __tablename__ = "species"

    id = Column(String, primary_key=True)               # e.g. "naked_mole_rat"
    taxid = Column(Integer, nullable=False, unique=True)
    scientific_name = Column(String, nullable=False)
    phenotypes = Column(ARRAY(String), nullable=False, default=list)
    genome_assembly = Column(String)
    proteome_path = Column(String)                      # local path or S3 URI
    lineage_group = Column(String)                      # Rodents, Cetaceans, Bats, etc.
    geo_search_terms = Column(ARRAY(String), default=list)
    is_control = Column(Boolean, default=False)         # True = negative control species

    orthologs = relationship("Ortholog", back_populates="species")

    def __repr__(self) -> str:
        return f"<Species {self.id} (taxid={self.taxid})>"


class Gene(Base):
    __tablename__ = "gene"

    id = Column(String, primary_key=True, default=_uuid)
    human_gene_id = Column(String, nullable=False, unique=True)   # NCBI Gene ID
    gene_symbol = Column(String, nullable=False)
    human_protein = Column(String)                                # UniProt accession
    narrative = Column(Text, nullable=True)                        # LLM-generated research summary (cached)
    go_terms = Column(ARRAY(String), default=list)                # GO biological process/function terms
    pathway_ids = Column(ARRAY(String), default=list)              # Reactome/pathway identifiers

    orthologs = relationship("Ortholog", back_populates="gene")
    evolution_score = relationship("EvolutionScore", back_populates="gene", uselist=False)
    disease_annotation = relationship("DiseaseAnnotation", back_populates="gene", uselist=False)
    drug_target = relationship("DrugTarget", back_populates="gene", uselist=False)
    gene_therapy_score = relationship("GeneTherapyScore", back_populates="gene", uselist=False)
    safety_flag = relationship("SafetyFlag", back_populates="gene", uselist=False)
    candidate_score = relationship("CandidateScore", back_populates="gene", uselist=False)
    regulatory_divergences = relationship("RegulatoryDivergence", back_populates="gene")

    def __repr__(self) -> str:
        return f"<Gene {self.gene_symbol} ({self.human_gene_id})>"


class Ortholog(Base):
    __tablename__ = "ortholog"
    __table_args__ = (
        UniqueConstraint("gene_id", "species_id", name="uq_ortholog_gene_species"),
        Index("idx_ortholog_species_protein", "species_id", "protein_id"),
    )

    id = Column(String, primary_key=True, default=_uuid)
    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), nullable=False)
    species_id = Column(String, ForeignKey("species.id", ondelete="CASCADE"), nullable=False)
    protein_id = Column(String)                         # Source protein accession
    protein_seq = Column(Text)
    sequence_identity_pct = Column(Float)               # vs. human ortholog (0–100)
    orthofinder_og = Column(String)                     # OrthoGroup ID, e.g. "OG0000042"
    is_one_to_one = Column(Boolean, default=True, nullable=False)  # False if orthogroup has paralogs

    gene = relationship("Gene", back_populates="orthologs")
    species = relationship("Species", back_populates="orthologs")
    motifs = relationship("DivergentMotif", back_populates="ortholog")

    def __repr__(self) -> str:
        return f"<Ortholog gene={self.gene_id} species={self.species_id} id={self.sequence_identity_pct:.1f}%>"


class DivergentMotif(Base):
    __tablename__ = "divergent_motif"
    __table_args__ = (
        UniqueConstraint("ortholog_id", "start_pos", name="uq_motif_ortholog_pos"),
    )

    id = Column(String, primary_key=True, default=_uuid)
    ortholog_id = Column(String, ForeignKey("ortholog.id", ondelete="CASCADE"), nullable=False)
    start_pos = Column(Integer, nullable=False)         # 0-indexed position in alignment
    end_pos = Column(Integer, nullable=False)
    animal_seq = Column(String, nullable=False)         # 10–20 aa window
    human_seq = Column(String, nullable=False)
    divergence_score = Column(Float)                    # fraction of differing positions
    esm_distance = Column(Float)                        # cosine distance of ESM-2 embeddings

    # Layer 4 peptide tractability fields (populated in Phase 2)
    half_life_min = Column(Float)
    logp = Column(Float)
    immunogenic = Column(Boolean)
    synthesisable = Column(Boolean)

    # U2: Pfam domain annotation (populated in Step 4b)
    domain_name = Column(String, nullable=True)         # Pfam/InterPro domain name overlapping this motif
    in_functional_domain = Column(Boolean, default=False, nullable=False)

    # U3: AlphaMissense functional consequence (populated in Step 4b)
    consequence_score = Column(Float, nullable=True)    # Mean AlphaMissense pathogenicity [0–1]

    # Item 2: True convergent amino acid count (populated in Step 7b)
    convergent_aa_count = Column(Integer, default=0, nullable=False)  # lineages sharing the same biochemical substitution

    # Item 5: ESM-1v variant effect score (populated in Step 4c)
    esm1v_score = Column(Float, nullable=True)          # Mean ESM-1v log-likelihood ratio [negative = destabilising]

    # A6: Variant direction inference (populated in Step 4d)
    motif_direction = Column(String, nullable=True)     # 'gain_of_function', 'loss_of_function', 'likely_pathogenic', 'neutral'

    ortholog = relationship("Ortholog", back_populates="motifs")

    def __repr__(self) -> str:
        return f"<DivergentMotif ortholog={self.ortholog_id} pos={self.start_pos}-{self.end_pos}>"


# ---------------------------------------------------------------------------
# Layer 2: Evolution
# ---------------------------------------------------------------------------


class ExpressionResult(Base):
    """Per-GEO-dataset expression evidence for traceability."""

    __tablename__ = "expression_result"
    __table_args__ = (
        UniqueConstraint("gene_id", "geo_accession", "comparison", name="uq_expression_result"),
    )

    id = Column(String, primary_key=True, default=_uuid)
    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), nullable=False, index=True)
    geo_accession = Column(String, nullable=False)
    comparison = Column(String, nullable=True)  # e.g. species_id or condition
    log2fc = Column(Float, nullable=True)
    padj = Column(Float, nullable=True)
    n_samples = Column(Integer, nullable=True)


class EvolutionScore(Base):
    __tablename__ = "evolution_score"

    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), primary_key=True)
    dnds_ratio = Column(Float)
    dnds_pvalue = Column(Float)
    selection_model = Column(String)                    # e.g. "aBSREL", "BUSTED"
    branches_under_selection = Column(ARRAY(String))    # species labels from HyPhy output
    convergence_count = Column(Integer, default=0)      # independent lineages with same motif change
    phylop_score = Column(Float)                        # mean PhyloP across divergent positions

    # U4: FEL + BUSTED supplementary selection tests (populated in Step 6b)
    fel_sites = Column(Integer, nullable=True)          # sites under pervasive selection (FEL p<0.05)
    busted_pvalue = Column(Float, nullable=True)        # gene-wide episodic selection p-value (BUSTED)

    # A2: BH-corrected q-value for global FDR control (populated post-scoring)
    meme_qvalue = Column(Float, nullable=True)          # BH-adjusted q-value across all genes

    # A5: RELAX branch-specific rate acceleration (populated in Step 6c)
    relax_k = Column(Float, nullable=True)              # RELAX k parameter (>1 = acceleration)
    relax_pvalue = Column(Float, nullable=True)         # RELAX p-value for rate shift

    # Fix 3: dedicated convergence weight column (populated in Step 7a)
    # Previously the convergence weight was incorrectly written to phylop_score when
    # phylop_score was None (conflating two distinct signals).  phylop_score now holds
    # only the UCSC PhyloP conservation score; this column holds the phylogenetically-
    # weighted motif-convergence score from convergence.py.
    convergence_weight = Column(Float, nullable=True)   # phylogenetically-weighted convergence score

    # Fix 2: permutation-based null-model p-value for convergence (populated in Step 7a)
    convergence_pval = Column(Float, nullable=True)     # empirical p-value vs. 200-iteration label shuffle

    gene = relationship("Gene", back_populates="evolution_score")

    def __repr__(self) -> str:
        return f"<EvolutionScore gene={self.gene_id} dN/dS={self.dnds_ratio}>"


# ---------------------------------------------------------------------------
# Layer 3: Disease (Phase 2)
# ---------------------------------------------------------------------------


class DiseaseAnnotation(Base):
    __tablename__ = "disease_annotation"

    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), primary_key=True)
    disease_id = Column(String)                         # EFO ID
    disease_name = Column(String)
    opentargets_score = Column(Float)
    gwas_pvalue = Column(Float)
    gnomad_pli = Column(Float)                          # loss-of-function constraint score
    mouse_ko_phenotype = Column(String)
    tissue_expression = Column(JSON)                    # {tissue: tpm_value}

    # U6: Rare protective variant mapping (populated in Step 11b)
    protective_variant_count = Column(Integer, nullable=True)   # rare variants matching animal divergence direction
    best_protective_trait = Column(String, nullable=True)       # most significant protective phenotype from PheWAS
    protective_variant_pvalue = Column(Float, nullable=True)    # best GWAS p-value for protective association

    # Item 10: Literature validation score (populated in Step 11c)
    lit_score = Column(Float, nullable=True)            # PubMed citation density in resilience/longevity literature
    lit_pmid_count = Column(Integer, nullable=True)     # Number of relevant PubMed papers

    # Phase 2 upgrades: OpenTargets extended fields (populated in Step 11)
    gnomad_loeuf = Column(Float, nullable=True)         # gnomAD v4 LOEUF (< 0.6 = LoF intolerant)
    known_drug_name = Column(String, nullable=True)     # Best-phase known drug targeting this gene
    known_drug_phase = Column(Integer, nullable=True)   # Max clinical trial phase (4 = approved)
    ot_safety_liability = Column(String, nullable=True) # Comma-separated safety event labels from OT

    gene = relationship("Gene", back_populates="disease_annotation")


# ---------------------------------------------------------------------------
# Layer 4: Druggability (Phase 2)
# ---------------------------------------------------------------------------


class DrugTarget(Base):
    __tablename__ = "drug_target"

    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), primary_key=True)
    pocket_count = Column(Integer)
    top_pocket_score = Column(Float)
    chembl_target_id = Column(String)
    existing_drugs = Column(ARRAY(String))
    cansar_score = Column(Float)
    druggability_tier = Column(String)                  # "A", "B", "C", or "undruggable"

    # Item 7: P2Rank ML pocket prediction (populated in Step 12b)
    p2rank_score = Column(Float, nullable=True)         # Best P2Rank pocket probability [0–1]
    p2rank_pocket_count = Column(Integer, nullable=True)  # Number of predicted pockets

    # Phase 2 upgrades: tractability modalities + convergent position proximity (Step 12)
    tractability_sm = Column(Boolean, nullable=True)    # Small molecule tractable (OT / DrugEBIlity)
    tractability_ab = Column(Boolean, nullable=True)    # Antibody tractable (surface/extracellular)
    tractability_protac = Column(Boolean, nullable=True)  # PROTAC-tractable (PROTACtable genome)
    convergent_pocket_proximal = Column(Boolean, nullable=True)  # Any convergent AA within 6Å of top pocket

    gene = relationship("Gene", back_populates="drug_target")


# ---------------------------------------------------------------------------
# Layer 5: Gene Therapy (Phase 2)
# ---------------------------------------------------------------------------


class GeneTherapyScore(Base):
    __tablename__ = "gene_therapy_score"

    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), primary_key=True)
    gene_size_bp = Column(Integer)
    aav_compatible = Column(Boolean)
    tissue_tropism = Column(ARRAY(String))              # matching AAV serotypes
    crispr_sites = Column(Integer)                      # number of valid guide sites
    offtarget_risk = Column(String)                     # "low", "medium", "high"

    gene = relationship("Gene", back_populates="gene_therapy_score")


# ---------------------------------------------------------------------------
# Layer 6: Safety (Phase 2)
# ---------------------------------------------------------------------------


class SafetyFlag(Base):
    __tablename__ = "safety_flag"

    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), primary_key=True)
    is_essential = Column(Boolean)                      # pLI > 0.9
    phewas_hits = Column(JSON)                          # {trait: pvalue}
    network_degree = Column(Integer)                    # STRING interaction count
    hub_risk = Column(Boolean)                          # degree > 50
    family_size = Column(Integer)                       # protein family size

    # Item 8: DepMap essentiality + GTEx expression breadth (populated in Step 14b)
    depmap_score = Column(Float, nullable=True)         # DepMap CRISPR chronos score (more negative = more essential)
    gtex_tissue_count = Column(Integer, nullable=True)  # Number of tissues with TPM > 1 (GTEx)
    gtex_max_tpm = Column(Float, nullable=True)         # Maximum median TPM across all GTEx tissues

    gene = relationship("Gene", back_populates="safety_flag")


# ---------------------------------------------------------------------------
# Composite Score
# ---------------------------------------------------------------------------


class CandidateScore(Base):
    __tablename__ = "candidate_score"

    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), primary_key=True)
    trait_id = Column(String, primary_key=True, default="", nullable=False)  # "" = default/legacy; e.g. "cancer_resistance"
    convergence_score = Column(Float, default=0.0)
    selection_score = Column(Float, default=0.0)
    disease_score = Column(Float, default=0.0)          # 0.0 in Phase 1
    druggability_score = Column(Float, default=0.0)     # 0.0 in Phase 1
    expression_score = Column(Float, default=0.0)
    safety_score = Column(Float, default=0.0)           # 0.0 in Phase 1
    composite_score = Column(Float, default=0.0)
    tier = Column(String)                               # "Tier1", "Tier2", "Tier3"
    regulatory_score = Column(Float, default=0.0)       # Track B: AlphaGenome regulatory divergence
    structural_score = Column(Float, default=0.0)       # Step 9b: structural context of convergent changes
    control_divergence_fraction = Column(Float, nullable=True)  # A3: fraction of control species also divergent
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    gene = relationship("Gene", back_populates="candidate_score")

    def __repr__(self) -> str:
        return f"<CandidateScore gene={self.gene_id} composite={self.composite_score:.3f} {self.tier}>"


# ---------------------------------------------------------------------------
# Step 9b: Convergent position structural context
# ---------------------------------------------------------------------------


class ConvergentPositionAnnotation(Base):
    """Per-convergent-motif structural annotation produced by Step 9b.

    One row per DivergentMotif that shows convergent_aa_count >= threshold.
    Each row captures the structural context of a convergently changed region
    in the canonical human protein: UniProt functional site, AlphaMissense
    pathogenicity, AlphaFold confidence, and 3D proximity to the top druggable
    pocket.  The per-gene structural_score in CandidateScore aggregates these
    rows into a single [0,1] signal for Phase 2 composite scoring.
    """

    __tablename__ = "convergent_position_annotation"
    __table_args__ = (
        UniqueConstraint("gene_id", "divergent_motif_id", name="uq_cpa_gene_motif"),
    )

    id = Column(String, primary_key=True, default=_uuid)
    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), nullable=False)
    divergent_motif_id = Column(
        String, ForeignKey("divergent_motif.id", ondelete="CASCADE"), nullable=False
    )
    uniprot_accession = Column(String, nullable=True)   # Resolved UniProt accession (e.g. "P31749")

    # Protein-space coordinates of the motif (1-indexed, UniProt canonical sequence)
    protein_start_pos = Column(Integer, nullable=True)
    protein_end_pos = Column(Integer, nullable=True)

    # Convergence strength (copied from DivergentMotif for quick access)
    convergent_lineage_count = Column(Integer, nullable=True)

    # AlphaMissense — mean score for this specific region's divergent positions
    am_score = Column(Float, nullable=True)             # Mean AM pathogenicity [0–1]
    am_class = Column(String, nullable=True)            # 'likely_pathogenic', 'ambiguous', 'likely_benign'

    # UniProt functional site annotation at this region
    uniprot_feature = Column(String, nullable=True)     # 'active_site', 'binding_site', 'modified_residue', etc.
    uniprot_feature_description = Column(String, nullable=True)

    # AlphaFold structural confidence
    plddt_mean = Column(Float, nullable=True)           # Mean pLDDT for residues in this region
    is_disordered = Column(Boolean, nullable=True)      # plddt_mean < 70

    # 3D pocket proximity (from fpocket top-ranked pocket)
    pocket_distance_angstrom = Column(Float, nullable=True)  # Å from region Cα centroid to pocket centroid
    is_pocket_adjacent = Column(Boolean, nullable=True)      # distance <= 8.0 Å

    # Derived structural context classifier
    structural_context = Column(String, nullable=True)  # 'active_site', 'binding_site', 'pocket_adjacent',
                                                        # 'surface_near_pocket', 'disordered', 'distal'

    gene = relationship("Gene")
    divergent_motif = relationship("DivergentMotif")


# ---------------------------------------------------------------------------
# Pipeline run history
# ---------------------------------------------------------------------------


class PipelineRun(Base):
    __tablename__ = "pipeline_run"

    id = Column(String, primary_key=True, default=_uuid)
    started_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    finished_at = Column(DateTime, nullable=True)
    status = Column(String, nullable=False, default="running")  # running, completed, failed, stopped
    species_ids = Column(ARRAY(String), default=list)
    step_statuses = Column(JSON, default=dict)  # {step_name: {status, label, updated_at, elapsed_s?}}
    phase = Column(String, nullable=True)
    notes = Column(Text, nullable=True)
    pid = Column(Integer, nullable=True)  # process id of the pipeline subprocess

    def __repr__(self) -> str:
        return f"<PipelineRun {self.id} status={self.status} pid={self.pid}>"


# ---------------------------------------------------------------------------
# Track B: Regulatory Divergence (AlphaGenome)
# ---------------------------------------------------------------------------


class RegulatoryDivergence(Base):
    __tablename__ = "regulatory_divergence"
    __table_args__ = (
        UniqueConstraint("gene_id", "species_id", name="uq_regulatory_divergence_gene_species"),
    )

    id = Column(String, primary_key=True, default=_uuid)
    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), nullable=False)
    species_id = Column(String, ForeignKey("species.id", ondelete="CASCADE"), nullable=False)
    promoter_divergence = Column(Float)                 # AlphaGenome predicted effect magnitude
    expression_log2fc = Column(Float)                   # log2 fold-change from DESeq2
    lineage_count = Column(Integer, default=0)          # independent lineages with regulatory signal
    regulatory_score = Column(Float, default=0.0)       # combined score per gene × species

    gene = relationship("Gene", back_populates="regulatory_divergences")
    species = relationship("Species")

    def __repr__(self) -> str:
        return f"<RegulatoryDivergence gene={self.gene_id} species={self.species_id} score={self.regulatory_score}>"


# ---------------------------------------------------------------------------
# A4: Pathway-level convergence
# ---------------------------------------------------------------------------


class PathwayConvergence(Base):
    __tablename__ = "pathway_convergence"

    pathway_id = Column(String, primary_key=True)
    pathway_name = Column(String)
    gene_count = Column(Integer)                            # background gene count in pathway
    candidate_count = Column(Integer)                       # convergent candidates in pathway
    log_pvalue = Column(Float)                              # log10 hypergeometric p-value (negative = enriched)
    evolutionary_weight = Column(Float)                     # sum of convergence_scores
    pathway_score = Column(Float)                           # combined ranking score
    gene_symbols = Column(ARRAY(String), default=list)      # gene symbols in this pathway

    def __repr__(self) -> str:
        return f"<PathwayConvergence {self.pathway_id} score={self.pathway_score}>"


# ---------------------------------------------------------------------------
# Nucleotide conservation (Steps 3c / 3d)
# ---------------------------------------------------------------------------


class NucleotideRegion(Base):
    """Raw nucleotide sequences extracted per gene × species for CDS, promoter, and downstream regions."""

    __tablename__ = "nucleotide_region"
    __table_args__ = (
        UniqueConstraint("gene_id", "species_id", "region_type", name="uq_nucregion"),
    )

    id          = Column(String, primary_key=True, default=_uuid)
    gene_id     = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), nullable=False, index=True)
    species_id  = Column(String, ForeignKey("species.id", ondelete="CASCADE"), nullable=False, index=True)
    region_type = Column(String, nullable=False)   # "cds" | "promoter" | "downstream"
    sequence    = Column(Text)
    chrom       = Column(String)
    start       = Column(Integer)
    end         = Column(Integer)

    gene    = relationship("Gene")
    species = relationship("Species")

    def __repr__(self) -> str:
        return f"<NucleotideRegion gene={self.gene_id} species={self.species_id} type={self.region_type}>"


class NucleotideScore(Base):
    """Per-gene per-region alignment statistics and regulatory divergence/convergence counts."""

    __tablename__ = "nucleotide_score"
    __table_args__ = (
        UniqueConstraint("gene_id", "region_type", name="uq_nucscore"),
    )

    id                           = Column(String, primary_key=True, default=_uuid)
    gene_id                      = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), nullable=False, index=True)
    region_type                  = Column(String, nullable=False)
    conservation_score           = Column(Float)    # (pct_identity × aln_len) / region_len; mean across resilient species
    percent_identity             = Column(Float)    # mean pct identity vs. human across resilient species
    alignment_length             = Column(Integer)  # mean alignment length
    gap_fraction                 = Column(Float)    # mean gap fraction
    regulatory_divergence_count  = Column(Integer, default=0)   # mutations in ≥3 resilient, absent human + controls
    regulatory_convergence_count = Column(Integer, default=0)   # same mutation in ≥2 distinct lineage clusters

    gene = relationship("Gene")

    def __repr__(self) -> str:
        return f"<NucleotideScore gene={self.gene_id} type={self.region_type} conservation={self.conservation_score}>"


class PhyloConservationScore(Base):
    """phyloP / PhastCons scores per gene for CDS, promoter, and downstream regions."""

    __tablename__ = "phylo_conservation_score"

    gene_id                = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), primary_key=True)
    cds_phylo_score        = Column(Float)
    promoter_phylo_score   = Column(Float)
    downstream_phylo_score = Column(Float)

    gene = relationship("Gene")

    def __repr__(self) -> str:
        return f"<PhyloConservationScore gene={self.gene_id} cds={self.cds_phylo_score}>"


# ---------------------------------------------------------------------------
# Auth: Platform users
# ---------------------------------------------------------------------------


class User(Base):
    """BioResilient platform user accounts.

    Roles:
      admin    — full access, can manage other users
      operator — can view data and trigger pipeline runs
      viewer   — read-only access to candidates and research pages
    """

    __tablename__ = "platform_user"

    id              = Column(String, primary_key=True, default=_uuid)
    email           = Column(String, nullable=False, unique=True, index=True)
    name            = Column(String, nullable=False)
    hashed_password = Column(String, nullable=False)
    role            = Column(String, nullable=False, default="viewer")   # admin | operator | viewer
    is_active       = Column(Boolean, nullable=False, default=True)
    created_at      = Column(DateTime, default=datetime.utcnow, nullable=False)

    def __repr__(self) -> str:
        return f"<User {self.email} role={self.role} active={self.is_active}>"
