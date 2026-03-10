"""SQLAlchemy ORM models for the BioResilient AI database.

All six layers are defined here. Phase 1 populates:
  Species, Gene, Ortholog, DivergentMotif, EvolutionScore, CandidateScore

Phase 2 adds data to:
  DiseaseAnnotation, DrugTarget, GeneTherapyScore, SafetyFlag
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
    )

    id = Column(String, primary_key=True, default=_uuid)
    gene_id = Column(String, ForeignKey("gene.id", ondelete="CASCADE"), nullable=False)
    species_id = Column(String, ForeignKey("species.id", ondelete="CASCADE"), nullable=False)
    protein_id = Column(String)                         # Source protein accession
    protein_seq = Column(Text)
    sequence_identity_pct = Column(Float)               # vs. human ortholog (0–100)
    orthofinder_og = Column(String)                     # OrthoGroup ID, e.g. "OG0000042"

    gene = relationship("Gene", back_populates="orthologs")
    species = relationship("Species", back_populates="orthologs")
    motifs = relationship("DivergentMotif", back_populates="ortholog")

    def __repr__(self) -> str:
        return f"<Ortholog gene={self.gene_id} species={self.species_id} id={self.sequence_identity_pct:.1f}%>"


class DivergentMotif(Base):
    __tablename__ = "divergent_motif"

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

    ortholog = relationship("Ortholog", back_populates="motifs")

    def __repr__(self) -> str:
        return f"<DivergentMotif ortholog={self.ortholog_id} pos={self.start_pos}-{self.end_pos}>"


# ---------------------------------------------------------------------------
# Layer 2: Evolution
# ---------------------------------------------------------------------------


class ExpressionResult(Base):
    """Per-GEO-dataset expression evidence for traceability."""

    __tablename__ = "expression_result"

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
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    gene = relationship("Gene", back_populates="candidate_score")

    def __repr__(self) -> str:
        return f"<CandidateScore gene={self.gene_id} composite={self.composite_score:.3f} {self.tier}>"


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
