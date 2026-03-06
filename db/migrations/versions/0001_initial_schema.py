"""Initial schema — all six layers

Revision ID: 0001
Revises:
Create Date: 2026-03-01
"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects import postgresql

revision: str = "0001"
down_revision: Union[str, None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        "species",
        sa.Column("id", sa.String(), nullable=False),
        sa.Column("taxid", sa.Integer(), nullable=False),
        sa.Column("scientific_name", sa.String(), nullable=False),
        sa.Column("phenotypes", postgresql.ARRAY(sa.String()), nullable=False),
        sa.Column("genome_assembly", sa.String(), nullable=True),
        sa.Column("proteome_path", sa.String(), nullable=True),
        sa.Column("lineage_group", sa.String(), nullable=True),
        sa.Column("geo_search_terms", postgresql.ARRAY(sa.String()), nullable=True),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("taxid"),
    )

    op.create_table(
        "gene",
        sa.Column("id", sa.String(), nullable=False),
        sa.Column("human_gene_id", sa.String(), nullable=False),
        sa.Column("gene_symbol", sa.String(), nullable=False),
        sa.Column("human_protein", sa.String(), nullable=True),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("human_gene_id"),
    )

    op.create_table(
        "ortholog",
        sa.Column("id", sa.String(), nullable=False),
        sa.Column("gene_id", sa.String(), nullable=False),
        sa.Column("species_id", sa.String(), nullable=False),
        sa.Column("protein_id", sa.String(), nullable=True),
        sa.Column("protein_seq", sa.Text(), nullable=True),
        sa.Column("sequence_identity_pct", sa.Float(), nullable=True),
        sa.Column("orthofinder_og", sa.String(), nullable=True),
        sa.ForeignKeyConstraint(["gene_id"], ["gene.id"], ondelete="CASCADE"),
        sa.ForeignKeyConstraint(["species_id"], ["species.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("gene_id", "species_id", name="uq_ortholog_gene_species"),
    )

    op.create_index("ix_ortholog_gene_id", "ortholog", ["gene_id"])
    op.create_index("ix_ortholog_orthofinder_og", "ortholog", ["orthofinder_og"])

    op.create_table(
        "divergent_motif",
        sa.Column("id", sa.String(), nullable=False),
        sa.Column("ortholog_id", sa.String(), nullable=False),
        sa.Column("start_pos", sa.Integer(), nullable=False),
        sa.Column("end_pos", sa.Integer(), nullable=False),
        sa.Column("animal_seq", sa.String(), nullable=False),
        sa.Column("human_seq", sa.String(), nullable=False),
        sa.Column("divergence_score", sa.Float(), nullable=True),
        sa.Column("esm_distance", sa.Float(), nullable=True),
        sa.Column("half_life_min", sa.Float(), nullable=True),
        sa.Column("logp", sa.Float(), nullable=True),
        sa.Column("immunogenic", sa.Boolean(), nullable=True),
        sa.Column("synthesisable", sa.Boolean(), nullable=True),
        sa.ForeignKeyConstraint(["ortholog_id"], ["ortholog.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("id"),
    )

    op.create_index("ix_divergent_motif_ortholog_id", "divergent_motif", ["ortholog_id"])

    op.create_table(
        "evolution_score",
        sa.Column("gene_id", sa.String(), nullable=False),
        sa.Column("dnds_ratio", sa.Float(), nullable=True),
        sa.Column("dnds_pvalue", sa.Float(), nullable=True),
        sa.Column("selection_model", sa.String(), nullable=True),
        sa.Column("branches_under_selection", postgresql.ARRAY(sa.String()), nullable=True),
        sa.Column("convergence_count", sa.Integer(), nullable=True),
        sa.Column("phylop_score", sa.Float(), nullable=True),
        sa.ForeignKeyConstraint(["gene_id"], ["gene.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("gene_id"),
    )

    op.create_table(
        "disease_annotation",
        sa.Column("gene_id", sa.String(), nullable=False),
        sa.Column("disease_id", sa.String(), nullable=True),
        sa.Column("disease_name", sa.String(), nullable=True),
        sa.Column("opentargets_score", sa.Float(), nullable=True),
        sa.Column("gwas_pvalue", sa.Float(), nullable=True),
        sa.Column("gnomad_pli", sa.Float(), nullable=True),
        sa.Column("mouse_ko_phenotype", sa.String(), nullable=True),
        sa.Column("tissue_expression", postgresql.JSON(), nullable=True),
        sa.ForeignKeyConstraint(["gene_id"], ["gene.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("gene_id"),
    )

    op.create_table(
        "drug_target",
        sa.Column("gene_id", sa.String(), nullable=False),
        sa.Column("pocket_count", sa.Integer(), nullable=True),
        sa.Column("top_pocket_score", sa.Float(), nullable=True),
        sa.Column("chembl_target_id", sa.String(), nullable=True),
        sa.Column("existing_drugs", postgresql.ARRAY(sa.String()), nullable=True),
        sa.Column("cansar_score", sa.Float(), nullable=True),
        sa.Column("druggability_tier", sa.String(), nullable=True),
        sa.ForeignKeyConstraint(["gene_id"], ["gene.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("gene_id"),
    )

    op.create_table(
        "gene_therapy_score",
        sa.Column("gene_id", sa.String(), nullable=False),
        sa.Column("gene_size_bp", sa.Integer(), nullable=True),
        sa.Column("aav_compatible", sa.Boolean(), nullable=True),
        sa.Column("tissue_tropism", postgresql.ARRAY(sa.String()), nullable=True),
        sa.Column("crispr_sites", sa.Integer(), nullable=True),
        sa.Column("offtarget_risk", sa.String(), nullable=True),
        sa.ForeignKeyConstraint(["gene_id"], ["gene.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("gene_id"),
    )

    op.create_table(
        "safety_flag",
        sa.Column("gene_id", sa.String(), nullable=False),
        sa.Column("is_essential", sa.Boolean(), nullable=True),
        sa.Column("phewas_hits", postgresql.JSON(), nullable=True),
        sa.Column("network_degree", sa.Integer(), nullable=True),
        sa.Column("hub_risk", sa.Boolean(), nullable=True),
        sa.Column("family_size", sa.Integer(), nullable=True),
        sa.ForeignKeyConstraint(["gene_id"], ["gene.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("gene_id"),
    )

    op.create_table(
        "candidate_score",
        sa.Column("gene_id", sa.String(), nullable=False),
        sa.Column("convergence_score", sa.Float(), nullable=True),
        sa.Column("selection_score", sa.Float(), nullable=True),
        sa.Column("disease_score", sa.Float(), nullable=True),
        sa.Column("druggability_score", sa.Float(), nullable=True),
        sa.Column("expression_score", sa.Float(), nullable=True),
        sa.Column("safety_score", sa.Float(), nullable=True),
        sa.Column("composite_score", sa.Float(), nullable=True),
        sa.Column("tier", sa.String(), nullable=True),
        sa.Column("updated_at", sa.DateTime(), nullable=True),
        sa.ForeignKeyConstraint(["gene_id"], ["gene.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("gene_id"),
    )

    op.create_index("ix_candidate_score_composite", "candidate_score", ["composite_score"])
    op.create_index("ix_candidate_score_tier", "candidate_score", ["tier"])


def downgrade() -> None:
    op.drop_table("candidate_score")
    op.drop_table("safety_flag")
    op.drop_table("gene_therapy_score")
    op.drop_table("drug_target")
    op.drop_table("disease_annotation")
    op.drop_table("evolution_score")
    op.drop_table("divergent_motif")
    op.drop_table("ortholog")
    op.drop_table("gene")
    op.drop_table("species")
