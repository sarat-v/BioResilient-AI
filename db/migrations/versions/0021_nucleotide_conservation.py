"""Add nucleotide conservation tables.

Adds:
  - nucleotide_region: CDS / promoter / downstream sequences per gene × species
  - nucleotide_score:  per-gene per-region alignment statistics + regulatory divergence counts
  - phylo_conservation_score: phyloP / PhastCons scores per gene

Revision ID: 0021
Revises: 0020
Create Date: 2026-03-11
"""
from alembic import op
import sqlalchemy as sa

revision = "0021"
down_revision = "0020"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "nucleotide_region",
        sa.Column("id",          sa.String(), primary_key=True),
        sa.Column("gene_id",     sa.String(), sa.ForeignKey("gene.id", ondelete="CASCADE"), nullable=False),
        sa.Column("species_id",  sa.String(), sa.ForeignKey("species.id", ondelete="CASCADE"), nullable=False),
        sa.Column("region_type", sa.String(), nullable=False),   # "cds" | "promoter" | "downstream"
        sa.Column("sequence",    sa.Text()),
        sa.Column("chrom",       sa.String()),
        sa.Column("start",       sa.Integer()),
        sa.Column("end",         sa.Integer()),
        sa.UniqueConstraint("gene_id", "species_id", "region_type", name="uq_nucregion"),
    )
    op.create_index("ix_nucregion_gene_id",    "nucleotide_region", ["gene_id"])
    op.create_index("ix_nucregion_species_id", "nucleotide_region", ["species_id"])

    op.create_table(
        "nucleotide_score",
        sa.Column("id",                            sa.String(), primary_key=True),
        sa.Column("gene_id",                       sa.String(), sa.ForeignKey("gene.id", ondelete="CASCADE"), nullable=False),
        sa.Column("region_type",                   sa.String(), nullable=False),
        sa.Column("conservation_score",            sa.Float()),
        sa.Column("percent_identity",              sa.Float()),
        sa.Column("alignment_length",              sa.Integer()),
        sa.Column("gap_fraction",                  sa.Float()),
        sa.Column("regulatory_divergence_count",   sa.Integer(), server_default="0"),
        sa.Column("regulatory_convergence_count",  sa.Integer(), server_default="0"),
        sa.UniqueConstraint("gene_id", "region_type", name="uq_nucscore"),
    )
    op.create_index("ix_nucscore_gene_id", "nucleotide_score", ["gene_id"])

    op.create_table(
        "phylo_conservation_score",
        sa.Column("gene_id",                sa.String(), sa.ForeignKey("gene.id", ondelete="CASCADE"), primary_key=True),
        sa.Column("cds_phylo_score",        sa.Float()),
        sa.Column("promoter_phylo_score",   sa.Float()),
        sa.Column("downstream_phylo_score", sa.Float()),
    )


def downgrade() -> None:
    op.drop_table("phylo_conservation_score")
    op.drop_index("ix_nucscore_gene_id",    table_name="nucleotide_score")
    op.drop_table("nucleotide_score")
    op.drop_index("ix_nucregion_gene_id",    table_name="nucleotide_region")
    op.drop_index("ix_nucregion_species_id", table_name="nucleotide_region")
    op.drop_table("nucleotide_region")
