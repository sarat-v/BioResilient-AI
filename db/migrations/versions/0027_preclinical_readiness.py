"""Add preclinical_readiness table and crispr_efficiency column.

Revision ID: 0027
Revises: 0026
Create Date: 2026-04-14
"""

from alembic import op
import sqlalchemy as sa

revision = "0027"
down_revision = "0026"
branch_labels = None
depends_on = None


def upgrade():
    # New preclinical_readiness table
    op.create_table(
        "preclinical_readiness",
        sa.Column("gene_id", sa.String(), sa.ForeignKey("gene.id", ondelete="CASCADE"),
                  nullable=False, primary_key=True),
        sa.Column("synthesis_score",      sa.Float(), nullable=True),
        sa.Column("deliverability_score", sa.Float(), nullable=True),
        sa.Column("model_score",          sa.Float(), nullable=True),
        sa.Column("assay_score",          sa.Float(), nullable=True),
        sa.Column("overall_readiness_score", sa.Float(), nullable=True),
        sa.Column("readiness_tier",       sa.String(), nullable=True),
        sa.Column("notes",                sa.String(), nullable=True),
        sa.Column("updated_at",           sa.DateTime(), nullable=True),
    )

    # CRISPR on-target efficiency score on existing gene_therapy_score table
    op.add_column(
        "gene_therapy_score",
        sa.Column("crispr_efficiency", sa.Float(), nullable=True),
    )


def downgrade():
    op.drop_column("gene_therapy_score", "crispr_efficiency")
    op.drop_table("preclinical_readiness")
