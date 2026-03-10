"""Create pathway_convergence table for pathway-level enrichment scoring.

Revision ID: 0019
Revises: 0018
Create Date: 2026-03-01
"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

revision = "0019"
down_revision = "0018"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "pathway_convergence",
        sa.Column("pathway_id", sa.String(), nullable=False),
        sa.Column("pathway_name", sa.String(), nullable=True),
        sa.Column("gene_count", sa.Integer(), nullable=True),
        sa.Column("candidate_count", sa.Integer(), nullable=True),
        sa.Column("log_pvalue", sa.Float(), nullable=True),
        sa.Column("evolutionary_weight", sa.Float(), nullable=True),
        sa.Column("pathway_score", sa.Float(), nullable=True),
        sa.Column("gene_symbols", postgresql.ARRAY(sa.String()), nullable=True),
        sa.PrimaryKeyConstraint("pathway_id"),
    )
    op.create_index("ix_pathway_convergence_score", "pathway_convergence", ["pathway_score"])


def downgrade() -> None:
    op.drop_index("ix_pathway_convergence_score", table_name="pathway_convergence")
    op.drop_table("pathway_convergence")
