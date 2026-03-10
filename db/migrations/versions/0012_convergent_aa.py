"""Add convergent_aa_count to divergent_motif.

Revision ID: 0012
Revises: 0011
Create Date: 2026-03-10
"""
from alembic import op
import sqlalchemy as sa

revision = "0012"
down_revision = "0011"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column(
        "divergent_motif",
        sa.Column("convergent_aa_count", sa.Integer(), nullable=False, server_default="0"),
    )
    op.create_index(
        "ix_divergent_motif_conv_aa",
        "divergent_motif",
        ["convergent_aa_count"],
    )


def downgrade() -> None:
    op.drop_index("ix_divergent_motif_conv_aa", table_name="divergent_motif")
    op.drop_column("divergent_motif", "convergent_aa_count")
