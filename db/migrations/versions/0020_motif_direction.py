"""Add motif_direction (GoF/LoF classification) to divergent_motif.

Revision ID: 0020
Revises: 0019
Create Date: 2026-03-01
"""
from alembic import op
import sqlalchemy as sa

revision = "0020"
down_revision = "0019"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("divergent_motif", sa.Column("motif_direction", sa.String(), nullable=True))
    op.create_index("ix_divergent_motif_direction", "divergent_motif", ["motif_direction"])


def downgrade() -> None:
    op.drop_index("ix_divergent_motif_direction", table_name="divergent_motif")
    op.drop_column("divergent_motif", "motif_direction")
