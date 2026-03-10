"""Add ESM-1v score to divergent_motif.

Revision ID: 0016
Revises: 0015
Create Date: 2026-03-10
"""
from alembic import op
import sqlalchemy as sa

revision = "0016"
down_revision = "0015"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("divergent_motif", sa.Column("esm1v_score", sa.Float(), nullable=True))


def downgrade() -> None:
    op.drop_column("divergent_motif", "esm1v_score")
