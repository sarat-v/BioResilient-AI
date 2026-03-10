"""Add DepMap + GTEx safety fields to safety_flag.

Revision ID: 0014
Revises: 0013
Create Date: 2026-03-10
"""
from alembic import op
import sqlalchemy as sa

revision = "0014"
down_revision = "0013"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("safety_flag", sa.Column("depmap_score", sa.Float(), nullable=True))
    op.add_column("safety_flag", sa.Column("gtex_tissue_count", sa.Integer(), nullable=True))
    op.add_column("safety_flag", sa.Column("gtex_max_tpm", sa.Float(), nullable=True))


def downgrade() -> None:
    op.drop_column("safety_flag", "gtex_max_tpm")
    op.drop_column("safety_flag", "gtex_tissue_count")
    op.drop_column("safety_flag", "depmap_score")
