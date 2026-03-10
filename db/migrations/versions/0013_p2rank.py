"""Add P2Rank fields to drug_target.

Revision ID: 0013
Revises: 0012
Create Date: 2026-03-10
"""
from alembic import op
import sqlalchemy as sa

revision = "0013"
down_revision = "0012"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("drug_target", sa.Column("p2rank_score", sa.Float(), nullable=True))
    op.add_column("drug_target", sa.Column("p2rank_pocket_count", sa.Integer(), nullable=True))


def downgrade() -> None:
    op.drop_column("drug_target", "p2rank_pocket_count")
    op.drop_column("drug_target", "p2rank_score")
