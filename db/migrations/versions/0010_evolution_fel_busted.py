"""Add FEL and BUSTED selection fields to evolution_score.

Revision ID: 0010
Revises: 0009
Create Date: 2026-03-10
"""
from alembic import op
import sqlalchemy as sa

revision = "0010"
down_revision = "0009"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("evolution_score", sa.Column("fel_sites", sa.Integer(), nullable=True))
    op.add_column("evolution_score", sa.Column("busted_pvalue", sa.Float(), nullable=True))


def downgrade() -> None:
    op.drop_column("evolution_score", "busted_pvalue")
    op.drop_column("evolution_score", "fel_sites")
