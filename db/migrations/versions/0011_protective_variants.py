"""Add rare protective variant fields to disease_annotation.

Revision ID: 0011
Revises: 0010
Create Date: 2026-03-10
"""
from alembic import op
import sqlalchemy as sa

revision = "0011"
down_revision = "0010"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("disease_annotation", sa.Column("protective_variant_count", sa.Integer(), nullable=True))
    op.add_column("disease_annotation", sa.Column("best_protective_trait", sa.String(), nullable=True))
    op.add_column("disease_annotation", sa.Column("protective_variant_pvalue", sa.Float(), nullable=True))


def downgrade() -> None:
    op.drop_column("disease_annotation", "protective_variant_pvalue")
    op.drop_column("disease_annotation", "best_protective_trait")
    op.drop_column("disease_annotation", "protective_variant_count")
