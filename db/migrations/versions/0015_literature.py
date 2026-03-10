"""Add literature validation fields to disease_annotation.

Revision ID: 0015
Revises: 0014
Create Date: 2026-03-10
"""
from alembic import op
import sqlalchemy as sa

revision = "0015"
down_revision = "0014"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("disease_annotation", sa.Column("lit_score", sa.Float(), nullable=True))
    op.add_column("disease_annotation", sa.Column("lit_pmid_count", sa.Integer(), nullable=True))


def downgrade() -> None:
    op.drop_column("disease_annotation", "lit_pmid_count")
    op.drop_column("disease_annotation", "lit_score")
