"""Add is_one_to_one to ortholog table.

Revision ID: 0008
Revises: 0007
Create Date: 2026-03-10
"""
from alembic import op
import sqlalchemy as sa

revision = "0008"
down_revision = "0007"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column(
        "ortholog",
        sa.Column("is_one_to_one", sa.Boolean(), nullable=False, server_default=sa.true()),
    )
    op.create_index("ix_ortholog_is_one_to_one", "ortholog", ["is_one_to_one"])


def downgrade() -> None:
    op.drop_index("ix_ortholog_is_one_to_one", table_name="ortholog")
    op.drop_column("ortholog", "is_one_to_one")
