"""Add trait_id to CandidateScore for multi-trait support

Revision ID: 0004
Revises: 0003
Create Date: 2026-03-01
"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op

revision: str = "0004"
down_revision: Union[str, None] = "0003"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column(
        "candidate_score",
        sa.Column("trait_id", sa.String(), nullable=True, server_default=""),
    )
    op.execute("UPDATE candidate_score SET trait_id = '' WHERE trait_id IS NULL")
    op.alter_column(
        "candidate_score",
        "trait_id",
        nullable=False,
        server_default="",
    )
    op.drop_constraint("candidate_score_pkey", "candidate_score", type_="primary")
    op.create_primary_key("candidate_score_pkey", "candidate_score", ["gene_id", "trait_id"])


def downgrade() -> None:
    op.drop_constraint("candidate_score_pkey", "candidate_score", type_="primary")
    op.create_primary_key("candidate_score_pkey", "candidate_score", ["gene_id"])
    op.alter_column("candidate_score", "trait_id", nullable=True)
    op.drop_column("candidate_score", "trait_id")
