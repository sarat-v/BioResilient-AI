"""Add is_control to species; control_divergence_fraction to candidate_score.

Revision ID: 0018
Revises: 0017
Create Date: 2026-03-01
"""
from alembic import op
import sqlalchemy as sa

revision = "0018"
down_revision = "0017"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("species", sa.Column("is_control", sa.Boolean(), nullable=True, server_default="false"))
    op.add_column("candidate_score", sa.Column("control_divergence_fraction", sa.Float(), nullable=True))
    op.create_index("ix_candidate_score_control_div", "candidate_score", ["control_divergence_fraction"])


def downgrade() -> None:
    op.drop_index("ix_candidate_score_control_div", table_name="candidate_score")
    op.drop_column("candidate_score", "control_divergence_fraction")
    op.drop_column("species", "is_control")
