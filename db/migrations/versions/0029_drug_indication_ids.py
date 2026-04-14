"""Add known_drug_indication_ids to disease_annotation for phenotype-aware drug phase scoring.

Revision ID: 0029
Revises: 0028
Create Date: 2026-04-14
"""

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSON

revision = "0029"
down_revision = "0028"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column(
        "disease_annotation",
        sa.Column("known_drug_indication_ids", JSON, nullable=True),
    )


def downgrade() -> None:
    op.drop_column("disease_annotation", "known_drug_indication_ids")
