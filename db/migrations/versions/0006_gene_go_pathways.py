"""Add go_terms and pathway_ids to Gene for pathway filtering

Revision ID: 0006
Revises: 0005
Create Date: 2026-03-01
"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects import postgresql

revision: str = "0006"
down_revision: Union[str, None] = "0005"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column(
        "gene",
        sa.Column("go_terms", postgresql.ARRAY(sa.String()), nullable=True),
    )
    op.add_column(
        "gene",
        sa.Column("pathway_ids", postgresql.ARRAY(sa.String()), nullable=True),
    )


def downgrade() -> None:
    op.drop_column("gene", "pathway_ids")
    op.drop_column("gene", "go_terms")
