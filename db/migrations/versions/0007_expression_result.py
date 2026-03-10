"""Add ExpressionResult table for expression evidence traceability

Revision ID: 0007
Revises: 0006
Create Date: 2026-03-01
"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op

revision: str = "0007"
down_revision: Union[str, None] = "0006"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        "expression_result",
        sa.Column("id", sa.String(), nullable=False),
        sa.Column("gene_id", sa.String(), nullable=False),
        sa.Column("geo_accession", sa.String(), nullable=False),
        sa.Column("comparison", sa.String(), nullable=True),
        sa.Column("log2fc", sa.Float(), nullable=True),
        sa.Column("padj", sa.Float(), nullable=True),
        sa.Column("n_samples", sa.Integer(), nullable=True),
        sa.ForeignKeyConstraint(["gene_id"], ["gene.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("id"),
    )
    op.create_index("ix_expression_result_gene_id", "expression_result", ["gene_id"])


def downgrade() -> None:
    op.drop_index("ix_expression_result_gene_id", table_name="expression_result")
    op.drop_table("expression_result")
