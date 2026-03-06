"""Add RegulatoryDivergence table and regulatory_score column on CandidateScore

Revision ID: 0002
Revises: 0001
Create Date: 2026-03-01
"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op
from sqlalchemy.dialects import postgresql

revision: str = "0002"
down_revision: Union[str, None] = "0001"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column(
        "candidate_score",
        sa.Column("regulatory_score", sa.Float(), nullable=True, server_default="0.0"),
    )

    op.create_table(
        "regulatory_divergence",
        sa.Column("id", sa.String(), nullable=False),
        sa.Column("gene_id", sa.String(), nullable=False),
        sa.Column("species_id", sa.String(), nullable=False),
        sa.Column("promoter_divergence", sa.Float(), nullable=True),
        sa.Column("expression_log2fc", sa.Float(), nullable=True),
        sa.Column("lineage_count", sa.Integer(), nullable=True),
        sa.Column("regulatory_score", sa.Float(), nullable=True),
        sa.ForeignKeyConstraint(["gene_id"], ["gene.id"], ondelete="CASCADE"),
        sa.ForeignKeyConstraint(["species_id"], ["species.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("gene_id", "species_id", name="uq_regulatory_divergence_gene_species"),
    )

    op.create_index("ix_regulatory_divergence_gene_id", "regulatory_divergence", ["gene_id"])


def downgrade() -> None:
    op.drop_index("ix_regulatory_divergence_gene_id", table_name="regulatory_divergence")
    op.drop_table("regulatory_divergence")
    op.drop_column("candidate_score", "regulatory_score")
