"""Add domain annotation and AlphaMissense consequence_score to divergent_motif.

Revision ID: 0009
Revises: 0008
Create Date: 2026-03-10
"""
from alembic import op
import sqlalchemy as sa

revision = "0009"
down_revision = "0008"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # U2: Pfam/InterPro domain columns
    op.add_column("divergent_motif", sa.Column("domain_name", sa.String(), nullable=True))
    op.add_column(
        "divergent_motif",
        sa.Column("in_functional_domain", sa.Boolean(), nullable=False, server_default=sa.false()),
    )
    # U3: AlphaMissense consequence score
    op.add_column("divergent_motif", sa.Column("consequence_score", sa.Float(), nullable=True))
    # Index for fast filtering on functional domain motifs
    op.create_index("ix_divergent_motif_functional", "divergent_motif", ["in_functional_domain"])


def downgrade() -> None:
    op.drop_index("ix_divergent_motif_functional", table_name="divergent_motif")
    op.drop_column("divergent_motif", "consequence_score")
    op.drop_column("divergent_motif", "in_functional_domain")
    op.drop_column("divergent_motif", "domain_name")
