"""Add meme_qvalue (BH-corrected FDR) and RELAX branch acceleration to evolution_score.

Revision ID: 0017
Revises: 0016
Create Date: 2026-03-01
"""
from alembic import op
import sqlalchemy as sa

revision = "0017"
down_revision = "0016"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column("evolution_score", sa.Column("meme_qvalue", sa.Float(), nullable=True))
    op.add_column("evolution_score", sa.Column("relax_k", sa.Float(), nullable=True))
    op.add_column("evolution_score", sa.Column("relax_pvalue", sa.Float(), nullable=True))
    op.create_index("ix_evolution_score_meme_qvalue", "evolution_score", ["meme_qvalue"])


def downgrade() -> None:
    op.drop_index("ix_evolution_score_meme_qvalue", table_name="evolution_score")
    op.drop_column("evolution_score", "relax_pvalue")
    op.drop_column("evolution_score", "relax_k")
    op.drop_column("evolution_score", "meme_qvalue")
