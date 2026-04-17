"""Step 9b — Structural context for convergent positions.

Adds:
  - convergent_position_annotation table (per-motif structural annotations)
  - candidate_score.structural_score column

Revision ID: 0025
Revises: 0024
Create Date: 2026-04-07
"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op

revision: str = "0025"
down_revision: Union[str, None] = "0024"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # New table: convergent_position_annotation
    op.create_table(
        "convergent_position_annotation",
        sa.Column("id", sa.String(), nullable=False),
        sa.Column("gene_id", sa.String(), nullable=False),
        sa.Column("divergent_motif_id", sa.String(), nullable=False),
        sa.Column("uniprot_accession", sa.String(), nullable=True),
        sa.Column("protein_start_pos", sa.Integer(), nullable=True),
        sa.Column("protein_end_pos", sa.Integer(), nullable=True),
        sa.Column("convergent_lineage_count", sa.Integer(), nullable=True),
        sa.Column("am_score", sa.Float(), nullable=True),
        sa.Column("am_class", sa.String(), nullable=True),
        sa.Column("uniprot_feature", sa.String(), nullable=True),
        sa.Column("uniprot_feature_description", sa.String(), nullable=True),
        sa.Column("plddt_mean", sa.Float(), nullable=True),
        sa.Column("is_disordered", sa.Boolean(), nullable=True),
        sa.Column("pocket_distance_angstrom", sa.Float(), nullable=True),
        sa.Column("is_pocket_adjacent", sa.Boolean(), nullable=True),
        sa.Column("structural_context", sa.String(), nullable=True),
        sa.ForeignKeyConstraint(["gene_id"], ["gene.id"], ondelete="CASCADE"),
        sa.ForeignKeyConstraint(
            ["divergent_motif_id"], ["divergent_motif.id"], ondelete="CASCADE"
        ),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("gene_id", "divergent_motif_id", name="uq_cpa_gene_motif"),
    )
    op.create_index(
        "ix_cpa_gene_id", "convergent_position_annotation", ["gene_id"]
    )

    # New column on candidate_score
    op.add_column(
        "candidate_score",
        sa.Column("structural_score", sa.Float(), nullable=True, server_default="0.0"),
    )


def downgrade() -> None:
    op.drop_column("candidate_score", "structural_score")
    op.drop_index("ix_cpa_gene_id", table_name="convergent_position_annotation")
    op.drop_table("convergent_position_annotation")
