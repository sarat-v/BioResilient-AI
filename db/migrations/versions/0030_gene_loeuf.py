"""Add loeuf column to gene table for variant direction classification.

Previously, variant_direction.py (step4d) downloaded gnomAD v4.1 LOEUF values
at runtime from S3/GCS into a throwaway Batch container's /tmp directory.
When that download failed or timed out, loeuf_map came back empty, causing all
ESM-destabilised + AM-high motifs to fall into 'functional_shift' instead of
'loss_of_function'. This resulted in loss_of_function = 0 for every gene.

Fix: store gnomAD LOEUF in the gene table once (populated by a backfill script
or step4b), then step4d reads it with a single SQL JOIN instead of a network
download. The loeuf column is nullable — NULL means 'unknown' and the
classifier already handles this gracefully by returning 'functional_shift'.

Revision ID: 0030
Revises: 0029
Create Date: 2026-04-16
"""
from alembic import op
import sqlalchemy as sa

revision = "0030"
down_revision = "0029"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column(
        "gene",
        sa.Column(
            "loeuf",
            sa.Float(),
            nullable=True,
            comment=(
                "gnomAD v4.1 loss-of-function observed/expected upper bound (LOEUF). "
                "Low values (<=0.5) indicate LoF intolerance. NULL = not available."
            ),
        ),
    )
    # Index for quick lookups by symbol when step4d reads loeuf_map
    op.create_index("ix_gene_loeuf", "gene", ["loeuf"])


def downgrade() -> None:
    op.drop_index("ix_gene_loeuf", table_name="gene")
    op.drop_column("gene", "loeuf")
