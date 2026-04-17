"""Add BH-FDR adjusted_pvalue column to pathway_convergence table.

Multiple-testing correction for pathway enrichment analysis (Step 11d).
Without this column, pathway p-values were raw hypergeometric values with
no adjustment for the number of pathways tested (typically 10-100), which
inflates false-positive pathway hits.

The BH-adjusted q-value is computed in pathway_convergence.py using the
shared apply_bh_correction() utility and stored here for downstream consumers
(API endpoint GET /research/pathway-convergence, docs, and visualization).

Reference: Benjamini & Hochberg (1995) J Royal Stat Soc B 57:289-300.

Revision ID: 0031
Revises: 0030
Create Date: 2026-04-17
"""
from alembic import op
import sqlalchemy as sa

revision = "0031"
down_revision = "0030"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column(
        "pathway_convergence",
        sa.Column(
            "adjusted_pvalue",
            sa.Float(),
            nullable=True,
            comment=(
                "Benjamini-Hochberg FDR-adjusted p-value across all pathways tested "
                "in a single run. NULL for legacy rows created before migration 0031. "
                "Pathways with adjusted_pvalue < 0.05 are considered statistically enriched."
            ),
        ),
    )


def downgrade() -> None:
    op.drop_column("pathway_convergence", "adjusted_pvalue")
