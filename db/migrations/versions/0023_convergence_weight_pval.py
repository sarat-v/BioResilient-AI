"""Add convergence_weight and convergence_pval columns to evolution_score.

Fix 3: convergence_weight stores the phylogenetically-weighted convergence score,
previously incorrectly written into phylop_score when that column was NULL.

Fix 2: convergence_pval stores the empirical p-value from the permutation-based
null model (200-iteration species-label shuffle in Step 7a).

Revision ID: 0023
Revises: 0022
Create Date: 2026-04-07
"""
from alembic import op
import sqlalchemy as sa

revision = "0023"
down_revision = "0022"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column(
        "evolution_score",
        sa.Column("convergence_weight", sa.Float(), nullable=True),
    )
    op.add_column(
        "evolution_score",
        sa.Column("convergence_pval", sa.Float(), nullable=True),
    )

    # Back-fill convergence_weight from phylop_score where phylop_score currently
    # holds the convergence weight (i.e. where convergence_count > 0 and the
    # UCSC PhyloP enrichment step has not yet been run on that gene).
    # Because we cannot distinguish the two sources reliably after-the-fact, we
    # copy the value so that existing scores are preserved until step 7a re-runs.
    op.execute("""
        UPDATE evolution_score
        SET convergence_weight = phylop_score
        WHERE convergence_count > 0
          AND phylop_score IS NOT NULL
          AND convergence_weight IS NULL
    """)


def downgrade() -> None:
    op.drop_column("evolution_score", "convergence_pval")
    op.drop_column("evolution_score", "convergence_weight")
