"""Add unique constraints for idempotent reruns.

Prevents duplicate rows in divergent_motif and expression_result when
pipeline steps are rerun. Existing duplicates must be cleaned first.

Revision ID: 0022
Revises: 0021
Create Date: 2026-03-24
"""
from alembic import op

revision = "0022"
down_revision = "0021"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Clean up existing duplicates in divergent_motif before adding constraint
    op.execute("""
        DELETE FROM divergent_motif
        WHERE id NOT IN (
            SELECT DISTINCT ON (ortholog_id, start_pos) id
            FROM divergent_motif
            ORDER BY ortholog_id, start_pos,
                CASE WHEN motif_direction IS NOT NULL THEN 0 ELSE 1 END,
                id DESC
        )
    """)

    # Clean up existing duplicates in expression_result before adding constraint
    op.execute("""
        DELETE FROM expression_result
        WHERE id NOT IN (
            SELECT DISTINCT ON (gene_id, geo_accession, COALESCE(comparison, '')) id
            FROM expression_result
            ORDER BY gene_id, geo_accession, COALESCE(comparison, ''), id DESC
        )
    """)

    op.create_unique_constraint(
        "uq_motif_ortholog_pos", "divergent_motif",
        ["ortholog_id", "start_pos"],
    )
    op.create_unique_constraint(
        "uq_expression_result", "expression_result",
        ["gene_id", "geo_accession", "comparison"],
    )


def downgrade() -> None:
    op.drop_constraint("uq_expression_result", "expression_result", type_="unique")
    op.drop_constraint("uq_motif_ortholog_pos", "divergent_motif", type_="unique")
