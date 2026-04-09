"""Phase 2 extended fields: gnomAD LOEUF, OT tractability, known drugs, pocket proximity.

Revision ID: 0024
Revises: 0023
Create Date: 2026-04-09
"""

from alembic import op
import sqlalchemy as sa

revision = "0024"
down_revision = "0023"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # disease_annotation: gnomAD LOEUF + OpenTargets extended fields
    op.add_column("disease_annotation", sa.Column("gnomad_loeuf", sa.Float(), nullable=True))
    op.add_column("disease_annotation", sa.Column("known_drug_name", sa.String(), nullable=True))
    op.add_column("disease_annotation", sa.Column("known_drug_phase", sa.Integer(), nullable=True))
    op.add_column("disease_annotation", sa.Column("ot_safety_liability", sa.String(), nullable=True))

    # drug_target: tractability modalities + convergent pocket proximity
    op.add_column("drug_target", sa.Column("tractability_sm", sa.Boolean(), nullable=True))
    op.add_column("drug_target", sa.Column("tractability_ab", sa.Boolean(), nullable=True))
    op.add_column("drug_target", sa.Column("tractability_protac", sa.Boolean(), nullable=True))
    op.add_column("drug_target", sa.Column("convergent_pocket_proximal", sa.Boolean(), nullable=True))


def downgrade() -> None:
    op.drop_column("drug_target", "convergent_pocket_proximal")
    op.drop_column("drug_target", "tractability_protac")
    op.drop_column("drug_target", "tractability_ab")
    op.drop_column("drug_target", "tractability_sm")
    op.drop_column("disease_annotation", "ot_safety_liability")
    op.drop_column("disease_annotation", "known_drug_phase")
    op.drop_column("disease_annotation", "known_drug_name")
    op.drop_column("disease_annotation", "gnomad_loeuf")
