"""Add Compound and CompoundScore tables for Phase 3/4 lead discovery.

Revision ID: 0028
Revises: 0027
Create Date: 2026-04-14
"""

from alembic import op
import sqlalchemy as sa

revision = "0028"
down_revision = "0027"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "compound",
        sa.Column("id",          sa.String(), nullable=False, primary_key=True),
        sa.Column("gene_id",     sa.String(), sa.ForeignKey("gene.id", ondelete="CASCADE"),
                  nullable=False, index=True),
        sa.Column("smiles",             sa.String(), nullable=False),
        sa.Column("name",               sa.String(), nullable=True),
        sa.Column("source",             sa.String(), nullable=True),
        sa.Column("zinc_id",            sa.String(), nullable=True),
        sa.Column("chembl_id",          sa.String(), nullable=True),
        sa.Column("docking_score",      sa.Float(),  nullable=True),
        sa.Column("docking_method",     sa.String(), nullable=True),
        sa.Column("tanimoto_to_known",  sa.Float(),  nullable=True),
        sa.Column("nearest_known_drug", sa.String(), nullable=True),
        sa.Column("hit_score",          sa.Float(),  nullable=True),
        sa.Column("hit_tier",           sa.String(), nullable=True),
        sa.Column("is_purchasable",     sa.Boolean(), nullable=True),
        sa.Column("created_at",         sa.DateTime(), nullable=True),
    )

    op.create_table(
        "compound_score",
        sa.Column("id",          sa.String(), nullable=False, primary_key=True),
        sa.Column("compound_id", sa.String(), sa.ForeignKey("compound.id", ondelete="CASCADE"),
                  nullable=False, unique=True),
        sa.Column("mw",                   sa.Float(),   nullable=True),
        sa.Column("logp",                 sa.Float(),   nullable=True),
        sa.Column("hbd",                  sa.Integer(), nullable=True),
        sa.Column("hba",                  sa.Integer(), nullable=True),
        sa.Column("tpsa",                 sa.Float(),   nullable=True),
        sa.Column("lipinski_pass",        sa.Boolean(), nullable=True),
        sa.Column("bbb_permeable",        sa.Boolean(), nullable=True),
        sa.Column("cyp3a4_risk",          sa.Float(),   nullable=True),
        sa.Column("tox21_score",          sa.Float(),   nullable=True),
        sa.Column("herg_risk",            sa.Float(),   nullable=True),
        sa.Column("pains_alerts",         sa.Integer(), nullable=True),
        sa.Column("structural_alerts",    sa.Integer(), nullable=True),
        sa.Column("selectivity_score",    sa.Float(),   nullable=True),
        sa.Column("n_offtarget_similar",  sa.Integer(), nullable=True),
        sa.Column("sa_score",             sa.Float(),   nullable=True),
        sa.Column("sa_normalized",        sa.Float(),   nullable=True),
        sa.Column("zinc_purchasable",     sa.Boolean(), nullable=True),
        sa.Column("overall_compound_score", sa.Float(), nullable=True),
        sa.Column("compound_tier",        sa.String(),  nullable=True),
    )


def downgrade():
    op.drop_table("compound_score")
    op.drop_table("compound")
