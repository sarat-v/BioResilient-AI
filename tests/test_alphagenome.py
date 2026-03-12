"""Tests for Track B — AlphaGenome regulatory divergence."""

import pytest

from db.models import RegulatoryDivergence, Gene
from pipeline.scoring import regulatory_score

_GENE_ID = "a0eebc99-9c0b-4ef8-bb6d-6bb9bd380a11"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _seed_gene(session):
    """Insert one Gene row; Species FK not enforced by SQLite without PRAGMA."""
    session.add(Gene(id=_GENE_ID, human_gene_id="G001", gene_symbol="TP53"))
    session.flush()


# ---------------------------------------------------------------------------
# regulatory_score() unit tests
# ---------------------------------------------------------------------------

def test_regulatory_score_no_rows(session):
    """Gene with no RegulatoryDivergence rows scores 0.0."""
    _seed_gene(session)
    score = regulatory_score(_GENE_ID, session)
    assert score == 0.0


def test_regulatory_score_single_row(session):
    """Single species divergence propagates to gene-level score."""
    _seed_gene(session)
    row = RegulatoryDivergence(
        gene_id=_GENE_ID,
        species_id="naked_mole_rat",
        promoter_divergence=0.45,
        expression_log2fc=1.2,
        lineage_count=1,
        regulatory_score=0.45,
    )
    session.add(row)
    session.flush()

    score = regulatory_score(_GENE_ID, session)
    # max_effect=0.45, lineage_count=1 → 0.45 + 0.05 = 0.50
    assert score == pytest.approx(0.50, abs=1e-4)


def test_regulatory_score_multiple_lineages(session):
    """Multiple independent lineages boost the gene-level score."""
    _seed_gene(session)
    rows = [
        RegulatoryDivergence(
            gene_id=_GENE_ID,
            species_id="naked_mole_rat",
            promoter_divergence=0.5,
            lineage_count=3,
            regulatory_score=0.5,
        ),
        RegulatoryDivergence(
            gene_id=_GENE_ID,
            species_id="bowhead_whale",
            promoter_divergence=0.4,
            lineage_count=3,
            regulatory_score=0.4,
        ),
        RegulatoryDivergence(
            gene_id=_GENE_ID,
            species_id="little_brown_bat",
            promoter_divergence=0.3,
            lineage_count=3,
            regulatory_score=0.3,
        ),
    ]
    session.add_all(rows)
    session.flush()

    score = regulatory_score(_GENE_ID, session)
    # max_effect=0.5, lineage_count=3 → 0.5 + 3*0.05 = 0.65
    assert score == pytest.approx(0.65, abs=1e-4)


def test_regulatory_score_capped_at_one(session):
    """Score is capped at 1.0 regardless of effect magnitude."""
    _seed_gene(session)
    row = RegulatoryDivergence(
        gene_id=_GENE_ID,
        species_id="naked_mole_rat",
        promoter_divergence=0.99,
        lineage_count=10,
        regulatory_score=0.99,
    )
    session.add(row)
    session.flush()

    score = regulatory_score(_GENE_ID, session)
    assert score == pytest.approx(1.0, abs=1e-4)


# ---------------------------------------------------------------------------
# AlphaGenome API key helper
# ---------------------------------------------------------------------------

def test_api_key_missing(monkeypatch):
    """When API key is absent, _get_api_key returns empty string."""
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    from pipeline.layer_regulatory.alphagenome import _get_api_key
    assert _get_api_key() == ""


def test_api_key_from_env(monkeypatch):
    """When ALPHAGENOME_API_KEY is set, _get_api_key returns it."""
    monkeypatch.setenv("ALPHAGENOME_API_KEY", "test-key-12345")
    from pipeline.layer_regulatory.alphagenome import _get_api_key
    assert _get_api_key() == "test-key-12345"


# ---------------------------------------------------------------------------
# Effect magnitude helpers
# ---------------------------------------------------------------------------

def test_alphagenome_effect_no_key():
    """Returns None when no API key is provided."""
    from pipeline.layer_regulatory.alphagenome import _alphagenome_effect
    result = _alphagenome_effect("ATCG" * 100, "ATCG" * 100, api_key="")
    assert result is None


def test_alphagenome_effect_empty_seq():
    """Returns None for empty sequences."""
    from pipeline.layer_regulatory.alphagenome import _alphagenome_effect
    result = _alphagenome_effect("", "", api_key="fake-key")
    assert result is None


# ---------------------------------------------------------------------------
# Lineage count helper
# ---------------------------------------------------------------------------

def test_compute_lineage_count_zero(session):
    """Returns 0 when no rows pass the effect threshold."""
    _seed_gene(session)
    from pipeline.layer_regulatory.alphagenome import _compute_lineage_count
    count = _compute_lineage_count(_GENE_ID, session)
    assert count == 0


def test_compute_lineage_count_two_lineages(session):
    """Correctly counts independent lineages (not species)."""
    _seed_gene(session)
    # Two rodent species → only 1 independent lineage; bat → second lineage
    session.add_all([
        RegulatoryDivergence(
            gene_id=_GENE_ID,
            species_id="naked_mole_rat",
            regulatory_score=0.8,
            lineage_count=0,
        ),
        RegulatoryDivergence(
            gene_id=_GENE_ID,
            species_id="damaraland_mole_rat",
            regulatory_score=0.7,
            lineage_count=0,
        ),
        RegulatoryDivergence(
            gene_id=_GENE_ID,
            species_id="little_brown_bat",
            regulatory_score=0.6,
            lineage_count=0,
        ),
    ])
    session.flush()

    from pipeline.layer_regulatory.alphagenome import _compute_lineage_count
    count = _compute_lineage_count(_GENE_ID, session)
    assert count == 2  # Rodents + Bats


def test_run_alphagenome_track_no_key(monkeypatch):
    """run_alphagenome_track returns 0 immediately when no API key is set."""
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    from pipeline.layer_regulatory.alphagenome import run_alphagenome_track
    result = run_alphagenome_track()
    assert result == 0
