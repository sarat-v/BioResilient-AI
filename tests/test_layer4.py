"""Tests for Layer 4 druggability modules."""

import pytest

from pipeline.scoring import druggability_score
from pipeline.layer4_druggability.peptide import boman_index, is_synthesisable
from db.models import DrugTarget


def test_druggability_score_none():
    assert druggability_score(None) == 0.0


def test_druggability_score_pockets():
    dt = DrugTarget(gene_id="a0eebc99-9c0b-4ef8-bb6d-6bb9bd380a11")
    dt.pocket_count = 3
    dt.top_pocket_score = 0.6
    dt.chembl_target_id = "CHEMBL123"
    dt.existing_drugs = ["DrugA"]
    dt.druggability_tier = "B"
    s = druggability_score(dt)
    assert s > 0.3 and s <= 1.0


def test_boman_index():
    assert boman_index("") == 0.0
    assert boman_index("L") > 0
    assert boman_index("LLL") > boman_index("GGG")


def test_is_synthesisable_length():
    assert is_synthesisable("A" * 15, "A" * 15, length_ok=True) is True
    assert is_synthesisable("A" * 5, "A" * 5, length_ok=False) is False
    assert is_synthesisable("A" * 25, "A" * 25, length_ok=True) is False
