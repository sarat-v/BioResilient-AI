"""Tests for Layer 6 safety modules."""

import pytest

from pipeline.scoring import safety_score
from db.models import SafetyFlag

_GENE_ID = "a0eebc99-9c0b-4ef8-bb6d-6bb9bd380a11"


def test_safety_score_none():
    assert safety_score(None) == 1.0


def test_safety_score_hub_risk():
    sf = SafetyFlag(gene_id=_GENE_ID)
    sf.hub_risk = True
    sf.is_essential = False
    sf.family_size = 10
    s = safety_score(sf)
    assert s < 1.0
    assert s >= 0.0


def test_safety_score_essential():
    sf = SafetyFlag(gene_id=_GENE_ID)
    sf.hub_risk = False
    sf.is_essential = True
    sf.family_size = 5
    s = safety_score(sf)
    assert s < 1.0


def test_safety_score_large_family():
    sf = SafetyFlag(gene_id=_GENE_ID)
    sf.hub_risk = False
    sf.is_essential = False
    sf.family_size = 200
    s = safety_score(sf)
    assert s < 1.0
