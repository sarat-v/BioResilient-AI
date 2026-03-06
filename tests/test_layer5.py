"""Tests for Layer 5 gene therapy modules."""

import pytest

from pipeline.layer5_gene_therapy.aav import AAV_PACKAGING_LIMIT_BP, AAV_TROPISM
from pipeline.layer5_gene_therapy.crispr import count_pam_sites, run_crispor_offtarget


def test_aav_packaging_limit():
    assert AAV_PACKAGING_LIMIT_BP == 4700


def test_aav_tropism_keys():
    assert "AAV9" in AAV_TROPISM
    assert "AAV8" in AAV_TROPISM
    assert "CNS" in AAV_TROPISM["AAV9"] or "liver" in AAV_TROPISM["AAV8"]


def test_count_pam_sites():
    # NGG pattern: 20 bp + GG
    seq = "A" * 20 + "GG" + "T" * 20 + "GG"
    assert count_pam_sites(seq) >= 2
    assert count_pam_sites("") == 0


def test_run_crispor_offtarget_returns_tuple():
    sites, risk = run_crispor_offtarget("")
    assert sites == 0
    assert risk == "unknown"
    sites2, risk2 = run_crispor_offtarget("A" * 100 + "GG")
    assert isinstance(sites2, int)
    assert risk2 in ("low", "medium", "high", "unknown")
