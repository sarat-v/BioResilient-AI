"""Tests for Layer 3 disease annotation modules (with mocked or offline logic)."""


from pipeline.scoring import disease_score
from db.models import DiseaseAnnotation

_GENE_ID = "a0eebc99-9c0b-4ef8-bb6d-6bb9bd380a11"


def test_disease_score_none():
    assert disease_score(None) == 0.0


def test_disease_score_opentargets_only():
    ann = DiseaseAnnotation(gene_id=_GENE_ID)
    ann.opentargets_score = 0.8
    ann.gwas_pvalue = None
    ann.gnomad_pli = None
    s = disease_score(ann)
    assert s > 0.2 and s <= 1.0


def test_disease_score_full():
    ann = DiseaseAnnotation(gene_id=_GENE_ID)
    ann.opentargets_score = 0.5
    ann.gwas_pvalue = 1e-8
    ann.gnomad_pli = 0.95
    s = disease_score(ann)
    assert s > 0.5 and s <= 1.0
