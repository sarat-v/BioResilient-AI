"""Tests for the composite scoring module."""

import pytest


class TestScoringFunctions:
    def test_convergence_score_zero_lineages(self):
        from pipeline.scoring import convergence_score

        assert convergence_score(0) == 0.0
        assert convergence_score(None) == 0.0

    def test_convergence_score_increases_with_lineages(self):
        from pipeline.scoring import convergence_score

        scores = [convergence_score(n) for n in range(7)]
        for i in range(1, len(scores)):
            assert scores[i] > scores[i - 1], f"Score should increase: {scores}"

    def test_convergence_score_bounded(self):
        from pipeline.scoring import convergence_score

        for n in range(10):
            s = convergence_score(n)
            assert 0.0 <= s <= 1.0

    def test_selection_score_neutral_evolution(self):
        from pipeline.scoring import selection_score

        # dN/dS = 1.0 (neutral), high p-value → low score
        s = selection_score(1.0, 0.5)
        assert s < 0.5

    def test_selection_score_strong_positive_selection(self):
        from pipeline.scoring import selection_score

        # dN/dS = 5.0 (max normalised), very significant (p < 0.001 → -log10 = 3, norm = 0.3)
        # score = 0.5 * 1.0 + 0.5 * 0.3 = 0.65+; definitely higher than neutral
        s = selection_score(5.0, 0.0001)
        assert s > 0.6
        # Also verify that higher dN/dS + lower p-value scores higher than neutral
        neutral = selection_score(1.0, 0.5)
        assert s > neutral

    def test_selection_score_purifying(self):
        from pipeline.scoring import selection_score

        # dN/dS = 0.1 (strong purifying), non-significant
        s = selection_score(0.1, 0.9)
        assert s < 0.2

    def test_selection_score_none_returns_zero(self):
        from pipeline.scoring import selection_score

        assert selection_score(None, None) == 0.0
        assert selection_score(2.0, None) == 0.0

    def test_selection_score_bounded(self):
        from pipeline.scoring import selection_score

        for dnds in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]:
            for pval in [0.001, 0.05, 0.1, 0.5, 1.0]:
                s = selection_score(dnds, pval)
                assert 0.0 <= s <= 1.0

    def test_composite_score_weighted(self):
        from pipeline.scoring import composite_score

        sub_scores = {
            "convergence": 0.8,
            "selection": 0.6,
            "expression": 0.4,
            "disease": 0.0,
            "druggability": 0.0,
            "safety": 0.0,
        }
        weights = {"convergence": 0.4, "selection": 0.35, "expression": 0.25,
                   "disease": 0.0, "druggability": 0.0, "safety": 0.0}
        result = composite_score(sub_scores, weights)
        expected = (0.8 * 0.4 + 0.6 * 0.35 + 0.4 * 0.25) / (0.4 + 0.35 + 0.25)
        assert result == pytest.approx(expected, rel=1e-4)

    def test_composite_score_zero_weights(self):
        from pipeline.scoring import composite_score

        result = composite_score({"convergence": 1.0}, {"convergence": 0.0})
        assert result == 0.0

    def test_assign_tier_tier1(self):
        from pipeline.scoring import assign_tier

        assert assign_tier(0.75, {"tier1": 0.70, "tier2": 0.40}) == "Tier1"

    def test_assign_tier_tier2(self):
        from pipeline.scoring import assign_tier

        assert assign_tier(0.55, {"tier1": 0.70, "tier2": 0.40}) == "Tier2"

    def test_assign_tier_tier3(self):
        from pipeline.scoring import assign_tier

        assert assign_tier(0.30, {"tier1": 0.70, "tier2": 0.40}) == "Tier3"

    def test_assign_tier_boundary_tier1(self):
        from pipeline.scoring import assign_tier

        assert assign_tier(0.70, {"tier1": 0.70, "tier2": 0.40}) == "Tier1"

    def test_assign_tier_boundary_tier2(self):
        from pipeline.scoring import assign_tier

        assert assign_tier(0.40, {"tier1": 0.70, "tier2": 0.40}) == "Tier2"


class TestScoringWeights:
    def test_phase1_weights_sum_to_one(self):
        from pipeline.config import get_scoring_weights

        weights = get_scoring_weights("phase1")
        total = sum(weights.values())
        assert total == pytest.approx(1.0, rel=1e-4)

    def test_phase1_disease_weight_zero(self):
        from pipeline.config import get_scoring_weights

        weights = get_scoring_weights("phase1")
        assert weights.get("disease", 0.0) == 0.0

    def test_phase2_weights_sum_to_one(self):
        from pipeline.config import get_scoring_weights

        weights = get_scoring_weights("phase2")
        total = sum(weights.values())
        assert total == pytest.approx(1.0, rel=1e-4)


class TestAPIModels:
    """Test Pydantic model validation for the API layer."""

    def test_candidate_list_item_fields(self):
        from api.routes.candidates import CandidateListItem

        item = CandidateListItem(
            gene_id="abc-123",
            gene_symbol="TP53",
            human_protein="P04637",
            composite_score=0.85,
            tier="Tier1",
            convergence_score=0.9,
            selection_score=0.8,
            expression_score=0.7,
            updated_at=None,
        )
        assert item.tier == "Tier1"
        assert item.composite_score == pytest.approx(0.85)

    def test_score_breakdown_model(self):
        from api.routes.scores import ScoreBreakdown

        breakdown = ScoreBreakdown(
            gene_id="abc",
            gene_symbol="BRCA1",
            tier="Tier2",
            composite_score=0.55,
            sub_scores={
                "convergence": 0.7,
                "selection": 0.5,
                "expression": 0.3,
                "disease": 0.0,
                "druggability": 0.0,
                "safety": 0.0,
            },
            evolution={"dnds_ratio": 2.1, "dnds_pvalue": 0.02},
        )
        assert breakdown.sub_scores["convergence"] == pytest.approx(0.7)
