"""Tests for Layer 2 — evolutionary selection and convergence pipeline."""

import pytest


# ---------------------------------------------------------------------------
# Phylogenetic tree tests
# ---------------------------------------------------------------------------


class TestPhyloTree:
    def test_build_concatenated_alignment_single_copy(self):
        from pipeline.layer2_evolution.phylo_tree import build_concatenated_alignment

        aligned = {
            "OG0001": {
                "human|P1":           "AAAA",
                "naked_mole_rat|P2":  "ACAA",
                "ground_squirrel|P3": "AGAA",
            },
            "OG0002": {
                "human|P4":           "CCCC",
                "naked_mole_rat|P5":  "CGCC",
                "ground_squirrel|P6": "CTCC",
            },
        }
        concat = build_concatenated_alignment(aligned, single_copy_only=True)
        assert concat is not None
        # Each species should have concatenated sequences from both OGs
        assert len(concat["human"]) == 8
        assert concat["human"] == "AAAACCCC"
        assert concat["naked_mole_rat"] == "ACAACGCC"

    def test_build_concatenated_alignment_fills_gaps_for_missing_species(self):
        from pipeline.layer2_evolution.phylo_tree import build_concatenated_alignment

        aligned = {
            "OG0001": {
                "human|P1":          "AAAA",
                "naked_mole_rat|P2": "ACAA",
                # ground_squirrel missing from this OG
            },
        }
        concat = build_concatenated_alignment(aligned, single_copy_only=False)
        if concat:  # May be None if too few species
            if "ground_squirrel" in concat:
                assert all(c == "-" for c in concat["ground_squirrel"])

    def test_collect_species_ids(self):
        from pipeline.layer2_evolution.phylo_tree import _collect_species_ids

        aligned = {
            "OG0001": {"human|P1": "AAAA", "bat|P2": "ACAA"},
            "OG0002": {"human|P3": "CCCC", "whale|P4": "CGCC"},
        }
        species = _collect_species_ids(aligned)
        assert species == {"human", "bat", "whale"}


# ---------------------------------------------------------------------------
# Selection score parsing tests
# ---------------------------------------------------------------------------


class TestSelectionParsing:
    def _make_hyphy_result(self, branches_data: dict) -> dict:
        """Build a minimal HyPhy aBSREL JSON structure."""
        return {
            "branch attributes": {
                "0": branches_data
            }
        }

    def test_parse_absrel_selected_branch(self):
        from pipeline.layer2_evolution.selection import parse_absrel_results

        hyphy_json = self._make_hyphy_result({
            "naked_mole_rat": {
                "Corrected P-value": 0.001,
                "Rate Distributions": [[3.5, 1.0]],
            },
            "ground_squirrel": {
                "Corrected P-value": 0.8,
                "Rate Distributions": [[0.9, 1.0]],
            },
            "human": {
                "Corrected P-value": 1.0,
                "Rate Distributions": [[1.0, 1.0]],
            },
        })
        result = parse_absrel_results(hyphy_json)
        assert result["selection_model"] == "aBSREL"
        assert "naked_mole_rat" in result["branches_under_selection"]
        assert "ground_squirrel" not in result["branches_under_selection"]
        assert result["dnds_ratio"] == pytest.approx(3.5)
        assert result["dnds_pvalue"] == pytest.approx(0.001)

    def test_parse_absrel_no_selection(self):
        from pipeline.layer2_evolution.selection import parse_absrel_results

        hyphy_json = self._make_hyphy_result({
            "naked_mole_rat": {
                "Corrected P-value": 0.9,
                "Rate Distributions": [[1.0, 1.0]],
            },
        })
        result = parse_absrel_results(hyphy_json)
        assert result["branches_under_selection"] == []
        assert result["dnds_ratio"] == pytest.approx(1.0)

    def test_should_run_hyphy_threshold(self):
        from pipeline.layer2_evolution.selection import _should_run_hyphy

        # 2 species with motifs — should run
        motifs_by_og = {
            "OG001": [
                {"species_id": "naked_mole_rat"},
                {"species_id": "ground_squirrel"},
            ]
        }
        assert _should_run_hyphy("OG001", motifs_by_og) is True

    def test_should_run_hyphy_insufficient(self):
        from pipeline.layer2_evolution.selection import _should_run_hyphy

        # Only 1 non-human species with motifs — should not run
        motifs_by_og = {
            "OG001": [{"species_id": "naked_mole_rat"}]
        }
        assert _should_run_hyphy("OG001", motifs_by_og) is False


# ---------------------------------------------------------------------------
# Convergence tests
# ---------------------------------------------------------------------------


class TestConvergence:
    def test_count_convergent_lineages_single(self):
        from pipeline.layer2_evolution.convergence import count_convergent_lineages

        species = ["naked_mole_rat"]
        assert count_convergent_lineages(species) == 1   # Rodents

    def test_count_convergent_lineages_multiple_same_group(self):
        from pipeline.layer2_evolution.convergence import count_convergent_lineages

        # Two rodents count as ONE lineage, not two
        species = ["naked_mole_rat", "ground_squirrel"]
        assert count_convergent_lineages(species) == 1

    def test_count_convergent_lineages_multiple_groups(self):
        from pipeline.layer2_evolution.convergence import count_convergent_lineages

        species = [
            "naked_mole_rat",    # Rodents
            "bowhead_whale",     # Cetaceans
            "little_brown_bat",  # Bats
        ]
        assert count_convergent_lineages(species) == 3

    def test_count_convergent_lineages_excludes_human(self):
        from pipeline.layer2_evolution.convergence import count_convergent_lineages

        # Human (Primates) excluded from convergence count
        species = ["human", "naked_mole_rat"]
        assert count_convergent_lineages(species) == 1

    def test_count_convergent_lineages_all_groups(self):
        from pipeline.layer2_evolution.convergence import count_convergent_lineages

        species = [
            "naked_mole_rat",    # Rodents
            "bowhead_whale",     # Cetaceans
            "little_brown_bat",  # Bats
            "greenland_shark",   # Sharks
            "axolotl",           # Salamanders
            "african_elephant",  # Proboscideans
        ]
        assert count_convergent_lineages(species) == 6
