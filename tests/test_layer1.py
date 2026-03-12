"""Tests for Layer 1 — sequence divergence pipeline.

Uses fixture dataset (3 species, 5 genes). Runs in under 5 minutes on any hardware.
"""

import pytest


# ---------------------------------------------------------------------------
# Alignment tests
# ---------------------------------------------------------------------------


class TestAlignment:
    def test_sequence_identity_identical(self):
        from pipeline.layer1_sequence.alignment import calculate_sequence_identity

        seq = "MEPQSDPSVEPPL"
        assert calculate_sequence_identity(seq, seq) == 100.0

    def test_sequence_identity_all_different(self):
        from pipeline.layer1_sequence.alignment import calculate_sequence_identity

        seq_a = "AAAAAAAAAA"
        seq_b = "CCCCCCCCCC"
        assert calculate_sequence_identity(seq_a, seq_b) == 0.0

    def test_sequence_identity_half_different(self):
        from pipeline.layer1_sequence.alignment import calculate_sequence_identity

        seq_a = "AAAAAACCCC"
        seq_b = "AAAAAAGGGG"
        result = calculate_sequence_identity(seq_a, seq_b)
        assert result == pytest.approx(60.0)

    def test_sequence_identity_with_gaps(self):
        from pipeline.layer1_sequence.alignment import calculate_sequence_identity

        # Gaps in both positions are excluded
        seq_a = "AAAA--AA"
        seq_b = "AAAA--CC"
        result = calculate_sequence_identity(seq_a, seq_b)
        assert result == pytest.approx(66.67, rel=1e-2)

    def test_sequence_identity_raises_on_different_length(self):
        from pipeline.layer1_sequence.alignment import calculate_sequence_identity

        with pytest.raises(ValueError, match="aligned"):
            calculate_sequence_identity("AAA", "AAAA")

    def test_align_orthogroup_requires_two_sequences(self):
        from pipeline.layer1_sequence.alignment import align_orthogroup

        result = align_orthogroup("OG0001", {"human|P1": "MEPQ"})
        assert result is None


# ---------------------------------------------------------------------------
# Divergence tests
# ---------------------------------------------------------------------------


class TestDivergence:
    def test_score_window_identical(self):
        from pipeline.layer1_sequence.divergence import score_window

        assert score_window("AAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAA") == 0.0

    def test_score_window_all_different(self):
        from pipeline.layer1_sequence.divergence import score_window

        result = score_window("AAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCC")
        assert result == pytest.approx(1.0)

    def test_score_window_partial(self):
        from pipeline.layer1_sequence.divergence import score_window

        animal = "AAAAAACCCCCAAAA"   # 5 diffs
        human  = "AAAAAAAAAAAAAAAA"  # 14 chars, but only 15 compared
        # rewrite: 5 out of 15 different
        animal = "AAAAAACCCCCAAAA"
        human  = "AAAAAAAAAAAAAAA"
        result = score_window(animal, human)
        assert result == pytest.approx(5 / 15)

    def test_score_window_gap_handling(self):
        from pipeline.layer1_sequence.divergence import score_window

        # Gap-only column excluded
        animal = "AAAA-AAAAAAAAAA"
        human  = "AAAA-AAAAAAAAAA"
        result = score_window(animal, human)
        assert result == 0.0

    def test_extract_motifs_finds_divergent_windows(self):
        from pipeline.layer1_sequence.divergence import extract_divergent_motifs

        # Human and animal differ by 5 AAs in window at position 0
        human_seq   = "A" * 60
        animal_seq  = "C" * 15 + "A" * 45   # first window is all different

        aligned_seqs = {
            "human|P1":           human_seq,
            "naked_mole_rat|P2":  animal_seq,
        }
        motifs = extract_divergent_motifs("OG0001", aligned_seqs, min_divergence=0.1)
        assert len(motifs) > 0
        assert motifs[0]["divergence_score"] == pytest.approx(1.0)
        assert motifs[0]["species_id"] == "naked_mole_rat"

    def test_extract_motifs_returns_empty_below_threshold(self):
        from pipeline.layer1_sequence.divergence import extract_divergent_motifs

        human_seq  = "A" * 60
        animal_seq = "A" * 59 + "C"   # Only 1 position different — below 15% threshold

        aligned_seqs = {
            "human|P1":           human_seq,
            "naked_mole_rat|P2":  animal_seq,
        }
        motifs = extract_divergent_motifs("OG0001", aligned_seqs, min_divergence=0.15)
        assert len(motifs) == 0

    def test_extract_motifs_skips_gap_heavy_windows(self):
        from pipeline.layer1_sequence.divergence import extract_divergent_motifs

        human_seq  = "A" * 60
        # First 8 of 15 positions are gaps (>50%) — first window should be skipped
        animal_seq = "-" * 8 + "C" * 7 + "A" * 45

        aligned_seqs = {
            "human|P1":           human_seq,
            "naked_mole_rat|P2":  animal_seq,
        }
        motifs = extract_divergent_motifs("OG0001", aligned_seqs, min_divergence=0.0)
        # The window starting at 0 has 8 gaps (>50%) — should be skipped
        # Windows at step 5 (pos 5) also have >50% gaps — also skipped
        # First valid window is at pos 10 (8 gaps in first 8 positions, step 10 has only 3 gaps)
        for m in motifs:
            # No window where >50% are gaps should appear
            gap_count = m["animal_seq"].count("-")
            assert gap_count <= len(m["animal_seq"]) * 0.5


# ---------------------------------------------------------------------------
# Download validation tests
# ---------------------------------------------------------------------------


class TestDownloadValidation:
    def test_count_sequences_empty_file(self, tmp_path):
        from pipeline.layer1_sequence.download import _count_sequences

        empty = tmp_path / "empty.faa"
        empty.write_text("")
        assert _count_sequences(empty) == 0

    def test_count_sequences_valid_fasta(self, tmp_path):
        from pipeline.layer1_sequence.download import _count_sequences

        fasta = tmp_path / "test.faa"
        fasta.write_text(">seq1\nMEPQ\n>seq2\nACGT\n>seq3\nLLLL\n")
        assert _count_sequences(fasta) == 3

    def test_validate_proteome_passes(self, tmp_path):
        from pipeline.layer1_sequence.download import validate_proteome

        fasta = tmp_path / "test.faa"
        records = "\n".join(f">seq{i}\nMEPQAAAAAAAAAAAAAAAAAAAAA\n" for i in range(200))
        fasta.write_text(records)
        assert validate_proteome(fasta, min_proteins=100) is True

    def test_validate_proteome_fails_below_min(self, tmp_path):
        from pipeline.layer1_sequence.download import validate_proteome

        fasta = tmp_path / "test.faa"
        fasta.write_text(">seq1\nMEPQ\n>seq2\nACGT\n")
        assert validate_proteome(fasta, min_proteins=100) is False

    def test_reheader_fasta(self, tmp_path):
        from pipeline.layer1_sequence.download import reheader_fasta

        fasta = tmp_path / "input.faa"
        fasta.write_text(">NP_000123.1 some protein\nMEPQ\n>NP_000124.1 another\nACGT\n")
        out = reheader_fasta(fasta, "naked_mole_rat", tmp_path / "out.faa")
        content = out.read_text()
        assert ">naked_mole_rat|NP_000123.1" in content
        assert ">naked_mole_rat|NP_000124.1" in content
