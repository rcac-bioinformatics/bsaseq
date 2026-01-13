"""Tests for sliding window analysis functionality."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from bsaseq.analysis.windows import (
    calculate_g_statistic,
    calculate_windows,
    g_statistic_to_pvalue,
    tricube_weight,
)
from bsaseq.core.models import Variant, Window
from bsaseq.io.vcf import parse_vcf


class TestTricubeWeight:
    """Tests for tricube weight function."""

    def test_distance_zero_returns_one(self) -> None:
        """Test that distance 0 gives weight 1.0."""
        assert tricube_weight(0, 100) == 1.0

    def test_distance_equals_max_returns_zero(self) -> None:
        """Test that distance = max_distance gives weight 0.0."""
        assert tricube_weight(100, 100) == 0.0

    def test_distance_greater_than_max_returns_zero(self) -> None:
        """Test that distance > max_distance gives weight 0.0."""
        assert tricube_weight(150, 100) == 0.0

    def test_distance_half_max(self) -> None:
        """Test that distance = max/2 gives expected weight."""
        # (1 - 0.5^3)^3 = (1 - 0.125)^3 = 0.875^3 ≈ 0.6699
        weight = tricube_weight(50, 100)
        assert weight == pytest.approx(0.669921875, rel=1e-6)

    def test_zero_max_distance_returns_zero(self) -> None:
        """Test that max_distance = 0 gives weight 0.0."""
        assert tricube_weight(10, 0) == 0.0

    def test_negative_max_distance_returns_zero(self) -> None:
        """Test that negative max_distance gives weight 0.0."""
        assert tricube_weight(10, -100) == 0.0


class TestGStatistic:
    """Tests for G-statistic calculation."""

    def test_equal_frequencies_gives_zero(self) -> None:
        """Test that equal frequencies in both bulks give G ≈ 0."""
        g = calculate_g_statistic(
            obs_high_alt=50,
            obs_high_ref=50,
            obs_low_alt=50,
            obs_low_ref=50,
        )
        assert g == pytest.approx(0.0, abs=1e-6)

    def test_different_frequencies_gives_positive(self) -> None:
        """Test that different frequencies give positive G."""
        g = calculate_g_statistic(
            obs_high_alt=80,
            obs_high_ref=20,
            obs_low_alt=20,
            obs_low_ref=80,
        )
        assert g > 0

    def test_extreme_difference(self) -> None:
        """Test G-statistic with extreme frequency difference."""
        g = calculate_g_statistic(
            obs_high_alt=95,
            obs_high_ref=5,
            obs_low_alt=5,
            obs_low_ref=95,
        )
        # Should be a large value
        assert g > 100

    def test_all_zero_returns_zero(self) -> None:
        """Test that all-zero counts return 0."""
        g = calculate_g_statistic(0, 0, 0, 0)
        assert g == 0.0

    def test_one_cell_zero(self) -> None:
        """Test with one cell having zero counts."""
        g = calculate_g_statistic(
            obs_high_alt=50,
            obs_high_ref=50,
            obs_low_alt=0,
            obs_low_ref=100,
        )
        assert g > 0


class TestGStatisticToPvalue:
    """Tests for G-statistic to p-value conversion."""

    def test_zero_g_gives_pvalue_one(self) -> None:
        """Test that G=0 gives p-value = 1.0."""
        p = g_statistic_to_pvalue(0.0)
        assert p == 1.0

    def test_large_g_gives_small_pvalue(self) -> None:
        """Test that large G gives small p-value."""
        p = g_statistic_to_pvalue(100.0)
        assert p < 0.001

    def test_critical_value_3_84(self) -> None:
        """Test that G = 3.84 gives p ≈ 0.05 (chi-squared critical value)."""
        p = g_statistic_to_pvalue(3.84)
        assert p == pytest.approx(0.05, abs=0.01)


class TestWindowCalculation:
    """Tests for sliding window calculation."""

    def test_window_calculation_basic(self, sample_variants: list[Variant]) -> None:
        """Test basic window statistics on synthetic data."""
        windows = list(
            calculate_windows(
                sample_variants,
                window_size=10000,
                step_size=5000,
                min_variants=2,
            )
        )

        # Should have at least some windows
        assert len(windows) > 0

        # Check window properties
        for window in windows:
            assert isinstance(window, Window)
            assert window.n_variants >= 2
            assert window.start < window.end
            assert -1.0 <= window.mean_delta_af <= 1.0

    def test_windows_have_z_scores(self, sample_variants: list[Variant]) -> None:
        """Test that windows have z_scores computed."""
        windows = list(
            calculate_windows(
                sample_variants,
                window_size=10000,
                step_size=5000,
                min_variants=2,
            )
        )

        # All windows should have z_scores
        for window in windows:
            assert window.z_score is not None

    def test_windows_have_g_statistics(self, sample_variants: list[Variant]) -> None:
        """Test that windows have G-statistics and p-values."""
        windows = list(
            calculate_windows(
                sample_variants,
                window_size=10000,
                step_size=5000,
                min_variants=2,
            )
        )

        for window in windows:
            assert window.g_statistic >= 0
            assert 0 <= window.p_value <= 1

    def test_min_variants_filter(self, sample_variants: list[Variant]) -> None:
        """Test that windows with too few variants are filtered."""
        # Use small window with high min_variants - should get no windows
        windows = list(
            calculate_windows(
                sample_variants,
                window_size=100,  # Very small window
                step_size=50,
                min_variants=10,  # More than any window will have
            )
        )

        assert len(windows) == 0

    def test_empty_variants_returns_empty(self) -> None:
        """Test that empty variant list returns no windows."""
        windows = list(calculate_windows([], window_size=1000, step_size=500))
        assert len(windows) == 0

    def test_multiple_chromosomes(self, sample_variants: list[Variant]) -> None:
        """Test that windows are calculated for multiple chromosomes."""
        windows = list(
            calculate_windows(
                sample_variants,
                window_size=10000,
                step_size=2500,
                min_variants=1,
            )
        )

        # Should have windows from both chr1 and chr2
        chroms = {w.chrom for w in windows}
        assert "chr1" in chroms
        assert "chr2" in chroms

    def test_window_midpoint(self) -> None:
        """Test Window midpoint property."""
        window = Window(
            chrom="chr1",
            start=1000,
            end=2000,
            n_variants=5,
            mean_delta_af=0.3,
            median_delta_af=0.3,
            mean_af_high=0.6,
            mean_af_low=0.3,
            tricube_delta_af=0.3,
            g_statistic=10.0,
            p_value=0.001,
        )
        assert window.midpoint == 1500

    def test_window_size(self) -> None:
        """Test Window size property."""
        window = Window(
            chrom="chr1",
            start=1000,
            end=2000,
            n_variants=5,
            mean_delta_af=0.3,
            median_delta_af=0.3,
            mean_af_high=0.6,
            mean_af_low=0.3,
            tricube_delta_af=0.3,
            g_statistic=10.0,
            p_value=0.001,
        )
        assert window.size == 1001


class TestMultiSamplePooling:
    """Tests for multi-sample bulk pooling."""

    def test_multi_sample_pooling_comma_separated(
        self, multi_sample_vcf_path: Path
    ) -> None:
        """Test that multi-sample AD counts are correctly pooled (comma-separated)."""
        variants = list(
            parse_vcf(
                vcf_path=multi_sample_vcf_path,
                high_bulk="MUT1,MUT2",
                low_bulk="WT1,WT2",
                min_dp=10,
                max_dp=500,
            )
        )

        assert len(variants) == 3

        # Check first variant
        # MUT1: 5,45 + MUT2: 3,47 = 8 ref, 92 alt, total=100
        # AF = 92/100 = 0.92
        v1 = variants[0]
        assert v1.dp_high == 100  # 50 + 50 pooled
        assert v1.af_high == pytest.approx(0.92, rel=1e-6)
        assert v1.n_samples_high == 2
        assert v1.n_samples_low == 2

    def test_multi_sample_pooling_list(self, multi_sample_vcf_path: Path) -> None:
        """Test that multi-sample parsing works with list input."""
        variants = list(
            parse_vcf(
                vcf_path=multi_sample_vcf_path,
                high_bulk=["MUT1", "MUT2"],
                low_bulk=["WT1", "WT2"],
                min_dp=10,
                max_dp=500,
            )
        )

        assert len(variants) == 3
        # Same expected results as comma-separated
        v1 = variants[0]
        assert v1.dp_high == 100
        assert v1.af_high == pytest.approx(0.92, rel=1e-6)

    def test_single_sample_still_works(self, multi_sample_vcf_path: Path) -> None:
        """Test that single sample parsing still works."""
        variants = list(
            parse_vcf(
                vcf_path=multi_sample_vcf_path,
                high_bulk="MUT1",
                low_bulk="WT1",
                min_dp=10,
                max_dp=200,
            )
        )

        assert len(variants) == 3

        # Single sample - no pooling
        v1 = variants[0]
        assert v1.dp_high == 50
        # MUT1 at pos 1000: 5 ref, 45 alt -> AF = 45/50 = 0.9
        assert v1.af_high == pytest.approx(0.9, rel=1e-6)
        assert v1.n_samples_high == 1
        assert v1.n_samples_low == 1

    def test_missing_sample_raises_error(self, multi_sample_vcf_path: Path) -> None:
        """Test that missing sample names raise ValueError."""
        with pytest.raises(ValueError, match="not found in VCF"):
            list(
                parse_vcf(
                    vcf_path=multi_sample_vcf_path,
                    high_bulk="MUT1,NONEXISTENT",
                    low_bulk="WT1",
                )
            )

    def test_pooled_depth_calculation(self, multi_sample_vcf_path: Path) -> None:
        """Test that pooled depth is sum of sample depths."""
        variants = list(
            parse_vcf(
                vcf_path=multi_sample_vcf_path,
                high_bulk="MUT1,MUT2",
                low_bulk="WT1,WT2",
                min_dp=10,
                max_dp=500,
            )
        )

        # Each sample has DP=50, so pooled should be 100
        for v in variants:
            assert v.dp_high == 100
            assert v.dp_low == 100
