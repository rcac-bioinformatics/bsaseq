"""Tests for candidate region detection and variant filtering."""

from __future__ import annotations

import pytest

from bsaseq.analysis.candidates import (
    CandidateRegion,
    CandidateVariant,
    InheritanceMode,
    count_variants_in_regions,
    filter_candidate_variants,
    get_variant_filter_thresholds,
    identify_candidate_regions,
    identify_candidate_regions_percentile,
)
from bsaseq.core.models import Variant, Window


def make_variant(
    chrom: str = "chr1",
    pos: int = 1000,
    ref: str = "A",
    alt: str = "G",
    af_high: float = 0.5,
    af_low: float = 0.5,
    dp_high: int = 50,
    dp_low: int = 50,
    n_samples_high: int = 1,
    n_samples_low: int = 1,
) -> Variant:
    """Helper to create test variants."""
    return Variant(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        af_high=af_high,
        af_low=af_low,
        dp_high=dp_high,
        dp_low=dp_low,
        n_samples_high=n_samples_high,
        n_samples_low=n_samples_low,
    )


def make_window(
    chrom: str = "chr1",
    start: int = 1,
    end: int = 1000000,
    n_variants: int = 10,
    mean_delta_af: float = 0.1,
    z_score: float | None = 1.0,
) -> Window:
    """Helper to create test windows."""
    return Window(
        chrom=chrom,
        start=start,
        end=end,
        n_variants=n_variants,
        mean_delta_af=mean_delta_af,
        median_delta_af=mean_delta_af,
        tricube_delta_af=mean_delta_af,
        mean_af_high=0.6,
        mean_af_low=0.5,
        g_statistic=10.0,
        p_value=0.001,
        z_score=z_score,
    )


class TestInheritanceMode:
    """Tests for InheritanceMode enum."""

    def test_recessive_mode(self) -> None:
        """Test recessive inheritance mode."""
        mode = InheritanceMode("recessive")
        assert mode == InheritanceMode.RECESSIVE
        assert mode.value == "recessive"

    def test_dominant_mode(self) -> None:
        """Test dominant inheritance mode."""
        mode = InheritanceMode("dominant")
        assert mode == InheritanceMode.DOMINANT
        assert mode.value == "dominant"


class TestVariantFilterThresholds:
    """Tests for get_variant_filter_thresholds."""

    def test_recessive_thresholds(self) -> None:
        """Test thresholds for recessive mode."""
        thresholds = get_variant_filter_thresholds(InheritanceMode.RECESSIVE)
        assert thresholds["min_delta_af"] == 0.8
        assert thresholds["min_af_high"] == 0.9
        assert thresholds["max_af_low"] == 0.1

    def test_dominant_thresholds(self) -> None:
        """Test thresholds for dominant mode."""
        thresholds = get_variant_filter_thresholds(InheritanceMode.DOMINANT)
        assert thresholds["min_delta_af"] == 0.3
        assert thresholds["min_af_high"] == 0.4
        assert thresholds["max_af_low"] == 0.1


class TestIdentifyCandidateRegions:
    """Tests for identify_candidate_regions function."""

    def test_region_detection_single_peak(self) -> None:
        """Test detection of a single candidate region."""
        # Create windows with one peak region (3 consecutive significant windows)
        windows = [
            make_window(start=1, end=1000000, z_score=1.0),
            make_window(start=250001, end=1250000, z_score=2.0),
            make_window(start=500001, end=1500000, z_score=4.0),  # Significant
            make_window(start=750001, end=1750000, z_score=5.0),  # Significant
            make_window(start=1000001, end=2000000, z_score=4.5),  # Significant
            make_window(start=1250001, end=2250000, z_score=2.0),
            make_window(start=1500001, end=2500000, z_score=1.0),
        ]

        regions = identify_candidate_regions(
            windows=windows,
            z_threshold=3.0,
            min_consecutive=2,
            merge_distance=500_000,
        )

        assert len(regions) == 1
        region = regions[0]
        assert region.chrom == "chr1"
        assert region.n_windows == 3
        assert region.max_z_score == 5.0
        assert region.start == 500001
        assert region.end == 2000000

    def test_region_detection_multiple_peaks(self) -> None:
        """Test detection of multiple separate candidate regions."""
        # Create windows with two separate peak regions
        windows = [
            # First peak
            make_window(chrom="chr1", start=1, end=1000000, z_score=4.0),
            make_window(chrom="chr1", start=250001, end=1250000, z_score=5.0),
            make_window(chrom="chr1", start=500001, end=1500000, z_score=4.0),
            # Gap (non-significant)
            make_window(chrom="chr1", start=5000001, end=6000000, z_score=1.0),
            make_window(chrom="chr1", start=5250001, end=6250000, z_score=1.0),
            # Second peak
            make_window(chrom="chr1", start=10000001, end=11000000, z_score=3.5),
            make_window(chrom="chr1", start=10250001, end=11250000, z_score=4.0),
        ]

        regions = identify_candidate_regions(
            windows=windows,
            z_threshold=3.0,
            min_consecutive=2,
            merge_distance=500_000,
        )

        assert len(regions) == 2
        # Should be sorted by max_z_score descending
        assert regions[0].max_z_score == 5.0
        assert regions[1].max_z_score == 4.0

    def test_region_merging(self) -> None:
        """Test merging of nearby regions on same chromosome."""
        # Create two close regions that should be merged
        windows = [
            make_window(start=1, end=1000000, z_score=4.0),
            make_window(start=250001, end=1250000, z_score=4.5),
            # Small gap (less than merge_distance)
            make_window(start=1300001, end=2300000, z_score=2.0),
            # Second significant region
            make_window(start=1500001, end=2500000, z_score=4.0),
            make_window(start=1750001, end=2750000, z_score=3.5),
        ]

        regions = identify_candidate_regions(
            windows=windows,
            z_threshold=3.0,
            min_consecutive=2,
            merge_distance=500_000,
        )

        # Should merge into one region
        assert len(regions) == 1
        assert regions[0].start == 1
        assert regions[0].end == 2750000
        assert regions[0].max_z_score == 4.5

    def test_region_no_merge_different_chrom(self) -> None:
        """Test that regions on different chromosomes are not merged."""
        windows = [
            # Region on chr1
            make_window(chrom="chr1", start=1, end=1000000, z_score=4.0),
            make_window(chrom="chr1", start=250001, end=1250000, z_score=4.5),
            # Region on chr2
            make_window(chrom="chr2", start=1, end=1000000, z_score=3.5),
            make_window(chrom="chr2", start=250001, end=1250000, z_score=4.0),
        ]

        regions = identify_candidate_regions(
            windows=windows,
            z_threshold=3.0,
            min_consecutive=2,
            merge_distance=500_000,
        )

        assert len(regions) == 2
        chroms = {r.chrom for r in regions}
        assert chroms == {"chr1", "chr2"}

    def test_no_regions_below_threshold(self) -> None:
        """Test that no regions are returned when all windows below threshold."""
        windows = [
            make_window(start=1, end=1000000, z_score=1.0),
            make_window(start=250001, end=1250000, z_score=2.0),
            make_window(start=500001, end=1500000, z_score=2.5),
        ]

        regions = identify_candidate_regions(
            windows=windows, z_threshold=3.0, min_consecutive=2
        )

        assert len(regions) == 0

    def test_single_significant_window_not_enough(self) -> None:
        """Test that single significant window doesn't form region with min_consecutive=2."""
        windows = [
            make_window(start=1, end=1000000, z_score=1.0),
            make_window(start=250001, end=1250000, z_score=5.0),  # Single peak
            make_window(start=500001, end=1500000, z_score=1.0),
        ]

        regions = identify_candidate_regions(
            windows=windows, z_threshold=3.0, min_consecutive=2
        )

        assert len(regions) == 0

    def test_region_peak_position(self) -> None:
        """Test that peak_position is at midpoint of highest z-score window."""
        windows = [
            make_window(start=1, end=1000000, z_score=4.0),
            make_window(start=250001, end=1250000, z_score=6.0),  # Max z-score
            make_window(start=500001, end=1500000, z_score=4.0),
        ]

        regions = identify_candidate_regions(
            windows=windows, z_threshold=3.0, min_consecutive=2
        )

        assert len(regions) == 1
        # Peak should be at midpoint of second window
        expected_peak = (250001 + 1250000) // 2
        assert regions[0].peak_position == expected_peak

    def test_none_z_scores_ignored(self) -> None:
        """Test that windows with None z-scores are ignored."""
        windows = [
            make_window(start=1, end=1000000, z_score=None),
            make_window(start=250001, end=1250000, z_score=4.0),
            make_window(start=500001, end=1500000, z_score=4.0),
        ]

        regions = identify_candidate_regions(
            windows=windows, z_threshold=3.0, min_consecutive=2
        )

        assert len(regions) == 1


class TestIdentifyCandidateRegionsPercentile:
    """Tests for identify_candidate_regions_percentile function."""

    def test_percentile_detection(self) -> None:
        """Test detection using percentile threshold."""
        # Create 100 windows, top 5 should be significant (95th percentile)
        windows = []
        for i in range(100):
            z = float(i) / 10  # Z-scores from 0.0 to 9.9
            delta_af = float(i) / 100  # delta_af from 0.0 to 0.99
            windows.append(
                make_window(
                    start=i * 100000 + 1,
                    end=(i + 1) * 100000,
                    z_score=z,
                    mean_delta_af=delta_af,
                )
            )

        regions = identify_candidate_regions_percentile(
            windows=windows,
            percentile=95.0,
            min_consecutive=2,
            merge_distance=500_000,
        )

        # Should find regions in the high delta_af area
        assert len(regions) >= 1
        # All regions should have high delta_af
        assert all(r.max_delta_af >= 0.9 for r in regions)


class TestCountVariantsInRegions:
    """Tests for count_variants_in_regions function."""

    def test_count_variants(self) -> None:
        """Test counting variants within regions."""
        variants = [
            make_variant(chrom="chr1", pos=500000),
            make_variant(chrom="chr1", pos=750000),
            make_variant(chrom="chr1", pos=1500000),  # Outside
            make_variant(chrom="chr2", pos=100000),
        ]

        regions = [
            CandidateRegion(
                chrom="chr1",
                start=1,
                end=1000000,
                n_windows=3,
                n_variants=0,
                max_z_score=5.0,
                mean_z_score=4.0,
                max_delta_af=0.8,
                mean_delta_af=0.6,
                peak_position=500000,
            ),
            CandidateRegion(
                chrom="chr2",
                start=1,
                end=500000,
                n_windows=2,
                n_variants=0,
                max_z_score=4.0,
                mean_z_score=3.5,
                max_delta_af=0.7,
                mean_delta_af=0.5,
                peak_position=250000,
            ),
        ]

        counts = count_variants_in_regions(variants, regions)

        # Returns dict mapping region index to count
        assert counts[0] == 2  # chr1 has 2 variants in region
        assert counts[1] == 1  # chr2 has 1 variant in region


class TestFilterCandidateVariants:
    """Tests for filter_candidate_variants function."""

    def test_candidate_variant_filtering_recessive(self) -> None:
        """Test filtering for recessive inheritance pattern."""
        variants = [
            # Good candidate - high AF in mutant, low in WT
            make_variant(chrom="chr1", pos=500000, af_high=0.95, af_low=0.05),
            # Not a candidate - AF too low in high bulk
            make_variant(chrom="chr1", pos=600000, af_high=0.7, af_low=0.05),
            # Not a candidate - AF too high in low bulk
            make_variant(chrom="chr1", pos=700000, af_high=0.95, af_low=0.3),
            # Not in region
            make_variant(chrom="chr1", pos=1500000, af_high=0.95, af_low=0.05),
        ]

        regions = [
            CandidateRegion(
                chrom="chr1",
                start=1,
                end=1000000,
                n_windows=3,
                n_variants=3,
                max_z_score=5.0,
                mean_z_score=4.0,
                max_delta_af=0.8,
                mean_delta_af=0.6,
                peak_position=500000,
            ),
        ]

        candidates = filter_candidate_variants(
            variants=variants,
            regions=regions,
            min_delta_af=0.8,
            min_af_high=0.9,
            max_af_low=0.1,
        )

        assert len(candidates) == 1
        assert candidates[0].pos == 500000
        assert candidates[0].region_rank == 1

    def test_candidate_variant_filtering_dominant(self) -> None:
        """Test filtering for dominant inheritance pattern."""
        variants = [
            # Good candidate for dominant - moderate AF difference
            make_variant(chrom="chr1", pos=500000, af_high=0.55, af_low=0.1),
            # Not a candidate - AF too low in high bulk
            make_variant(chrom="chr1", pos=600000, af_high=0.3, af_low=0.1),
            # Not a candidate - AF too high in low bulk
            make_variant(chrom="chr1", pos=700000, af_high=0.55, af_low=0.3),
        ]

        regions = [
            CandidateRegion(
                chrom="chr1",
                start=1,
                end=1000000,
                n_windows=3,
                n_variants=3,
                max_z_score=5.0,
                mean_z_score=4.0,
                max_delta_af=0.5,
                mean_delta_af=0.4,
                peak_position=500000,
            ),
        ]

        thresholds = get_variant_filter_thresholds(InheritanceMode.DOMINANT)
        candidates = filter_candidate_variants(
            variants=variants,
            regions=regions,
            min_delta_af=thresholds["min_delta_af"],
            min_af_high=thresholds["min_af_high"],
            max_af_low=thresholds["max_af_low"],
        )

        assert len(candidates) == 1
        assert candidates[0].pos == 500000

    def test_candidate_variant_sorting(self) -> None:
        """Test that candidates are sorted by region rank then distance_to_peak."""
        variants = [
            # Closest to peak in region 1
            make_variant(chrom="chr1", pos=500000, af_high=0.95, af_low=0.05),
            # Further from peak in region 1
            make_variant(chrom="chr1", pos=700000, af_high=0.98, af_low=0.02),
            # In region 2
            make_variant(chrom="chr2", pos=250000, af_high=0.95, af_low=0.05),
        ]

        regions = [
            CandidateRegion(
                chrom="chr1",
                start=1,
                end=1000000,
                n_windows=3,
                n_variants=2,
                max_z_score=5.0,
                mean_z_score=4.0,
                max_delta_af=0.9,
                mean_delta_af=0.8,
                peak_position=500000,  # Peak at 500000
            ),
            CandidateRegion(
                chrom="chr2",
                start=1,
                end=500000,
                n_windows=2,
                n_variants=1,
                max_z_score=4.0,
                mean_z_score=3.5,
                max_delta_af=0.8,
                mean_delta_af=0.7,
                peak_position=250000,
            ),
        ]

        candidates = filter_candidate_variants(
            variants=variants,
            regions=regions,
            min_delta_af=0.8,
            min_af_high=0.9,
            max_af_low=0.1,
        )

        assert len(candidates) == 3
        # First should be from region 1, closest to peak (distance=0)
        assert candidates[0].region_rank == 1
        assert candidates[0].pos == 500000
        assert candidates[0].distance_to_peak == 0
        # Second from region 1, further from peak (distance=200000)
        assert candidates[1].region_rank == 1
        assert candidates[1].pos == 700000
        assert candidates[1].distance_to_peak == 200000
        # Third from region 2 (distance=0)
        assert candidates[2].region_rank == 2

    def test_distance_to_peak(self) -> None:
        """Test that distance_to_peak is calculated correctly."""
        variants = [
            make_variant(chrom="chr1", pos=500000, af_high=0.95, af_low=0.05),
            make_variant(chrom="chr1", pos=700000, af_high=0.95, af_low=0.05),
        ]

        regions = [
            CandidateRegion(
                chrom="chr1",
                start=1,
                end=1000000,
                n_windows=3,
                n_variants=2,
                max_z_score=5.0,
                mean_z_score=4.0,
                max_delta_af=0.9,
                mean_delta_af=0.8,
                peak_position=600000,
            ),
        ]

        candidates = filter_candidate_variants(
            variants=variants,
            regions=regions,
            min_delta_af=0.8,
            min_af_high=0.9,
            max_af_low=0.1,
        )

        # Check distances
        distances = {c.pos: c.distance_to_peak for c in candidates}
        assert distances[500000] == 100000  # |500000 - 600000|
        assert distances[700000] == 100000  # |700000 - 600000|

    def test_no_regions_returns_empty(self) -> None:
        """Test that empty regions list returns empty candidates."""
        variants = [
            make_variant(chrom="chr1", pos=500000, af_high=0.95, af_low=0.05),
        ]

        candidates = filter_candidate_variants(
            variants=variants,
            regions=[],
            min_delta_af=0.8,
            min_af_high=0.9,
            max_af_low=0.1,
        )

        assert len(candidates) == 0


class TestCandidateVariantDataclass:
    """Tests for CandidateVariant dataclass properties."""

    def test_candidate_variant_creation(self) -> None:
        """Test creating a CandidateVariant."""
        cv = CandidateVariant(
            chrom="chr1",
            pos=500000,
            ref="A",
            alt="G",
            af_high=0.95,
            af_low=0.05,
            delta_af=0.9,
            dp_high=100,
            dp_low=100,
            region_rank=1,
            distance_to_peak=50000,
        )

        assert cv.chrom == "chr1"
        assert cv.pos == 500000
        assert cv.delta_af == pytest.approx(0.9)
        assert cv.region_rank == 1
        assert cv.distance_to_peak == 50000

    def test_candidate_variant_from_variant(self) -> None:
        """Test creating CandidateVariant from Variant."""
        v = make_variant(chrom="chr1", pos=500000, af_high=0.95, af_low=0.05)
        region = CandidateRegion(
            chrom="chr1",
            start=1,
            end=1000000,
            n_windows=3,
            n_variants=2,
            max_z_score=5.0,
            mean_z_score=4.0,
            max_delta_af=0.9,
            mean_delta_af=0.8,
            peak_position=600000,
        )

        cv = CandidateVariant.from_variant(v, region, region_rank=1)

        assert cv.chrom == "chr1"
        assert cv.pos == 500000
        assert cv.delta_af == pytest.approx(0.9)
        assert cv.region_rank == 1
        assert cv.distance_to_peak == 100000  # |500000 - 600000|


class TestCandidateRegionDataclass:
    """Tests for CandidateRegion dataclass properties."""

    def test_region_length(self) -> None:
        """Test region length calculation."""
        region = CandidateRegion(
            chrom="chr1",
            start=1000,
            end=2000,
            n_windows=3,
            n_variants=10,
            max_z_score=5.0,
            mean_z_score=4.0,
            max_delta_af=0.8,
            mean_delta_af=0.6,
            peak_position=1500,
        )

        # Implementation uses end - start (not +1)
        assert region.length == 1000

    def test_region_contains(self) -> None:
        """Test region contains method."""
        region = CandidateRegion(
            chrom="chr1",
            start=1000,
            end=2000,
            n_windows=3,
            n_variants=10,
            max_z_score=5.0,
            mean_z_score=4.0,
            max_delta_af=0.8,
            mean_delta_af=0.6,
            peak_position=1500,
        )

        # Inside region
        assert region.contains("chr1", 1500) is True
        assert region.contains("chr1", 1000) is True
        assert region.contains("chr1", 2000) is True

        # Outside region
        assert region.contains("chr1", 999) is False
        assert region.contains("chr1", 2001) is False

        # Different chromosome
        assert region.contains("chr2", 1500) is False
