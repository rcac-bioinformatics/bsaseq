"""Tests for plotting modules."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend for testing

import matplotlib.pyplot as plt
import pytest

from bsaseq.analysis.candidates import CandidateRegion, CandidateVariant
from bsaseq.core.models import Variant, Window
from bsaseq.plotting.diagnostics import plot_af_distribution, plot_depth_distribution
from bsaseq.plotting.genome import plot_all_regions, plot_genome_wide, plot_region
from bsaseq.plotting.style import set_publication_style
from bsaseq.utils.sorting import (
    natural_sort_key,
    simplify_chromosome_label,
    sort_chromosomes,
)


# ============ Fixtures ============


@pytest.fixture
def sample_variants() -> list[Variant]:
    """Generate sample variants for testing."""
    variants = []
    # Create variants across 2 chromosomes
    for chrom in ["chr1", "chr2"]:
        for pos in range(100_000, 2_000_000, 50_000):
            af_high = 0.3 + (pos / 5_000_000)  # Varies 0.3 to 0.7
            af_low = 0.2 + (pos / 10_000_000)  # Varies 0.2 to 0.4
            variants.append(
                Variant(
                    chrom=chrom,
                    pos=pos,
                    ref="A",
                    alt="G",
                    af_high=min(1.0, af_high),
                    af_low=min(1.0, af_low),
                    dp_high=50,
                    dp_low=50,
                )
            )
    return variants


@pytest.fixture
def sample_windows() -> list[Window]:
    """Generate sample windows for testing."""
    windows = []
    # Create windows across 2 chromosomes
    for chrom in ["chr1", "chr2"]:
        for start in range(1, 3_000_000, 250_000):
            end = start + 1_000_000 - 1
            # Create a peak in chr1 around position 1.5M
            base_delta = 0.2
            if chrom == "chr1" and 1_000_000 <= start <= 2_000_000:
                base_delta = 0.6
            z_score = base_delta * 10  # Scale to z-score
            windows.append(
                Window(
                    chrom=chrom,
                    start=start,
                    end=end,
                    n_variants=10,
                    mean_delta_af=base_delta,
                    median_delta_af=base_delta,
                    tricube_delta_af=base_delta,
                    mean_af_high=0.6,
                    mean_af_low=0.4,
                    g_statistic=10.0,
                    p_value=0.01,
                    z_score=z_score,
                )
            )
    return windows


@pytest.fixture
def sample_region() -> CandidateRegion:
    """Generate sample candidate region."""
    return CandidateRegion(
        chrom="chr1",
        start=1_000_000,
        end=2_000_000,
        n_windows=4,
        n_variants=20,
        max_z_score=6.0,
        mean_z_score=5.0,
        max_delta_af=0.7,
        mean_delta_af=0.6,
        peak_position=1_500_000,
    )


@pytest.fixture
def sample_regions() -> list[CandidateRegion]:
    """Generate multiple sample candidate regions."""
    return [
        CandidateRegion(
            chrom="chr1",
            start=1_000_000,
            end=2_000_000,
            n_windows=4,
            n_variants=20,
            max_z_score=6.0,
            mean_z_score=5.0,
            max_delta_af=0.7,
            mean_delta_af=0.6,
            peak_position=1_500_000,
        ),
        CandidateRegion(
            chrom="chr2",
            start=500_000,
            end=1_500_000,
            n_windows=3,
            n_variants=15,
            max_z_score=4.0,
            mean_z_score=3.5,
            max_delta_af=0.5,
            mean_delta_af=0.4,
            peak_position=1_000_000,
        ),
    ]


@pytest.fixture
def sample_candidate_variants() -> list[CandidateVariant]:
    """Generate sample candidate variants."""
    return [
        CandidateVariant(
            chrom="chr1",
            pos=1_500_000,
            ref="A",
            alt="G",
            af_high=0.95,
            af_low=0.05,
            delta_af=0.9,
            dp_high=50,
            dp_low=50,
            region_rank=1,
            distance_to_peak=0,
        ),
        CandidateVariant(
            chrom="chr1",
            pos=1_600_000,
            ref="C",
            alt="T",
            af_high=0.92,
            af_low=0.08,
            delta_af=0.84,
            dp_high=55,
            dp_low=45,
            region_rank=1,
            distance_to_peak=100_000,
        ),
    ]


# ============ Sorting Tests ============


class TestNaturalSorting:
    """Tests for chromosome sorting utilities."""

    def test_natural_sort_key_basic(self) -> None:
        """Test basic natural sort key generation."""
        assert natural_sort_key("chr1") < natural_sort_key("chr2")
        assert natural_sort_key("chr2") < natural_sort_key("chr10")
        assert natural_sort_key("chr10") < natural_sort_key("chrX")

    def test_natural_sort_chromosomes(self) -> None:
        """Test chromosome sorting."""
        chroms = ["Chr10", "Chr2", "Chr1", "ChrX", "Chr3"]
        expected = ["Chr1", "Chr2", "Chr3", "Chr10", "ChrX"]
        assert sort_chromosomes(chroms) == expected

    def test_natural_sort_lowercase(self) -> None:
        """Test case-insensitive sorting."""
        chroms = ["chr10", "chr2", "chr1"]
        expected = ["chr1", "chr2", "chr10"]
        assert sort_chromosomes(chroms) == expected

    def test_natural_sort_numeric_only(self) -> None:
        """Test sorting numeric-only chromosome names."""
        chroms = ["10", "2", "1", "X"]
        expected = ["1", "2", "10", "X"]
        assert sort_chromosomes(chroms) == expected

    def test_simplify_chromosome_label(self) -> None:
        """Test chromosome label simplification."""
        assert simplify_chromosome_label("Chr01") == "1"
        assert simplify_chromosome_label("chromosome10") == "10"
        assert simplify_chromosome_label("chrX") == "X"
        assert simplify_chromosome_label("X") == "X"


# ============ Style Tests ============


class TestPlottingStyle:
    """Tests for plotting style configuration."""

    def test_set_publication_style(self) -> None:
        """Test that publication style can be set without error."""
        set_publication_style()
        # Check that some key parameters were set
        assert plt.rcParams["axes.spines.top"] is False
        assert plt.rcParams["axes.spines.right"] is False


# ============ Genome-wide Plot Tests ============


class TestGenomeWidePlot:
    """Tests for genome-wide Manhattan plot."""

    def test_genome_wide_plot_creates_file(
        self, tmp_path: Path, sample_windows: list[Window]
    ) -> None:
        """Test that genome-wide plot is created."""
        out = tmp_path / "test_genome"
        fig = plot_genome_wide(sample_windows, output_path=out)

        assert (tmp_path / "test_genome.png").exists()
        assert (tmp_path / "test_genome.pdf").exists()
        plt.close(fig)

    def test_genome_wide_plot_with_regions(
        self,
        tmp_path: Path,
        sample_windows: list[Window],
        sample_regions: list[CandidateRegion],
    ) -> None:
        """Test genome plot with candidate regions highlighted."""
        out = tmp_path / "test_genome_regions"
        fig = plot_genome_wide(
            sample_windows, regions=sample_regions, output_path=out
        )

        assert (tmp_path / "test_genome_regions.png").exists()
        plt.close(fig)

    def test_genome_wide_plot_without_regions(
        self, tmp_path: Path, sample_windows: list[Window]
    ) -> None:
        """Test genome plot works without candidate regions."""
        fig = plot_genome_wide(sample_windows, regions=None, output_path=tmp_path / "test")
        assert fig is not None
        plt.close(fig)

    def test_genome_wide_plot_custom_threshold(
        self, tmp_path: Path, sample_windows: list[Window]
    ) -> None:
        """Test genome plot with custom z-threshold."""
        fig = plot_genome_wide(
            sample_windows, z_threshold=5.0, output_path=tmp_path / "test"
        )
        assert fig is not None
        plt.close(fig)

    def test_genome_wide_plot_empty_windows(self, tmp_path: Path) -> None:
        """Test genome plot with empty windows list."""
        fig = plot_genome_wide([], output_path=tmp_path / "test")
        assert fig is not None
        plt.close(fig)

    def test_genome_wide_plot_custom_figsize(
        self, tmp_path: Path, sample_windows: list[Window]
    ) -> None:
        """Test genome plot with custom figure size."""
        fig = plot_genome_wide(
            sample_windows, figsize=(10, 6), output_path=tmp_path / "test"
        )
        assert fig is not None
        plt.close(fig)


# ============ Region Plot Tests ============


class TestRegionPlot:
    """Tests for regional zoom plots."""

    def test_region_plot_creates_file(
        self,
        tmp_path: Path,
        sample_variants: list[Variant],
        sample_region: CandidateRegion,
    ) -> None:
        """Test that region plot is created."""
        out = tmp_path / "test_region"
        fig = plot_region(sample_variants, sample_region, output_path=out)

        assert (tmp_path / "test_region.png").exists()
        assert (tmp_path / "test_region.pdf").exists()
        plt.close(fig)

    def test_region_plot_with_windows(
        self,
        tmp_path: Path,
        sample_variants: list[Variant],
        sample_region: CandidateRegion,
        sample_windows: list[Window],
    ) -> None:
        """Test region plot with window overlay."""
        fig = plot_region(
            sample_variants,
            sample_region,
            windows=sample_windows,
            output_path=tmp_path / "test",
        )
        assert fig is not None
        plt.close(fig)

    def test_region_plot_with_highlights(
        self,
        tmp_path: Path,
        sample_variants: list[Variant],
        sample_region: CandidateRegion,
        sample_candidate_variants: list[CandidateVariant],
    ) -> None:
        """Test region plot with highlighted variants."""
        fig = plot_region(
            sample_variants,
            sample_region,
            highlight_variants=sample_candidate_variants,
            output_path=tmp_path / "test",
        )
        assert fig is not None
        plt.close(fig)

    def test_region_plot_empty_region(
        self, tmp_path: Path, sample_variants: list[Variant]
    ) -> None:
        """Test region plot with region containing no variants."""
        # Create region in a different location with no variants
        empty_region = CandidateRegion(
            chrom="chr3",  # No variants on chr3
            start=1_000_000,
            end=2_000_000,
            n_windows=2,
            n_variants=0,
            max_z_score=4.0,
            mean_z_score=3.5,
            max_delta_af=0.5,
            mean_delta_af=0.4,
            peak_position=1_500_000,
        )
        fig = plot_region(sample_variants, empty_region, output_path=tmp_path / "test")
        assert fig is not None
        plt.close(fig)


# ============ Batch Plot Tests ============


class TestBatchPlots:
    """Tests for batch plot generation."""

    def test_plot_all_regions(
        self,
        tmp_path: Path,
        sample_variants: list[Variant],
        sample_regions: list[CandidateRegion],
        sample_windows: list[Window],
        sample_candidate_variants: list[CandidateVariant],
    ) -> None:
        """Test batch generation of region plots."""
        paths = plot_all_regions(
            variants=sample_variants,
            regions=sample_regions,
            windows=sample_windows,
            candidate_variants=sample_candidate_variants,
            output_dir=tmp_path,
            max_regions=2,
        )

        # Should generate plots for both regions
        assert len(paths) >= 2  # At least 2 files (png + pdf for each)

    def test_plot_all_regions_max_limit(
        self,
        tmp_path: Path,
        sample_variants: list[Variant],
        sample_regions: list[CandidateRegion],
        sample_windows: list[Window],
        sample_candidate_variants: list[CandidateVariant],
    ) -> None:
        """Test that max_regions limits output."""
        paths = plot_all_regions(
            variants=sample_variants,
            regions=sample_regions,
            windows=sample_windows,
            candidate_variants=sample_candidate_variants,
            output_dir=tmp_path,
            max_regions=1,
        )

        # Should only generate for 1 region (2 files: png + pdf)
        assert len(paths) == 2

    def test_plot_all_regions_empty(self, tmp_path: Path) -> None:
        """Test batch plot with no regions."""
        paths = plot_all_regions(
            variants=[],
            regions=[],
            windows=[],
            candidate_variants=[],
            output_dir=tmp_path,
        )

        assert len(paths) == 0


# ============ Diagnostic Plot Tests ============


class TestDiagnosticPlots:
    """Tests for diagnostic plots."""

    def test_af_distribution_creates_file(
        self, tmp_path: Path, sample_variants: list[Variant]
    ) -> None:
        """Test that AF distribution plot is created."""
        out = tmp_path / "test_af"
        fig = plot_af_distribution(sample_variants, output_path=out)

        assert (tmp_path / "test_af.png").exists()
        plt.close(fig)

    def test_af_distribution_empty_variants(self, tmp_path: Path) -> None:
        """Test AF distribution with empty variants."""
        fig = plot_af_distribution([], output_path=tmp_path / "test")
        assert fig is not None
        plt.close(fig)

    def test_depth_distribution_creates_file(
        self, tmp_path: Path, sample_variants: list[Variant]
    ) -> None:
        """Test that depth distribution plot is created."""
        out = tmp_path / "test_depth"
        fig = plot_depth_distribution(sample_variants, output_path=out)

        assert (tmp_path / "test_depth.png").exists()
        plt.close(fig)

    def test_depth_distribution_empty_variants(self, tmp_path: Path) -> None:
        """Test depth distribution with empty variants."""
        fig = plot_depth_distribution([], output_path=tmp_path / "test")
        assert fig is not None
        plt.close(fig)


# ============ TSV Reader Tests ============


class TestTSVReaders:
    """Tests for TSV file readers."""

    def test_read_variants_tsv(self, tmp_path: Path, sample_variants: list[Variant]) -> None:
        """Test reading variants from TSV."""
        from bsaseq.io.readers import read_variants_tsv
        from bsaseq.io.writers import write_variants_tsv

        # Write variants to file
        tsv_path = tmp_path / "variants.tsv"
        write_variants_tsv(sample_variants, tsv_path)

        # Read them back
        read_variants = read_variants_tsv(tsv_path)

        assert len(read_variants) == len(sample_variants)
        assert read_variants[0].chrom == sample_variants[0].chrom
        assert read_variants[0].pos == sample_variants[0].pos

    def test_read_windows_tsv(self, tmp_path: Path, sample_windows: list[Window]) -> None:
        """Test reading windows from TSV."""
        from bsaseq.io.readers import read_windows_tsv
        from bsaseq.io.writers import write_windows_tsv

        # Write windows to file
        tsv_path = tmp_path / "windows.tsv"
        write_windows_tsv(sample_windows, tsv_path)

        # Read them back
        read_windows = read_windows_tsv(tsv_path)

        assert len(read_windows) == len(sample_windows)
        assert read_windows[0].chrom == sample_windows[0].chrom
        assert read_windows[0].start == sample_windows[0].start

    def test_read_regions_tsv(
        self, tmp_path: Path, sample_regions: list[CandidateRegion]
    ) -> None:
        """Test reading regions from TSV."""
        from bsaseq.io.readers import read_regions_tsv
        from bsaseq.io.writers import write_regions_tsv

        # Write regions to file
        tsv_path = tmp_path / "regions.tsv"
        write_regions_tsv(sample_regions, tsv_path)

        # Read them back
        read_regions = read_regions_tsv(tsv_path)

        assert len(read_regions) == len(sample_regions)
        assert read_regions[0].chrom == sample_regions[0].chrom
        assert read_regions[0].start == sample_regions[0].start

    def test_read_missing_file(self, tmp_path: Path) -> None:
        """Test reading non-existent file raises error."""
        from bsaseq.io.readers import read_variants_tsv

        with pytest.raises(FileNotFoundError):
            read_variants_tsv(tmp_path / "nonexistent.tsv")
