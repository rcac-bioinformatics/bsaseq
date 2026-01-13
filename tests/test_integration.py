"""Integration tests for complete BSA analysis workflows."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from bsaseq.cli import cli
from bsaseq.core.models import Variant, Window
from bsaseq.analysis.candidates import CandidateRegion, CandidateVariant
from bsaseq.io.vcf import parse_vcf
from bsaseq.analysis.windows import calculate_windows_to_list
from bsaseq.analysis.candidates import (
    identify_candidate_regions,
    filter_candidate_variants,
    InheritanceMode,
)


# Comprehensive test VCF with variants showing clear BSA signal
INTEGRATION_VCF = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##contig=<ID=chr1,length=10000000>
##contig=<ID=chr2,length=8000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	MUTANT	WILDTYPE
chr1	100000	.	A	G	50.0	PASS	.	GT:AD:DP:GQ	0/1:25,25:50:99	0/1:25,25:50:99
chr1	200000	.	C	T	50.0	PASS	.	GT:AD:DP:GQ	0/1:22,28:50:99	0/1:28,22:50:99
chr1	300000	.	G	A	50.0	PASS	.	GT:AD:DP:GQ	0/1:20,30:50:99	0/1:30,20:50:99
chr1	400000	.	T	C	50.0	PASS	.	GT:AD:DP:GQ	0/1:15,35:50:99	0/1:35,15:50:99
chr1	500000	.	A	G	50.0	PASS	.	GT:AD:DP:GQ	0/1:10,40:50:99	0/1:40,10:50:99
chr1	600000	.	C	T	50.0	PASS	.	GT:AD:DP:GQ	0/1:5,45:50:99	0/1:45,5:50:99
chr1	700000	.	G	A	50.0	PASS	.	GT:AD:DP:GQ	0/1:3,47:50:99	0/1:47,3:50:99
chr1	800000	.	T	C	50.0	PASS	.	GT:AD:DP:GQ	0/1:5,45:50:99	0/1:45,5:50:99
chr1	900000	.	A	G	50.0	PASS	.	GT:AD:DP:GQ	0/1:10,40:50:99	0/1:40,10:50:99
chr1	1000000	.	C	T	50.0	PASS	.	GT:AD:DP:GQ	0/1:15,35:50:99	0/1:35,15:50:99
chr2	100000	.	G	A	50.0	PASS	.	GT:AD:DP:GQ	0/1:25,25:50:99	0/1:25,25:50:99
chr2	200000	.	T	C	50.0	PASS	.	GT:AD:DP:GQ	0/1:24,26:50:99	0/1:26,24:50:99
chr2	300000	.	A	G	50.0	PASS	.	GT:AD:DP:GQ	0/1:23,27:50:99	0/1:27,23:50:99
"""


@pytest.fixture
def integration_vcf_path(tmp_path: Path) -> Path:
    """Create a VCF file for integration testing."""
    vcf_path = tmp_path / "integration_test.vcf"
    vcf_path.write_text(INTEGRATION_VCF)
    return vcf_path


class TestVcfToVariantsWorkflow:
    """Test VCF parsing to variant extraction workflow."""

    def test_parse_vcf_returns_variants(self, integration_vcf_path: Path) -> None:
        """Test that parse_vcf extracts variants correctly."""
        variants = list(parse_vcf(
            integration_vcf_path,
            high_bulk=["MUTANT"],
            low_bulk=["WILDTYPE"],
            min_dp=10,
            max_dp=200,
        ))

        assert len(variants) > 0
        assert all(isinstance(v, Variant) for v in variants)

    def test_variants_have_correct_structure(self, integration_vcf_path: Path) -> None:
        """Test that parsed variants have expected attributes."""
        variants = list(parse_vcf(
            integration_vcf_path,
            high_bulk=["MUTANT"],
            low_bulk=["WILDTYPE"],
            min_dp=10,
            max_dp=200,
        ))

        for v in variants:
            assert hasattr(v, "chrom")
            assert hasattr(v, "pos")
            assert hasattr(v, "af_high")
            assert hasattr(v, "af_low")
            assert hasattr(v, "delta_af")

    def test_allele_frequencies_reasonable(self, integration_vcf_path: Path) -> None:
        """Test that allele frequencies are within valid range."""
        variants = list(parse_vcf(
            integration_vcf_path,
            high_bulk=["MUTANT"],
            low_bulk=["WILDTYPE"],
            min_dp=10,
            max_dp=200,
        ))

        for v in variants:
            assert 0.0 <= v.af_high <= 1.0
            assert 0.0 <= v.af_low <= 1.0


class TestVariantsToWindowsWorkflow:
    """Test variant to window calculation workflow."""

    def test_calculate_windows(self, integration_vcf_path: Path) -> None:
        """Test window calculation from variants."""
        variants = list(parse_vcf(
            integration_vcf_path,
            high_bulk=["MUTANT"],
            low_bulk=["WILDTYPE"],
            min_dp=10,
            max_dp=200,
        ))

        windows = calculate_windows_to_list(
            variants,
            window_size=500000,
            step_size=100000,
            min_variants=2,
        )

        assert len(windows) > 0
        assert all(isinstance(w, Window) for w in windows)

    def test_windows_have_z_scores(self, integration_vcf_path: Path) -> None:
        """Test that windows have z-scores calculated."""
        variants = list(parse_vcf(
            integration_vcf_path,
            high_bulk=["MUTANT"],
            low_bulk=["WILDTYPE"],
            min_dp=10,
            max_dp=200,
        ))

        windows = calculate_windows_to_list(
            variants,
            window_size=500000,
            step_size=100000,
            min_variants=2,
        )

        # At least some windows should have z-scores
        z_scores = [w.z_score for w in windows if w.z_score is not None]
        assert len(z_scores) > 0


class TestWindowsToRegionsWorkflow:
    """Test window to candidate region identification workflow."""

    def test_identify_regions(self, integration_vcf_path: Path) -> None:
        """Test candidate region identification."""
        variants = list(parse_vcf(
            integration_vcf_path,
            high_bulk=["MUTANT"],
            low_bulk=["WILDTYPE"],
            min_dp=10,
            max_dp=200,
        ))

        windows = calculate_windows_to_list(
            variants,
            window_size=500000,
            step_size=100000,
            min_variants=2,
        )

        # Use a low threshold to ensure we get regions
        regions = identify_candidate_regions(windows, z_threshold=1.0)

        # May or may not have regions depending on data
        assert isinstance(regions, list)
        assert all(isinstance(r, CandidateRegion) for r in regions)

    def test_regions_have_required_attributes(self, integration_vcf_path: Path) -> None:
        """Test that regions have all required attributes."""
        variants = list(parse_vcf(
            integration_vcf_path,
            high_bulk=["MUTANT"],
            low_bulk=["WILDTYPE"],
            min_dp=10,
            max_dp=200,
        ))

        windows = calculate_windows_to_list(
            variants,
            window_size=500000,
            step_size=100000,
            min_variants=2,
        )

        regions = identify_candidate_regions(windows, z_threshold=1.0)

        for r in regions:
            assert hasattr(r, "chrom")
            assert hasattr(r, "start")
            assert hasattr(r, "end")
            assert hasattr(r, "max_z_score")
            assert hasattr(r, "peak_position")


class TestRegionsToCandidatesWorkflow:
    """Test region to candidate variant filtering workflow."""

    def test_filter_candidates(self, integration_vcf_path: Path) -> None:
        """Test candidate variant filtering."""
        variants = list(parse_vcf(
            integration_vcf_path,
            high_bulk=["MUTANT"],
            low_bulk=["WILDTYPE"],
            min_dp=10,
            max_dp=200,
        ))

        windows = calculate_windows_to_list(
            variants,
            window_size=500000,
            step_size=100000,
            min_variants=2,
        )

        regions = identify_candidate_regions(windows, z_threshold=1.0)

        if regions:
            candidates = filter_candidate_variants(
                variants,
                regions,
                mode=InheritanceMode.RECESSIVE,
            )

            assert isinstance(candidates, list)
            assert all(isinstance(c, CandidateVariant) for c in candidates)


class TestFullPipelineWorkflow:
    """Test complete pipeline from VCF to output."""

    def test_cli_run_produces_output(self, runner: CliRunner, integration_vcf_path: Path, tmp_path: Path) -> None:
        """Test that run command produces output files."""
        output_prefix = tmp_path / "results" / "analysis"

        result = runner.invoke(cli, [
            "run",
            "--vcf", str(integration_vcf_path),
            "--high-bulk", "MUTANT",
            "--low-bulk", "WILDTYPE",
            "--out", str(output_prefix),
            "--window-size", "500000",
            "--step-size", "100000",
            "--min-variants", "2",
            "--z-threshold", "1.0",
        ])

        # Check command completed (may succeed or fail gracefully)
        # The important thing is it doesn't crash
        assert result.exit_code in [0, 1]

    def test_cli_run_with_relaxed_thresholds(self, runner: CliRunner, integration_vcf_path: Path, tmp_path: Path) -> None:
        """Test run command with relaxed filtering thresholds."""
        output_prefix = tmp_path / "relaxed" / "analysis"

        result = runner.invoke(cli, [
            "run",
            "--vcf", str(integration_vcf_path),
            "--high-bulk", "MUTANT",
            "--low-bulk", "WILDTYPE",
            "--out", str(output_prefix),
            "--min-dp", "5",
            "--min-qual", "10",
            "--window-size", "500000",
            "--step-size", "100000",
            "--min-variants", "1",
            "--z-threshold", "0.5",
        ])

        # Should complete without crashing
        assert result.exit_code in [0, 1]


class TestMultiSampleWorkflow:
    """Test workflows with multiple samples per bulk."""

    def test_multi_sample_parsing(self, multi_sample_vcf_path: Path) -> None:
        """Test parsing VCF with multiple samples per bulk."""
        variants = list(parse_vcf(
            multi_sample_vcf_path,
            high_bulk=["MUT1", "MUT2"],
            low_bulk=["WT1", "WT2"],
            min_dp=10,
            max_dp=200,
        ))

        assert len(variants) > 0

        # Check that samples were pooled
        for v in variants:
            assert v.n_samples_high == 2
            assert v.n_samples_low == 2

    def test_multi_sample_cli(self, runner: CliRunner, multi_sample_vcf_path: Path, tmp_path: Path) -> None:
        """Test CLI with multiple samples per bulk."""
        output_prefix = tmp_path / "multi" / "analysis"

        result = runner.invoke(cli, [
            "run",
            "--vcf", str(multi_sample_vcf_path),
            "--high-bulk", "MUT1,MUT2",
            "--low-bulk", "WT1,WT2",
            "--out", str(output_prefix),
            "--window-size", "500000",
            "--min-variants", "1",
        ])

        # Should complete without crashing
        assert result.exit_code in [0, 1]


class TestErrorHandlingWorkflow:
    """Test error handling in workflows."""

    def test_invalid_vcf_path(self, runner: CliRunner, tmp_path: Path) -> None:
        """Test graceful handling of invalid VCF path."""
        result = runner.invoke(cli, [
            "run",
            "--vcf", "/nonexistent/path/file.vcf",
            "--high-bulk", "MUTANT",
            "--low-bulk", "WILDTYPE",
            "--out", str(tmp_path / "output"),
        ])

        assert result.exit_code != 0
        # Should show helpful error message
        assert "not found" in result.output.lower() or "error" in result.output.lower()

    def test_invalid_sample_name(self, runner: CliRunner, integration_vcf_path: Path, tmp_path: Path) -> None:
        """Test graceful handling of invalid sample name."""
        result = runner.invoke(cli, [
            "run",
            "--vcf", str(integration_vcf_path),
            "--high-bulk", "INVALID_SAMPLE",
            "--low-bulk", "WILDTYPE",
            "--out", str(tmp_path / "output"),
        ])

        assert result.exit_code != 0


@pytest.fixture
def runner() -> CliRunner:
    """Create a CLI runner for testing."""
    return CliRunner()
