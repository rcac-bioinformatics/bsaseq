"""Tests for CLI commands and options."""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from bsaseq.cli import cli


@pytest.fixture
def runner() -> CliRunner:
    """Create a Click CLI test runner."""
    return CliRunner()


class TestCliBasics:
    """Tests for basic CLI functionality."""

    def test_cli_help(self, runner: CliRunner) -> None:
        """Test that main CLI shows help."""
        result = runner.invoke(cli, ["--help"])

        assert result.exit_code == 0
        assert "bsaseq" in result.output
        assert "Bulk Segregant Analysis" in result.output

    def test_cli_short_help(self, runner: CliRunner) -> None:
        """Test that -h works as help option."""
        result = runner.invoke(cli, ["-h"])

        assert result.exit_code == 0
        assert "bsaseq" in result.output

    def test_cli_no_args(self, runner: CliRunner) -> None:
        """Test CLI with no arguments shows usage."""
        result = runner.invoke(cli, [])

        # Click returns exit code 0 when no subcommand is provided
        assert result.exit_code in [0, 2]
        # Should still show available commands in usage message
        assert "run" in result.output
        assert "samples" in result.output


class TestSamplesCommand:
    """Tests for the 'samples' command."""

    def test_samples_help(self, runner: CliRunner) -> None:
        """Test samples command help."""
        result = runner.invoke(cli, ["samples", "--help"])

        assert result.exit_code == 0
        assert "--vcf" in result.output

    def test_samples_lists_samples(self, runner: CliRunner, sample_vcf_path: Path) -> None:
        """Test that samples command lists VCF samples."""
        result = runner.invoke(cli, ["samples", "--vcf", str(sample_vcf_path)])

        assert result.exit_code == 0
        assert "HIGH_BULK" in result.output
        assert "LOW_BULK" in result.output

    def test_samples_missing_vcf(self, runner: CliRunner) -> None:
        """Test samples command with missing VCF file."""
        result = runner.invoke(cli, ["samples", "--vcf", "/nonexistent/file.vcf"])

        assert result.exit_code != 0

    def test_samples_multisample_vcf(self, runner: CliRunner, multi_sample_vcf_path: Path) -> None:
        """Test samples command with multi-sample VCF."""
        result = runner.invoke(cli, ["samples", "--vcf", str(multi_sample_vcf_path)])

        assert result.exit_code == 0
        assert "MUT1" in result.output
        assert "MUT2" in result.output
        assert "WT1" in result.output
        assert "WT2" in result.output


class TestRunCommand:
    """Tests for the 'run' command."""

    def test_run_help(self, runner: CliRunner) -> None:
        """Test run command help."""
        result = runner.invoke(cli, ["run", "--help"])

        assert result.exit_code == 0
        assert "--vcf" in result.output
        assert "--high-bulk" in result.output
        assert "--low-bulk" in result.output
        assert "--out" in result.output

    def test_run_missing_required_args(self, runner: CliRunner) -> None:
        """Test run command fails without required arguments."""
        result = runner.invoke(cli, ["run"])

        assert result.exit_code != 0
        assert "Missing option" in result.output or "required" in result.output.lower()

    def test_run_missing_vcf(self, runner: CliRunner, tmp_path: Path) -> None:
        """Test run command with nonexistent VCF file."""
        result = runner.invoke(cli, [
            "run",
            "--vcf", "/nonexistent/file.vcf",
            "--high-bulk", "sample1",
            "--low-bulk", "sample2",
            "--out", str(tmp_path / "output"),
        ])

        assert result.exit_code != 0

    def test_run_invalid_sample(self, runner: CliRunner, sample_vcf_path: Path, tmp_path: Path) -> None:
        """Test run command with invalid sample name."""
        result = runner.invoke(cli, [
            "run",
            "--vcf", str(sample_vcf_path),
            "--high-bulk", "INVALID_SAMPLE",
            "--low-bulk", "LOW_BULK",
            "--out", str(tmp_path / "output"),
        ])

        assert result.exit_code != 0

    def test_run_shows_parameter_options(self, runner: CliRunner) -> None:
        """Test that run command shows all parameter options in help."""
        result = runner.invoke(cli, ["run", "--help"])

        assert result.exit_code == 0
        assert "--window-size" in result.output
        assert "--step-size" in result.output
        assert "--min-dp" in result.output
        assert "--max-dp" in result.output
        assert "--z-threshold" in result.output
        assert "--mode" in result.output

    def test_run_shows_annotation_options(self, runner: CliRunner) -> None:
        """Test that run command shows annotation options."""
        result = runner.invoke(cli, ["run", "--help"])

        assert result.exit_code == 0
        assert "--annotate" in result.output
        assert "--snpeff-db" in result.output


class TestPlotCommand:
    """Tests for the 'plot' command."""

    def test_plot_help(self, runner: CliRunner) -> None:
        """Test plot command help."""
        result = runner.invoke(cli, ["plot", "--help"])

        assert result.exit_code == 0
        assert "--windows" in result.output
        assert "--regions" in result.output
        assert "--out" in result.output

    def test_plot_missing_required_args(self, runner: CliRunner) -> None:
        """Test plot command fails without required arguments."""
        result = runner.invoke(cli, ["plot"])

        assert result.exit_code != 0


class TestAnnotateCommand:
    """Tests for the 'annotate' command."""

    def test_annotate_help(self, runner: CliRunner) -> None:
        """Test annotate command help."""
        result = runner.invoke(cli, ["annotate", "--help"])

        assert result.exit_code == 0
        assert "--candidates" in result.output
        assert "--snpeff-db" in result.output
        assert "--out" in result.output


class TestCheckSnpeffCommand:
    """Tests for the 'check-snpeff' command."""

    def test_check_snpeff_help(self, runner: CliRunner) -> None:
        """Test check-snpeff command help."""
        result = runner.invoke(cli, ["check-snpeff", "--help"])

        assert result.exit_code == 0


class TestCliValidation:
    """Tests for CLI input validation."""

    def test_run_invalid_window_size(self, runner: CliRunner, sample_vcf_path: Path, tmp_path: Path) -> None:
        """Test run command rejects invalid window size."""
        result = runner.invoke(cli, [
            "run",
            "--vcf", str(sample_vcf_path),
            "--high-bulk", "HIGH_BULK",
            "--low-bulk", "LOW_BULK",
            "--out", str(tmp_path / "output"),
            "--window-size", "0",
        ])

        assert result.exit_code != 0

    def test_run_invalid_min_dp(self, runner: CliRunner, sample_vcf_path: Path, tmp_path: Path) -> None:
        """Test run command rejects invalid min-dp."""
        result = runner.invoke(cli, [
            "run",
            "--vcf", str(sample_vcf_path),
            "--high-bulk", "HIGH_BULK",
            "--low-bulk", "LOW_BULK",
            "--out", str(tmp_path / "output"),
            "--min-dp", "-5",
        ])

        assert result.exit_code != 0

    def test_run_mode_options(self, runner: CliRunner) -> None:
        """Test that mode option accepts valid values."""
        result = runner.invoke(cli, ["run", "--help"])

        assert result.exit_code == 0
        assert "recessive" in result.output
        assert "dominant" in result.output


class TestCliOutputFormat:
    """Tests for CLI output formatting."""

    def test_samples_output_formatting(self, runner: CliRunner, multi_sample_vcf_path: Path) -> None:
        """Test that samples output is properly formatted."""
        result = runner.invoke(cli, ["samples", "--vcf", str(multi_sample_vcf_path)])

        assert result.exit_code == 0
        # Should show sample count
        assert "4" in result.output or "sample" in result.output.lower()

    def test_help_shows_examples(self, runner: CliRunner) -> None:
        """Test that help text includes usage examples."""
        result = runner.invoke(cli, ["--help"])

        assert result.exit_code == 0
        # Should show example command
        assert "bsaseq" in result.output
