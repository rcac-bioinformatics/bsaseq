"""Tests for input validation utilities."""

from __future__ import annotations

from pathlib import Path

import pytest

from bsaseq.utils.validation import (
    ValidationError,
    validate_output_path,
    validate_parameters,
    validate_samples,
    validate_vcf,
)


class TestValidateVcf:
    """Tests for validate_vcf function."""

    def test_valid_vcf(self, sample_vcf_path: Path) -> None:
        """Test validation of a valid VCF file."""
        result = validate_vcf(sample_vcf_path)

        assert "samples" in result
        assert "has_ad" in result
        assert "has_gq" in result
        assert "contigs" in result

        assert result["samples"] == ["HIGH_BULK", "LOW_BULK"]
        assert result["has_ad"] is True
        assert result["has_gq"] is True
        assert "chr1" in result["contigs"]
        assert "chr2" in result["contigs"]

    def test_vcf_file_not_found(self, tmp_path: Path) -> None:
        """Test error when VCF file does not exist."""
        nonexistent = tmp_path / "nonexistent.vcf"

        with pytest.raises(ValidationError, match="VCF file not found"):
            validate_vcf(nonexistent)

    def test_vcf_missing_ad_field(self, vcf_missing_ad_path: Path) -> None:
        """Test error when VCF is missing AD field."""
        with pytest.raises(ValidationError, match="AD.*allelic depth"):
            validate_vcf(vcf_missing_ad_path)

    def test_vcf_missing_gt_field(self, vcf_missing_gt_path: Path) -> None:
        """Test error when VCF is missing GT field."""
        with pytest.raises(ValidationError, match="GT.*genotype"):
            validate_vcf(vcf_missing_gt_path)

    def test_vcf_insufficient_samples(self, vcf_single_sample_path: Path) -> None:
        """Test error when VCF has fewer than 2 samples."""
        with pytest.raises(ValidationError, match="at least 2 samples"):
            validate_vcf(vcf_single_sample_path)

    def test_vcf_multi_sample(self, multi_sample_vcf_path: Path) -> None:
        """Test validation of multi-sample VCF."""
        result = validate_vcf(multi_sample_vcf_path)

        assert len(result["samples"]) == 4
        assert "MUT1" in result["samples"]
        assert "MUT2" in result["samples"]
        assert "WT1" in result["samples"]
        assert "WT2" in result["samples"]

    def test_vcf_returns_contigs(self, sample_vcf_path: Path) -> None:
        """Test that validate_vcf returns contig information."""
        result = validate_vcf(sample_vcf_path)

        assert len(result["contigs"]) == 2
        assert "chr1" in result["contigs"]
        assert "chr2" in result["contigs"]


class TestValidateSamples:
    """Tests for validate_samples function."""

    def test_valid_samples_as_strings(self, sample_vcf_path: Path) -> None:
        """Test validation with sample names as strings."""
        high_list, low_list = validate_samples(
            sample_vcf_path,
            high_bulk="HIGH_BULK",
            low_bulk="LOW_BULK",
        )

        assert high_list == ["HIGH_BULK"]
        assert low_list == ["LOW_BULK"]

    def test_valid_samples_as_lists(self, sample_vcf_path: Path) -> None:
        """Test validation with sample names as lists."""
        high_list, low_list = validate_samples(
            sample_vcf_path,
            high_bulk=["HIGH_BULK"],
            low_bulk=["LOW_BULK"],
        )

        assert high_list == ["HIGH_BULK"]
        assert low_list == ["LOW_BULK"]

    def test_comma_separated_samples(self, multi_sample_vcf_path: Path) -> None:
        """Test validation with comma-separated sample names."""
        high_list, low_list = validate_samples(
            multi_sample_vcf_path,
            high_bulk="MUT1,MUT2",
            low_bulk="WT1,WT2",
        )

        assert high_list == ["MUT1", "MUT2"]
        assert low_list == ["WT1", "WT2"]

    def test_sample_not_found(self, sample_vcf_path: Path) -> None:
        """Test error when sample name not found in VCF."""
        with pytest.raises(ValidationError, match="not found in VCF"):
            validate_samples(
                sample_vcf_path,
                high_bulk="NONEXISTENT",
                low_bulk="LOW_BULK",
            )

    def test_sample_overlap(self, sample_vcf_path: Path) -> None:
        """Test error when sample appears in both bulks."""
        with pytest.raises(ValidationError, match="cannot be in both"):
            validate_samples(
                sample_vcf_path,
                high_bulk="HIGH_BULK",
                low_bulk="HIGH_BULK",
            )

    def test_case_sensitivity(self, sample_vcf_path: Path) -> None:
        """Test that sample names are case-sensitive."""
        with pytest.raises(ValidationError, match="not found"):
            validate_samples(
                sample_vcf_path,
                high_bulk="high_bulk",  # lowercase
                low_bulk="LOW_BULK",
            )

    def test_whitespace_handling(self, multi_sample_vcf_path: Path) -> None:
        """Test that whitespace around sample names is trimmed."""
        high_list, low_list = validate_samples(
            multi_sample_vcf_path,
            high_bulk=" MUT1 , MUT2 ",
            low_bulk=" WT1 , WT2 ",
        )

        assert high_list == ["MUT1", "MUT2"]
        assert low_list == ["WT1", "WT2"]

    def test_partial_overlap(self, multi_sample_vcf_path: Path) -> None:
        """Test error when one sample appears in both bulks."""
        with pytest.raises(ValidationError, match="MUT1"):
            validate_samples(
                multi_sample_vcf_path,
                high_bulk="MUT1,MUT2",
                low_bulk="MUT1,WT1",  # MUT1 in both
            )


class TestValidateParameters:
    """Tests for validate_parameters function."""

    def test_valid_parameters(self) -> None:
        """Test validation of valid parameters."""
        # Should not raise any exception
        validate_parameters(
            window_size=1000000,
            step_size=250000,
            min_dp=10,
            max_dp=200,
            z_threshold=3.0,
        )

    def test_invalid_window_size_zero(self) -> None:
        """Test error when window_size is zero."""
        with pytest.raises(ValidationError, match="window_size must be positive"):
            validate_parameters(
                window_size=0,
                step_size=250000,
                min_dp=10,
                max_dp=200,
                z_threshold=3.0,
            )

    def test_invalid_window_size_negative(self) -> None:
        """Test error when window_size is negative."""
        with pytest.raises(ValidationError, match="window_size must be positive"):
            validate_parameters(
                window_size=-1000,
                step_size=250000,
                min_dp=10,
                max_dp=200,
                z_threshold=3.0,
            )

    def test_invalid_step_size_zero(self) -> None:
        """Test error when step_size is zero."""
        with pytest.raises(ValidationError, match="step_size must be positive"):
            validate_parameters(
                window_size=1000000,
                step_size=0,
                min_dp=10,
                max_dp=200,
                z_threshold=3.0,
            )

    def test_step_size_larger_than_window(self) -> None:
        """Test error when step_size exceeds window_size."""
        with pytest.raises(ValidationError, match="step_size.*cannot be larger"):
            validate_parameters(
                window_size=1000000,
                step_size=2000000,
                min_dp=10,
                max_dp=200,
                z_threshold=3.0,
            )

    def test_invalid_min_dp_zero(self) -> None:
        """Test error when min_dp is zero."""
        with pytest.raises(ValidationError, match="min_dp must be positive"):
            validate_parameters(
                window_size=1000000,
                step_size=250000,
                min_dp=0,
                max_dp=200,
                z_threshold=3.0,
            )

    def test_invalid_max_dp_zero(self) -> None:
        """Test error when max_dp is zero."""
        with pytest.raises(ValidationError, match="max_dp must be positive"):
            validate_parameters(
                window_size=1000000,
                step_size=250000,
                min_dp=10,
                max_dp=0,
                z_threshold=3.0,
            )

    def test_min_dp_greater_than_max_dp(self) -> None:
        """Test error when min_dp >= max_dp."""
        with pytest.raises(ValidationError, match="min_dp.*must be less than max_dp"):
            validate_parameters(
                window_size=1000000,
                step_size=250000,
                min_dp=200,
                max_dp=100,
                z_threshold=3.0,
            )

    def test_min_dp_equal_to_max_dp(self) -> None:
        """Test error when min_dp equals max_dp."""
        with pytest.raises(ValidationError, match="min_dp.*must be less than max_dp"):
            validate_parameters(
                window_size=1000000,
                step_size=250000,
                min_dp=100,
                max_dp=100,
                z_threshold=3.0,
            )

    def test_invalid_z_threshold_zero(self) -> None:
        """Test error when z_threshold is zero."""
        with pytest.raises(ValidationError, match="z_threshold must be positive"):
            validate_parameters(
                window_size=1000000,
                step_size=250000,
                min_dp=10,
                max_dp=200,
                z_threshold=0.0,
            )

    def test_invalid_z_threshold_negative(self) -> None:
        """Test error when z_threshold is negative."""
        with pytest.raises(ValidationError, match="z_threshold must be positive"):
            validate_parameters(
                window_size=1000000,
                step_size=250000,
                min_dp=10,
                max_dp=200,
                z_threshold=-1.0,
            )

    def test_step_size_equal_to_window_size(self) -> None:
        """Test that step_size equal to window_size is allowed."""
        # Non-overlapping windows should be valid
        validate_parameters(
            window_size=1000000,
            step_size=1000000,
            min_dp=10,
            max_dp=200,
            z_threshold=3.0,
        )


class TestValidateOutputPath:
    """Tests for validate_output_path function."""

    def test_valid_output_in_existing_dir(self, output_dir: Path) -> None:
        """Test validation of output path in existing directory."""
        out_prefix = output_dir / "my_analysis"
        result = validate_output_path(out_prefix)

        assert result.parent == output_dir
        assert result.name == "my_analysis"

    def test_creates_parent_directory(self, tmp_path: Path) -> None:
        """Test that parent directory is created if it doesn't exist."""
        new_dir = tmp_path / "new_subdir" / "another"
        out_prefix = new_dir / "analysis"

        result = validate_output_path(out_prefix)

        assert new_dir.exists()
        assert result.parent == new_dir

    def test_nested_directory_creation(self, tmp_path: Path) -> None:
        """Test creation of deeply nested directories."""
        nested = tmp_path / "a" / "b" / "c" / "d"
        out_prefix = nested / "output"

        result = validate_output_path(out_prefix)

        assert nested.exists()
        assert result.parent == nested

    def test_returns_resolved_path(self, tmp_path: Path) -> None:
        """Test that returned path is resolved (absolute)."""
        out_prefix = tmp_path / "output"
        result = validate_output_path(out_prefix)

        assert result.is_absolute()

    def test_output_in_current_directory(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Test output path in current working directory."""
        monkeypatch.chdir(tmp_path)
        result = validate_output_path("my_output")

        assert result.is_absolute()
        assert result.name == "my_output"


class TestValidationError:
    """Tests for ValidationError exception."""

    def test_exception_message(self) -> None:
        """Test that ValidationError contains the correct message."""
        error = ValidationError("Test error message")
        assert str(error) == "Test error message"

    def test_exception_can_be_raised(self) -> None:
        """Test that ValidationError can be raised and caught."""
        with pytest.raises(ValidationError) as exc_info:
            raise ValidationError("Custom validation error")

        assert "Custom validation error" in str(exc_info.value)

    def test_exception_inheritance(self) -> None:
        """Test that ValidationError inherits from Exception."""
        assert issubclass(ValidationError, Exception)
