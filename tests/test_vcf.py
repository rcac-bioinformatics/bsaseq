"""Tests for VCF parsing functionality."""

from __future__ import annotations

from pathlib import Path

import pytest

from bsaseq.core.models import Variant
from bsaseq.io.vcf import get_sample_names, parse_vcf


class TestGetSampleNames:
    """Tests for get_sample_names function."""

    def test_returns_sample_names(self, sample_vcf_path: Path) -> None:
        """Test that sample names are correctly extracted from VCF."""
        samples = get_sample_names(sample_vcf_path)
        assert samples == ["HIGH_BULK", "LOW_BULK"]

    def test_empty_vcf_has_samples(self, empty_vcf_path: Path) -> None:
        """Test that sample names are extracted even from VCF with no variants."""
        samples = get_sample_names(empty_vcf_path)
        assert samples == ["HIGH_BULK", "LOW_BULK"]


class TestParseVcf:
    """Tests for parse_vcf function."""

    def test_parses_variants(self, sample_vcf_path: Path) -> None:
        """Test that variants are correctly parsed from VCF."""
        variants = list(
            parse_vcf(
                vcf_path=sample_vcf_path,
                high_bulk="HIGH_BULK",
                low_bulk="LOW_BULK",
                min_dp=10,
                max_dp=200,
                min_gq=20,
                min_qual=30.0,
            )
        )

        assert len(variants) == 4

        # Check first variant
        v1 = variants[0]
        assert v1.chrom == "chr1"
        assert v1.pos == 1000
        assert v1.ref == "A"
        assert v1.alt == "G"
        assert v1.dp_high == 50
        assert v1.dp_low == 50

    def test_calculates_allele_frequencies(self, sample_vcf_path: Path) -> None:
        """Test that allele frequencies are correctly calculated."""
        variants = list(
            parse_vcf(
                vcf_path=sample_vcf_path,
                high_bulk="HIGH_BULK",
                low_bulk="LOW_BULK",
            )
        )

        # First variant: HIGH_BULK has 20 ref, 30 alt -> AF = 30/50 = 0.6
        # LOW_BULK has 40 ref, 10 alt -> AF = 10/50 = 0.2
        v1 = variants[0]
        assert v1.af_high == pytest.approx(0.6, rel=1e-6)
        assert v1.af_low == pytest.approx(0.2, rel=1e-6)
        assert v1.delta_af == pytest.approx(0.4, rel=1e-6)

    def test_filters_low_qual(self, sample_vcf_with_filters_path: Path) -> None:
        """Test that low QUAL variants are filtered out."""
        variants = list(
            parse_vcf(
                vcf_path=sample_vcf_with_filters_path,
                high_bulk="HIGH_BULK",
                low_bulk="LOW_BULK",
                min_qual=30.0,
            )
        )

        # Second variant has QUAL=10, should be filtered
        positions = [v.pos for v in variants]
        assert 2000 not in positions

    def test_filters_low_depth(self, sample_vcf_with_filters_path: Path) -> None:
        """Test that low depth variants are filtered out."""
        variants = list(
            parse_vcf(
                vcf_path=sample_vcf_with_filters_path,
                high_bulk="HIGH_BULK",
                low_bulk="LOW_BULK",
                min_dp=10,
            )
        )

        # Third variant has DP=5 in high bulk, should be filtered
        positions = [v.pos for v in variants]
        assert 3000 not in positions

    def test_filters_low_gq(self, sample_vcf_with_filters_path: Path) -> None:
        """Test that low GQ variants are filtered out."""
        variants = list(
            parse_vcf(
                vcf_path=sample_vcf_with_filters_path,
                high_bulk="HIGH_BULK",
                low_bulk="LOW_BULK",
                min_gq=20,
            )
        )

        # Fourth variant has GQ=5 in high bulk, should be filtered
        positions = [v.pos for v in variants]
        assert 4000 not in positions

    def test_filters_multiallelic(self, sample_vcf_with_filters_path: Path) -> None:
        """Test that multiallelic variants are filtered out."""
        variants = list(
            parse_vcf(
                vcf_path=sample_vcf_with_filters_path,
                high_bulk="HIGH_BULK",
                low_bulk="LOW_BULK",
                min_qual=0,
                min_dp=1,
                min_gq=0,
            )
        )

        # Fifth variant is multiallelic (A>C,G), should be filtered
        positions = [v.pos for v in variants]
        assert 5000 not in positions

    def test_filters_indels(self, sample_vcf_with_filters_path: Path) -> None:
        """Test that indels are filtered out."""
        variants = list(
            parse_vcf(
                vcf_path=sample_vcf_with_filters_path,
                high_bulk="HIGH_BULK",
                low_bulk="LOW_BULK",
                min_qual=0,
                min_dp=1,
                min_gq=0,
            )
        )

        # Sixth variant is an indel (AT>A), should be filtered
        positions = [v.pos for v in variants]
        assert 6000 not in positions

    def test_raises_on_missing_sample(self, sample_vcf_path: Path) -> None:
        """Test that ValueError is raised for missing sample names."""
        with pytest.raises(ValueError, match="not found in VCF"):
            list(
                parse_vcf(
                    vcf_path=sample_vcf_path,
                    high_bulk="NONEXISTENT",
                    low_bulk="LOW_BULK",
                )
            )

    def test_raises_on_missing_file(self, tmp_path: Path) -> None:
        """Test that FileNotFoundError is raised for missing files."""
        with pytest.raises(FileNotFoundError):
            list(
                parse_vcf(
                    vcf_path=tmp_path / "nonexistent.vcf",
                    high_bulk="HIGH_BULK",
                    low_bulk="LOW_BULK",
                )
            )

    def test_returns_iterator(self, sample_vcf_path: Path) -> None:
        """Test that parse_vcf returns an iterator (memory efficient)."""
        result = parse_vcf(
            vcf_path=sample_vcf_path,
            high_bulk="HIGH_BULK",
            low_bulk="LOW_BULK",
        )

        # Should be an iterator, not a list
        assert hasattr(result, "__iter__")
        assert hasattr(result, "__next__")

    def test_empty_vcf_returns_no_variants(self, empty_vcf_path: Path) -> None:
        """Test that empty VCF returns no variants."""
        variants = list(
            parse_vcf(
                vcf_path=empty_vcf_path,
                high_bulk="HIGH_BULK",
                low_bulk="LOW_BULK",
            )
        )

        assert len(variants) == 0


class TestVariantProperties:
    """Tests for Variant dataclass properties."""

    def test_delta_af_calculation(self) -> None:
        """Test delta_af property calculation."""
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="G",
            af_high=0.8,
            af_low=0.2,
            dp_high=50,
            dp_low=50,
        )

        assert variant.delta_af == pytest.approx(0.6, rel=1e-6)

    def test_snp_index_properties(self) -> None:
        """Test SNP index properties."""
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="G",
            af_high=0.7,
            af_low=0.3,
            dp_high=50,
            dp_low=50,
        )

        assert variant.snp_index_high == pytest.approx(0.7, rel=1e-6)
        assert variant.snp_index_low == pytest.approx(0.3, rel=1e-6)

    def test_negative_delta_af(self) -> None:
        """Test that delta_af can be negative."""
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="G",
            af_high=0.2,
            af_low=0.8,
            dp_high=50,
            dp_low=50,
        )

        assert variant.delta_af == pytest.approx(-0.6, rel=1e-6)

    def test_repr(self) -> None:
        """Test string representation of Variant."""
        variant = Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="G",
            af_high=0.6,
            af_low=0.2,
            dp_high=50,
            dp_low=50,
        )

        repr_str = repr(variant)
        assert "chr1:1000" in repr_str
        assert "A>G" in repr_str
