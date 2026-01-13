"""Pytest configuration and fixtures for bsaseq tests."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from bsaseq.core.models import Variant, Window
from bsaseq.analysis.candidates import CandidateRegion, CandidateVariant

if TYPE_CHECKING:
    pass


# Sample VCF content for testing
SAMPLE_VCF_CONTENT = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##contig=<ID=chr1,length=10000000>
##contig=<ID=chr2,length=8000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HIGH_BULK	LOW_BULK
chr1	1000	.	A	G	50.0	PASS	.	GT:AD:DP:GQ	0/1:20,30:50:99	0/1:40,10:50:99
chr1	2000	.	C	T	45.0	PASS	.	GT:AD:DP:GQ	0/1:15,35:50:99	0/1:35,15:50:99
chr1	3000	.	G	A	55.0	PASS	.	GT:AD:DP:GQ	0/1:25,25:50:99	0/1:30,20:50:99
chr2	1500	.	T	C	60.0	PASS	.	GT:AD:DP:GQ	0/1:10,40:50:99	0/1:45,5:50:99
"""

# VCF with various edge cases for filtering tests
SAMPLE_VCF_WITH_FILTERS = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##contig=<ID=chr1,length=10000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HIGH_BULK	LOW_BULK
chr1	1000	.	A	G	50.0	PASS	.	GT:AD:DP:GQ	0/1:20,30:50:99	0/1:40,10:50:99
chr1	2000	.	C	T	10.0	PASS	.	GT:AD:DP:GQ	0/1:20,30:50:99	0/1:40,10:50:99
chr1	3000	.	G	A	50.0	PASS	.	GT:AD:DP:GQ	0/1:2,3:5:99	0/1:40,10:50:99
chr1	4000	.	T	C	50.0	PASS	.	GT:AD:DP:GQ	0/1:20,30:50:5	0/1:40,10:50:99
chr1	5000	.	A	C,G	50.0	PASS	.	GT:AD:DP:GQ	0/1:20,30:50:99	0/1:40,10:50:99
chr1	6000	.	AT	A	50.0	PASS	.	GT:AD:DP:GQ	0/1:20,30:50:99	0/1:40,10:50:99
"""

# VCF with multiple samples for multi-bulk testing
MULTI_SAMPLE_VCF_CONTENT = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##contig=<ID=chr1,length=10000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	MUT1	MUT2	WT1	WT2
chr1	1000	.	A	G	50.0	PASS	.	GT:AD:DP:GQ	0/1:5,45:50:99	0/1:3,47:50:99	0/1:40,10:50:99	0/1:42,8:50:99
chr1	2000	.	C	T	50.0	PASS	.	GT:AD:DP:GQ	0/1:10,40:50:99	0/1:8,42:50:99	0/1:35,15:50:99	0/1:38,12:50:99
chr1	3000	.	G	A	50.0	PASS	.	GT:AD:DP:GQ	0/1:15,35:50:99	0/1:12,38:50:99	0/1:30,20:50:99	0/1:32,18:50:99
"""


@pytest.fixture
def sample_vcf_path(tmp_path: Path) -> Path:
    """Create a temporary VCF file with sample data.

    Returns:
        Path to the temporary VCF file.
    """
    vcf_path = tmp_path / "test.vcf"
    vcf_path.write_text(SAMPLE_VCF_CONTENT)
    return vcf_path


@pytest.fixture
def sample_vcf_with_filters_path(tmp_path: Path) -> Path:
    """Create a temporary VCF file with edge cases for filter testing.

    Returns:
        Path to the temporary VCF file.
    """
    vcf_path = tmp_path / "test_filters.vcf"
    vcf_path.write_text(SAMPLE_VCF_WITH_FILTERS)
    return vcf_path


@pytest.fixture
def multi_sample_vcf_path(tmp_path: Path) -> Path:
    """Create a temporary VCF file with multiple samples per bulk.

    Returns:
        Path to the temporary VCF file.
    """
    vcf_path = tmp_path / "multi_sample.vcf"
    vcf_path.write_text(MULTI_SAMPLE_VCF_CONTENT)
    return vcf_path


@pytest.fixture
def empty_vcf_path(tmp_path: Path) -> Path:
    """Create a temporary VCF file with no variants.

    Returns:
        Path to the temporary VCF file.
    """
    content = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HIGH_BULK	LOW_BULK
"""
    vcf_path = tmp_path / "empty.vcf"
    vcf_path.write_text(content)
    return vcf_path


@pytest.fixture
def sample_variants() -> list[Variant]:
    """Create a list of sample variants for testing.

    Returns:
        List of Variant objects.
    """
    return [
        Variant(
            chrom="chr1",
            pos=1000,
            ref="A",
            alt="G",
            af_high=0.6,
            af_low=0.2,
            dp_high=50,
            dp_low=50,
        ),
        Variant(
            chrom="chr1",
            pos=2000,
            ref="C",
            alt="T",
            af_high=0.7,
            af_low=0.3,
            dp_high=50,
            dp_low=50,
        ),
        Variant(
            chrom="chr1",
            pos=3000,
            ref="G",
            alt="A",
            af_high=0.5,
            af_low=0.4,
            dp_high=50,
            dp_low=50,
        ),
        Variant(
            chrom="chr1",
            pos=4000,
            ref="T",
            alt="C",
            af_high=0.8,
            af_low=0.1,
            dp_high=50,
            dp_low=50,
        ),
        Variant(
            chrom="chr1",
            pos=5000,
            ref="A",
            alt="G",
            af_high=0.65,
            af_low=0.25,
            dp_high=50,
            dp_low=50,
        ),
        Variant(
            chrom="chr2",
            pos=1000,
            ref="C",
            alt="T",
            af_high=0.75,
            af_low=0.15,
            dp_high=50,
            dp_low=50,
        ),
        Variant(
            chrom="chr2",
            pos=2000,
            ref="G",
            alt="A",
            af_high=0.55,
            af_low=0.35,
            dp_high=50,
            dp_low=50,
        ),
    ]


@pytest.fixture
def sample_windows() -> list[Window]:
    """Create a list of sample windows for testing.

    Returns:
        List of Window objects spanning chr1 and chr2.
    """
    return [
        Window(
            chrom="chr1",
            start=1,
            end=1000000,
            n_variants=50,
            mean_delta_af=0.3,
            median_delta_af=0.28,
            mean_af_high=0.65,
            mean_af_low=0.35,
            tricube_delta_af=0.32,
            g_statistic=45.2,
            p_value=1.8e-11,
            z_score=2.1,
        ),
        Window(
            chrom="chr1",
            start=250001,
            end=1250000,
            n_variants=65,
            mean_delta_af=0.45,
            median_delta_af=0.42,
            mean_af_high=0.72,
            mean_af_low=0.27,
            tricube_delta_af=0.48,
            g_statistic=78.5,
            p_value=8.3e-19,
            z_score=3.5,
        ),
        Window(
            chrom="chr1",
            start=500001,
            end=1500000,
            n_variants=72,
            mean_delta_af=0.52,
            median_delta_af=0.50,
            mean_af_high=0.78,
            mean_af_low=0.26,
            tricube_delta_af=0.55,
            g_statistic=92.1,
            p_value=1.2e-21,
            z_score=4.2,
        ),
        Window(
            chrom="chr1",
            start=750001,
            end=1750000,
            n_variants=58,
            mean_delta_af=0.38,
            median_delta_af=0.36,
            mean_af_high=0.69,
            mean_af_low=0.31,
            tricube_delta_af=0.40,
            g_statistic=62.3,
            p_value=2.9e-15,
            z_score=2.8,
        ),
        Window(
            chrom="chr2",
            start=1,
            end=1000000,
            n_variants=45,
            mean_delta_af=0.15,
            median_delta_af=0.14,
            mean_af_high=0.55,
            mean_af_low=0.40,
            tricube_delta_af=0.16,
            g_statistic=18.5,
            p_value=1.7e-5,
            z_score=0.8,
        ),
    ]


@pytest.fixture
def sample_regions() -> list[CandidateRegion]:
    """Create a list of sample candidate regions for testing.

    Returns:
        List of CandidateRegion objects.
    """
    return [
        CandidateRegion(
            chrom="chr1",
            start=250001,
            end=1750000,
            n_windows=3,
            n_variants=195,
            max_z_score=4.2,
            mean_z_score=3.5,
            max_delta_af=0.52,
            mean_delta_af=0.45,
            peak_position=750000,
        ),
    ]


@pytest.fixture
def sample_candidates(sample_variants: list[Variant], sample_regions: list[CandidateRegion]) -> list[CandidateVariant]:
    """Create a list of sample candidate variants for testing.

    Args:
        sample_variants: Fixture providing sample variants.
        sample_regions: Fixture providing sample regions.

    Returns:
        List of CandidateVariant objects.
    """
    region = sample_regions[0]
    candidates = []
    for v in sample_variants:
        if region.contains(v.chrom, v.pos):
            candidates.append(CandidateVariant.from_variant(v, region, region_rank=1))
    return candidates


# VCF for validation testing (missing required fields)
VCF_MISSING_AD = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##contig=<ID=chr1,length=10000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	1000	.	A	G	50.0	PASS	.	GT:DP	0/1:50	0/1:50
"""

VCF_MISSING_GT = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##contig=<ID=chr1,length=10000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1	SAMPLE2
chr1	1000	.	A	G	50.0	PASS	.	AD:DP	20,30:50	40,10:50
"""

VCF_SINGLE_SAMPLE = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##contig=<ID=chr1,length=10000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ONLY_SAMPLE
chr1	1000	.	A	G	50.0	PASS	.	GT:AD	0/1:20,30
"""


@pytest.fixture
def vcf_missing_ad_path(tmp_path: Path) -> Path:
    """Create a VCF file missing the AD field.

    Returns:
        Path to the temporary VCF file.
    """
    vcf_path = tmp_path / "missing_ad.vcf"
    vcf_path.write_text(VCF_MISSING_AD)
    return vcf_path


@pytest.fixture
def vcf_missing_gt_path(tmp_path: Path) -> Path:
    """Create a VCF file missing the GT field.

    Returns:
        Path to the temporary VCF file.
    """
    vcf_path = tmp_path / "missing_gt.vcf"
    vcf_path.write_text(VCF_MISSING_GT)
    return vcf_path


@pytest.fixture
def vcf_single_sample_path(tmp_path: Path) -> Path:
    """Create a VCF file with only one sample.

    Returns:
        Path to the temporary VCF file.
    """
    vcf_path = tmp_path / "single_sample.vcf"
    vcf_path.write_text(VCF_SINGLE_SAMPLE)
    return vcf_path


@pytest.fixture
def output_dir(tmp_path: Path) -> Path:
    """Create a temporary output directory.

    Returns:
        Path to the temporary output directory.
    """
    out_dir = tmp_path / "output"
    out_dir.mkdir()
    return out_dir


@pytest.fixture
def readonly_dir(tmp_path: Path) -> Path:
    """Create a read-only directory for testing permission errors.

    Returns:
        Path to the read-only directory.
    """
    ro_dir = tmp_path / "readonly"
    ro_dir.mkdir()
    ro_dir.chmod(0o444)
    yield ro_dir
    # Restore permissions for cleanup
    ro_dir.chmod(0o755)
