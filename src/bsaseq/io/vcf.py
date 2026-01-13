"""VCF parsing and writing utilities for bsaseq.

This module provides functions for parsing VCF files and extracting
variant information for bulk segregant analysis, as well as writing
candidate variants to VCF format for annotation.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from cyvcf2 import VCF

from bsaseq.core.models import Variant
from bsaseq.core.stats import calculate_allele_frequency
from bsaseq.utils.logging import get_logger, print_stats

if TYPE_CHECKING:
    from collections.abc import Iterator

    from bsaseq.analysis.candidates import CandidateVariant

logger = get_logger(__name__)


class FilterStats:
    """Track filtering statistics during VCF parsing."""

    def __init__(self) -> None:
        """Initialize filter statistics."""
        self.total_records = 0
        self.not_snp = 0
        self.multiallelic = 0
        self.low_qual = 0
        self.missing_ad = 0
        self.depth_filter_high = 0
        self.depth_filter_low = 0
        self.gq_filter_high = 0
        self.gq_filter_low = 0
        self.passed = 0

    def to_dict(self) -> dict[str, int]:
        """Convert stats to dictionary for display."""
        return {
            "Total records": self.total_records,
            "Not SNP (filtered)": self.not_snp,
            "Multiallelic (filtered)": self.multiallelic,
            "Low QUAL (filtered)": self.low_qual,
            "Missing AD (filtered)": self.missing_ad,
            "Depth filter (high bulk)": self.depth_filter_high,
            "Depth filter (low bulk)": self.depth_filter_low,
            "GQ filter (high bulk)": self.gq_filter_high,
            "GQ filter (low bulk)": self.gq_filter_low,
            "Variants passed": self.passed,
        }


def get_sample_names(vcf_path: str | Path) -> list[str]:
    """Get list of sample names from a VCF file.

    Args:
        vcf_path: Path to VCF file (can be gzipped).

    Returns:
        List of sample names in the VCF.
    """
    vcf = VCF(str(vcf_path))
    return list(vcf.samples)


def _parse_sample_names(samples: str | list[str]) -> list[str]:
    """Parse sample names from string or list.

    Accepts either a single sample name, comma-separated string,
    or list of sample names.

    Args:
        samples: Sample name(s) as string or list.

    Returns:
        List of sample names.

    Examples:
        >>> _parse_sample_names("sample1")
        ['sample1']
        >>> _parse_sample_names("sample1,sample2")
        ['sample1', 'sample2']
        >>> _parse_sample_names(["sample1", "sample2"])
        ['sample1', 'sample2']
    """
    if isinstance(samples, list):
        return [s.strip() for s in samples]
    return [s.strip() for s in samples.split(",")]


def _validate_sample_names(
    requested: list[str],
    available: list[str],
    bulk_name: str,
) -> list[int]:
    """Validate sample names and return their indices.

    Args:
        requested: List of requested sample names.
        available: List of available sample names in VCF.
        bulk_name: Name of the bulk (for error messages).

    Returns:
        List of indices for the requested samples.

    Raises:
        ValueError: If any requested samples are not found.
    """
    missing = [s for s in requested if s not in available]
    if missing:
        raise ValueError(
            f"{bulk_name} sample(s) not found in VCF: {', '.join(missing)}. "
            f"Available samples: {', '.join(available)}"
        )
    return [available.index(s) for s in requested]


def parse_vcf(
    vcf_path: str | Path,
    high_bulk: str | list[str],
    low_bulk: str | list[str],
    min_dp: int = 10,
    max_dp: int = 200,
    min_gq: int = 20,
    min_qual: float = 30.0,
) -> Iterator[Variant]:
    """Parse VCF and yield Variant objects.

    Extracts biallelic SNPs from a VCF file, calculating allele
    frequencies for both bulk samples. Supports multiple samples
    per bulk, with allele frequencies calculated by pooling read
    counts across samples.

    Args:
        vcf_path: Path to VCF file (can be gzipped).
        high_bulk: Sample name(s) for high bulk (e.g., mutant pool).
            Can be a single name, comma-separated string ("mut1,mut2"),
            or list of names.
        low_bulk: Sample name(s) for low bulk (e.g., wild-type pool).
            Can be a single name, comma-separated string, or list.
        min_dp: Minimum read depth per bulk (pooled across samples).
        max_dp: Maximum read depth per bulk (pooled across samples).
        min_gq: Minimum genotype quality. ALL samples in a bulk must
            pass this threshold.
        min_qual: Minimum variant QUAL score.

    Yields:
        Variant objects passing all filters.

    Raises:
        ValueError: If sample names not found in VCF.
        FileNotFoundError: If VCF file does not exist.

    Example:
        >>> # Single sample per bulk
        >>> variants = parse_vcf("data.vcf", "MUT", "WT")

        >>> # Multiple samples per bulk (comma-separated)
        >>> variants = parse_vcf("data.vcf", "mut1,mut2", "wt1,wt2")

        >>> # Multiple samples per bulk (list)
        >>> variants = parse_vcf("data.vcf", ["mut1", "mut2"], ["wt1", "wt2"])
    """
    vcf_path = Path(vcf_path)

    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    logger.info(f"Opening VCF file: {vcf_path}")

    vcf = VCF(str(vcf_path))
    available_samples = vcf.samples

    # Parse and validate sample names
    high_samples = _parse_sample_names(high_bulk)
    low_samples = _parse_sample_names(low_bulk)

    high_indices = _validate_sample_names(high_samples, available_samples, "High bulk")
    low_indices = _validate_sample_names(low_samples, available_samples, "Low bulk")

    n_samples_high = len(high_samples)
    n_samples_low = len(low_samples)

    logger.info(f"High bulk samples ({n_samples_high}): {', '.join(high_samples)}")
    logger.info(f"Low bulk samples ({n_samples_low}): {', '.join(low_samples)}")
    logger.info(
        f"Filters: min_dp={min_dp}, max_dp={max_dp}, "
        f"min_gq={min_gq}, min_qual={min_qual}"
    )

    stats = FilterStats()

    for record in vcf:
        stats.total_records += 1

        # Skip non-SNPs
        if not record.is_snp:
            stats.not_snp += 1
            continue

        # Skip multiallelic sites (keep biallelic only)
        if len(record.ALT) != 1:
            stats.multiallelic += 1
            continue

        # Check QUAL score
        if record.QUAL is not None and record.QUAL < min_qual:
            stats.low_qual += 1
            continue

        # Get AD (allelic depth) for all samples
        ad = record.format("AD")
        if ad is None:
            stats.missing_ad += 1
            continue

        # Pool allelic depths across samples in each bulk
        try:
            # High bulk pooling
            ref_dp_high = 0
            alt_dp_high = 0
            for idx in high_indices:
                ad_sample = ad[idx]
                # Check for missing values (-1 in cyvcf2)
                if ad_sample[0] < 0 or ad_sample[1] < 0:
                    raise ValueError("Missing AD value")
                ref_dp_high += int(ad_sample[0])
                alt_dp_high += int(ad_sample[1])

            # Low bulk pooling
            ref_dp_low = 0
            alt_dp_low = 0
            for idx in low_indices:
                ad_sample = ad[idx]
                if ad_sample[0] < 0 or ad_sample[1] < 0:
                    raise ValueError("Missing AD value")
                ref_dp_low += int(ad_sample[0])
                alt_dp_low += int(ad_sample[1])

        except (IndexError, TypeError, ValueError):
            stats.missing_ad += 1
            continue

        # Calculate pooled total depths
        dp_high = ref_dp_high + alt_dp_high
        dp_low = ref_dp_low + alt_dp_low

        # Apply depth filters (on pooled depths)
        if dp_high < min_dp or dp_high > max_dp:
            stats.depth_filter_high += 1
            continue
        if dp_low < min_dp or dp_low > max_dp:
            stats.depth_filter_low += 1
            continue

        # Get GQ (genotype quality) if available
        # ALL samples in a bulk must pass the GQ filter
        gq = record.format("GQ")
        if gq is not None:
            try:
                # Check all high bulk samples
                gq_filter_failed = False
                for idx in high_indices:
                    gq_val = gq[idx]
                    if gq_val is not None and gq_val >= 0 and gq_val < min_gq:
                        gq_filter_failed = True
                        break

                if gq_filter_failed:
                    stats.gq_filter_high += 1
                    continue

                # Check all low bulk samples
                for idx in low_indices:
                    gq_val = gq[idx]
                    if gq_val is not None and gq_val >= 0 and gq_val < min_gq:
                        gq_filter_failed = True
                        break

                if gq_filter_failed:
                    stats.gq_filter_low += 1
                    continue

            except (IndexError, TypeError):
                # GQ not available for these samples, continue without filtering
                pass

        # Calculate pooled allele frequencies
        af_high = calculate_allele_frequency(ref_dp_high, alt_dp_high)
        af_low = calculate_allele_frequency(ref_dp_low, alt_dp_low)

        stats.passed += 1

        yield Variant(
            chrom=record.CHROM,
            pos=record.POS,
            ref=record.REF,
            alt=record.ALT[0],
            af_high=af_high,
            af_low=af_low,
            dp_high=dp_high,
            dp_low=dp_low,
            n_samples_high=n_samples_high,
            n_samples_low=n_samples_low,
        )

    # Log filtering statistics
    logger.info("VCF parsing complete")
    print_stats(stats.to_dict(), title="Filtering Statistics")


def parse_vcf_to_list(
    vcf_path: str | Path,
    high_bulk: str | list[str],
    low_bulk: str | list[str],
    min_dp: int = 10,
    max_dp: int = 200,
    min_gq: int = 20,
    min_qual: float = 30.0,
) -> list[Variant]:
    """Parse VCF and return list of Variant objects.

    Convenience wrapper around parse_vcf() that returns a list
    instead of an iterator. Use with caution for large VCF files
    as this loads all variants into memory.

    Args:
        vcf_path: Path to VCF file (can be gzipped).
        high_bulk: Sample name(s) for high bulk. Can be comma-separated.
        low_bulk: Sample name(s) for low bulk. Can be comma-separated.
        min_dp: Minimum read depth per bulk (pooled).
        max_dp: Maximum read depth per bulk (pooled).
        min_gq: Minimum genotype quality (all samples must pass).
        min_qual: Minimum variant QUAL score.

    Returns:
        List of Variant objects passing all filters.

    Raises:
        ValueError: If sample names not found in VCF.
        FileNotFoundError: If VCF file does not exist.
    """
    return list(
        parse_vcf(
            vcf_path=vcf_path,
            high_bulk=high_bulk,
            low_bulk=low_bulk,
            min_dp=min_dp,
            max_dp=max_dp,
            min_gq=min_gq,
            min_qual=min_qual,
        )
    )


def write_candidates_vcf(
    candidates: list[CandidateVariant],
    output_path: str | Path,
) -> Path:
    """Write candidate variants to VCF format for snpEff annotation.

    Creates a minimal VCF file containing only candidate variants,
    suitable for input to snpEff.

    Args:
        candidates: List of CandidateVariant objects.
        output_path: Output VCF file path.

    Returns:
        Path to written VCF file.
    """
    output_path = Path(output_path)

    logger.info(f"Writing {len(candidates)} candidate variants to VCF: {output_path}")

    with open(output_path, "w") as f:
        # Write VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=bsaseq\n")
        f.write('##INFO=<ID=DELTA_AF,Number=1,Type=Float,Description="Delta allele frequency between bulks">\n')
        f.write('##INFO=<ID=AF_HIGH,Number=1,Type=Float,Description="Allele frequency in high bulk">\n')
        f.write('##INFO=<ID=AF_LOW,Number=1,Type=Float,Description="Allele frequency in low bulk">\n')
        f.write('##INFO=<ID=DP_HIGH,Number=1,Type=Integer,Description="Read depth in high bulk">\n')
        f.write('##INFO=<ID=DP_LOW,Number=1,Type=Integer,Description="Read depth in low bulk">\n')
        f.write('##INFO=<ID=REGION_RANK,Number=1,Type=Integer,Description="Rank of containing candidate region">\n')
        f.write('##INFO=<ID=DIST_PEAK,Number=1,Type=Integer,Description="Distance to region peak">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Write variant records
        for candidate in candidates:
            info = (
                f"DELTA_AF={candidate.delta_af:.4f};"
                f"AF_HIGH={candidate.af_high:.4f};"
                f"AF_LOW={candidate.af_low:.4f};"
                f"DP_HIGH={candidate.dp_high};"
                f"DP_LOW={candidate.dp_low};"
                f"REGION_RANK={candidate.region_rank};"
                f"DIST_PEAK={candidate.distance_to_peak}"
            )
            f.write(
                f"{candidate.chrom}\t{candidate.pos}\t.\t"
                f"{candidate.ref}\t{candidate.alt}\t.\tPASS\t{info}\n"
            )

    logger.info(f"VCF written: {output_path}")
    return output_path
