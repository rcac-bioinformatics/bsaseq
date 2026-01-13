"""Input validation utilities for bsaseq.

This module provides functions for validating VCF files, sample names,
and analysis parameters before running the analysis pipeline.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from bsaseq.utils.logging import get_logger

logger = get_logger(__name__)


class ValidationError(Exception):
    """Raised when input validation fails."""

    pass


def validate_vcf(vcf_path: str | Path) -> dict:
    """Validate VCF file and return metadata.

    Checks:
    - File exists and is readable
    - Valid VCF format
    - Contains required FORMAT fields (GT, AD)
    - Has at least 2 samples

    Args:
        vcf_path: Path to VCF file.

    Returns:
        dict with keys: samples, has_ad, has_gq, contigs

    Raises:
        ValidationError: If validation fails.
    """
    from cyvcf2 import VCF

    vcf_path = Path(vcf_path)

    # Check file exists
    if not vcf_path.exists():
        raise ValidationError(f"VCF file not found: {vcf_path}")

    # Try to open VCF
    try:
        vcf = VCF(str(vcf_path))
    except Exception as e:
        raise ValidationError(f"Cannot read VCF file: {e}") from e

    # Get samples
    samples = list(vcf.samples)
    if len(samples) < 2:
        raise ValidationError(
            f"VCF must contain at least 2 samples for BSA analysis. "
            f"Found {len(samples)} sample(s)."
        )

    # Check for required FORMAT fields
    has_ad = False
    has_gq = False
    has_gt = False

    for header in vcf.header_iter():
        if header["HeaderType"] == "FORMAT":
            field_id = header["ID"]
            if field_id == "AD":
                has_ad = True
            elif field_id == "GQ":
                has_gq = True
            elif field_id == "GT":
                has_gt = True

    if not has_gt:
        raise ValidationError(
            "VCF does not contain GT (genotype) FORMAT field. "
            "This is required for variant analysis."
        )

    if not has_ad:
        raise ValidationError(
            "VCF does not contain AD (allelic depth) FORMAT field.\n\n"
            "The AD field is required to calculate allele frequencies. Ensure your VCF "
            "was generated with a variant caller that outputs per-allele read depths "
            "(e.g., GATK HaplotypeCaller, bcftools mpileup)."
        )

    # Get contigs
    contigs = []
    for header in vcf.header_iter():
        if header["HeaderType"] == "CONTIG":
            contigs.append(header["ID"])

    vcf.close()

    return {
        "samples": samples,
        "has_ad": has_ad,
        "has_gq": has_gq,
        "contigs": contigs,
    }


def validate_samples(
    vcf_path: str | Path,
    high_bulk: str | Sequence[str],
    low_bulk: str | Sequence[str],
) -> tuple[list[str], list[str]]:
    """Validate sample names exist in VCF.

    Args:
        vcf_path: Path to VCF file.
        high_bulk: Sample name(s) for high bulk.
        low_bulk: Sample name(s) for low bulk.

    Returns:
        Tuple of (high_bulk_list, low_bulk_list) with validated names.

    Raises:
        ValidationError: If any sample not found or samples overlap.
    """
    from cyvcf2 import VCF

    # Parse sample names
    if isinstance(high_bulk, str):
        high_list = [s.strip() for s in high_bulk.split(",")]
    else:
        high_list = [s.strip() for s in high_bulk]

    if isinstance(low_bulk, str):
        low_list = [s.strip() for s in low_bulk.split(",")]
    else:
        low_list = [s.strip() for s in low_bulk]

    # Get available samples
    vcf = VCF(str(vcf_path))
    available = set(vcf.samples)
    vcf.close()

    # Check for overlap
    overlap = set(high_list) & set(low_list)
    if overlap:
        raise ValidationError(
            f"Sample(s) cannot be in both high and low bulk: {', '.join(overlap)}"
        )

    # Check high bulk samples
    missing_high = [s for s in high_list if s not in available]
    if missing_high:
        raise ValidationError(
            f"High bulk sample(s) not found in VCF: {', '.join(missing_high)}\n\n"
            f"Available samples: {', '.join(sorted(available))}\n\n"
            f"Hint: Sample names are case-sensitive. "
            f"Use 'bsaseq samples --vcf {vcf_path}' to list all samples."
        )

    # Check low bulk samples
    missing_low = [s for s in low_list if s not in available]
    if missing_low:
        raise ValidationError(
            f"Low bulk sample(s) not found in VCF: {', '.join(missing_low)}\n\n"
            f"Available samples: {', '.join(sorted(available))}\n\n"
            f"Hint: Sample names are case-sensitive. "
            f"Use 'bsaseq samples --vcf {vcf_path}' to list all samples."
        )

    return high_list, low_list


def validate_parameters(
    window_size: int,
    step_size: int,
    min_dp: int,
    max_dp: int,
    z_threshold: float,
) -> None:
    """Validate analysis parameters.

    Checks:
    - window_size > 0
    - step_size > 0 and <= window_size
    - min_dp > 0 and < max_dp
    - z_threshold > 0

    Args:
        window_size: Window size in bp.
        step_size: Step size in bp.
        min_dp: Minimum read depth.
        max_dp: Maximum read depth.
        z_threshold: Z-score threshold.

    Raises:
        ValidationError: If any parameter is invalid.
    """
    if window_size <= 0:
        raise ValidationError(f"window_size must be positive, got {window_size}")

    if step_size <= 0:
        raise ValidationError(f"step_size must be positive, got {step_size}")

    if step_size > window_size:
        raise ValidationError(
            f"step_size ({step_size:,}) cannot be larger than "
            f"window_size ({window_size:,})"
        )

    if min_dp <= 0:
        raise ValidationError(f"min_dp must be positive, got {min_dp}")

    if max_dp <= 0:
        raise ValidationError(f"max_dp must be positive, got {max_dp}")

    if min_dp >= max_dp:
        raise ValidationError(
            f"min_dp ({min_dp}) must be less than max_dp ({max_dp}). "
            f"Check your depth filter settings."
        )

    if z_threshold <= 0:
        raise ValidationError(f"z_threshold must be positive, got {z_threshold}")


def validate_output_path(out_prefix: str | Path) -> Path:
    """Validate output path and create parent directories.

    Checks:
    - Parent directory exists or can be created
    - Path is writable

    Args:
        out_prefix: Output prefix path.

    Returns:
        Resolved Path object.

    Raises:
        ValidationError: If path is invalid or not writable.
    """
    out_path = Path(out_prefix).resolve()
    parent = out_path.parent

    # Handle current directory case
    if str(parent) == ".":
        parent = Path.cwd()

    # Try to create parent directory
    if not parent.exists():
        try:
            parent.mkdir(parents=True, exist_ok=True)
            logger.debug(f"Created output directory: {parent}")
        except PermissionError as e:
            raise ValidationError(
                f"Cannot create output directory: {parent}\n"
                f"Permission denied."
            ) from e
        except OSError as e:
            raise ValidationError(
                f"Cannot create output directory: {parent}\n"
                f"Error: {e}"
            ) from e

    # Check if writable
    if not parent.is_dir():
        raise ValidationError(f"Output path parent is not a directory: {parent}")

    # Try to test write access
    test_file = parent / f".bsaseq_write_test_{out_path.name}"
    try:
        test_file.touch()
        test_file.unlink()
    except PermissionError as e:
        raise ValidationError(
            f"Output directory is not writable: {parent}\n"
            f"Check file permissions."
        ) from e
    except OSError as e:
        raise ValidationError(
            f"Cannot write to output directory: {parent}\n"
            f"Error: {e}"
        ) from e

    return out_path
