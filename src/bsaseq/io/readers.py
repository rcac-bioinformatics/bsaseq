"""TSV file readers for bsaseq.

This module provides functions for reading analysis output files
back into Python objects, useful for re-generating plots from
existing analysis results.
"""

from __future__ import annotations

from pathlib import Path

from bsaseq.analysis.candidates import CandidateRegion, CandidateVariant
from bsaseq.core.models import Variant, Window
from bsaseq.utils.logging import get_logger

logger = get_logger(__name__)


def read_variants_tsv(path: str | Path) -> list[Variant]:
    """Read variants from TSV file.

    Expects columns:
        chrom, pos, ref, alt, af_high, af_low, delta_af,
        dp_high, dp_low, n_samples_high, n_samples_low

    Args:
        path: Path to variants TSV file.

    Returns:
        List of Variant objects.

    Raises:
        FileNotFoundError: If file doesn't exist.
        ValueError: If file format is invalid.
    """
    path = Path(path)
    logger.info(f"Reading variants from {path}")

    if not path.exists():
        raise FileNotFoundError(f"Variants file not found: {path}")

    variants = []

    with open(path) as f:
        # Read header
        header = f.readline().strip().split("\t")

        # Validate required columns
        required = ["chrom", "pos", "ref", "alt", "af_high", "af_low", "dp_high", "dp_low"]
        missing = [col for col in required if col not in header]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Get column indices
        col_idx = {col: i for i, col in enumerate(header)}

        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")

            try:
                variant = Variant(
                    chrom=fields[col_idx["chrom"]],
                    pos=int(fields[col_idx["pos"]]),
                    ref=fields[col_idx["ref"]],
                    alt=fields[col_idx["alt"]],
                    af_high=float(fields[col_idx["af_high"]]),
                    af_low=float(fields[col_idx["af_low"]]),
                    dp_high=int(fields[col_idx["dp_high"]]),
                    dp_low=int(fields[col_idx["dp_low"]]),
                    n_samples_high=int(fields[col_idx.get("n_samples_high", len(fields))])
                        if "n_samples_high" in col_idx else 1,
                    n_samples_low=int(fields[col_idx.get("n_samples_low", len(fields))])
                        if "n_samples_low" in col_idx else 1,
                )
                variants.append(variant)
            except (IndexError, ValueError) as e:
                logger.warning(f"Skipping invalid line {line_num}: {e}")
                continue

    logger.info(f"Read {len(variants):,} variants from {path}")
    return variants


def read_windows_tsv(path: str | Path) -> list[Window]:
    """Read windows from TSV file.

    Expects columns:
        chrom, start, end, midpoint, n_variants, mean_delta_af,
        median_delta_af, tricube_delta_af, mean_af_high, mean_af_low,
        g_statistic, p_value, z_score

    Args:
        path: Path to windows TSV file.

    Returns:
        List of Window objects.

    Raises:
        FileNotFoundError: If file doesn't exist.
        ValueError: If file format is invalid.
    """
    path = Path(path)
    logger.info(f"Reading windows from {path}")

    if not path.exists():
        raise FileNotFoundError(f"Windows file not found: {path}")

    windows = []

    with open(path) as f:
        # Read header
        header = f.readline().strip().split("\t")

        # Validate required columns
        required = ["chrom", "start", "end", "n_variants", "mean_delta_af"]
        missing = [col for col in required if col not in header]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Get column indices
        col_idx = {col: i for i, col in enumerate(header)}

        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")

            try:
                # Parse z_score (may be "NA")
                z_score_str = fields[col_idx.get("z_score", len(fields))] if "z_score" in col_idx else "NA"
                z_score = None if z_score_str.upper() == "NA" else float(z_score_str)

                window = Window(
                    chrom=fields[col_idx["chrom"]],
                    start=int(fields[col_idx["start"]]),
                    end=int(fields[col_idx["end"]]),
                    n_variants=int(fields[col_idx["n_variants"]]),
                    mean_delta_af=float(fields[col_idx["mean_delta_af"]]),
                    median_delta_af=float(fields[col_idx.get("median_delta_af", col_idx["mean_delta_af"])]),
                    tricube_delta_af=float(fields[col_idx.get("tricube_delta_af", col_idx["mean_delta_af"])]),
                    mean_af_high=float(fields[col_idx.get("mean_af_high", 0)]) if "mean_af_high" in col_idx else 0.5,
                    mean_af_low=float(fields[col_idx.get("mean_af_low", 0)]) if "mean_af_low" in col_idx else 0.5,
                    g_statistic=float(fields[col_idx.get("g_statistic", 0)]) if "g_statistic" in col_idx else 0.0,
                    p_value=float(fields[col_idx.get("p_value", 0)]) if "p_value" in col_idx else 1.0,
                    z_score=z_score,
                )
                windows.append(window)
            except (IndexError, ValueError) as e:
                logger.warning(f"Skipping invalid line {line_num}: {e}")
                continue

    logger.info(f"Read {len(windows):,} windows from {path}")
    return windows


def read_regions_tsv(path: str | Path) -> list[CandidateRegion]:
    """Read candidate regions from TSV file.

    Expects columns:
        rank, chrom, start, end, length, n_windows, n_variants,
        max_z_score, mean_z_score, max_delta_af, mean_delta_af, peak_position

    Args:
        path: Path to regions TSV file.

    Returns:
        List of CandidateRegion objects.

    Raises:
        FileNotFoundError: If file doesn't exist.
        ValueError: If file format is invalid.
    """
    path = Path(path)
    logger.info(f"Reading regions from {path}")

    if not path.exists():
        raise FileNotFoundError(f"Regions file not found: {path}")

    regions = []

    with open(path) as f:
        # Read header
        header = f.readline().strip().split("\t")

        # Validate required columns
        required = ["chrom", "start", "end"]
        missing = [col for col in required if col not in header]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Get column indices
        col_idx = {col: i for i, col in enumerate(header)}

        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")

            try:
                region = CandidateRegion(
                    chrom=fields[col_idx["chrom"]],
                    start=int(fields[col_idx["start"]]),
                    end=int(fields[col_idx["end"]]),
                    n_windows=int(fields[col_idx.get("n_windows", 0)]) if "n_windows" in col_idx else 1,
                    n_variants=int(fields[col_idx.get("n_variants", 0)]) if "n_variants" in col_idx else 0,
                    max_z_score=float(fields[col_idx.get("max_z_score", 0)]) if "max_z_score" in col_idx else 0.0,
                    mean_z_score=float(fields[col_idx.get("mean_z_score", 0)]) if "mean_z_score" in col_idx else 0.0,
                    max_delta_af=float(fields[col_idx.get("max_delta_af", 0)]) if "max_delta_af" in col_idx else 0.0,
                    mean_delta_af=float(fields[col_idx.get("mean_delta_af", 0)]) if "mean_delta_af" in col_idx else 0.0,
                    peak_position=int(fields[col_idx.get("peak_position", col_idx["start"])])
                        if "peak_position" in col_idx else int(fields[col_idx["start"]]),
                )
                regions.append(region)
            except (IndexError, ValueError) as e:
                logger.warning(f"Skipping invalid line {line_num}: {e}")
                continue

    logger.info(f"Read {len(regions):,} regions from {path}")
    return regions


def read_candidate_variants_tsv(path: str | Path) -> list[CandidateVariant]:
    """Read candidate variants from TSV file.

    Expects columns:
        chrom, pos, ref, alt, af_high, af_low, delta_af,
        dp_high, dp_low, region_rank, distance_to_peak

    Args:
        path: Path to candidate variants TSV file.

    Returns:
        List of CandidateVariant objects.

    Raises:
        FileNotFoundError: If file doesn't exist.
        ValueError: If file format is invalid.
    """
    path = Path(path)
    logger.info(f"Reading candidate variants from {path}")

    if not path.exists():
        raise FileNotFoundError(f"Candidate variants file not found: {path}")

    variants = []

    with open(path) as f:
        # Read header
        header = f.readline().strip().split("\t")

        # Validate required columns
        required = ["chrom", "pos", "ref", "alt", "af_high", "af_low"]
        missing = [col for col in required if col not in header]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Get column indices
        col_idx = {col: i for i, col in enumerate(header)}

        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")

            try:
                af_high = float(fields[col_idx["af_high"]])
                af_low = float(fields[col_idx["af_low"]])

                variant = CandidateVariant(
                    chrom=fields[col_idx["chrom"]],
                    pos=int(fields[col_idx["pos"]]),
                    ref=fields[col_idx["ref"]],
                    alt=fields[col_idx["alt"]],
                    af_high=af_high,
                    af_low=af_low,
                    delta_af=float(fields[col_idx.get("delta_af", 0)])
                        if "delta_af" in col_idx else af_high - af_low,
                    dp_high=int(fields[col_idx.get("dp_high", 0)]) if "dp_high" in col_idx else 0,
                    dp_low=int(fields[col_idx.get("dp_low", 0)]) if "dp_low" in col_idx else 0,
                    region_rank=int(fields[col_idx.get("region_rank", 0)]) if "region_rank" in col_idx else 1,
                    distance_to_peak=int(fields[col_idx.get("distance_to_peak", 0)])
                        if "distance_to_peak" in col_idx else 0,
                )
                variants.append(variant)
            except (IndexError, ValueError) as e:
                logger.warning(f"Skipping invalid line {line_num}: {e}")
                continue

    logger.info(f"Read {len(variants):,} candidate variants from {path}")
    return variants
