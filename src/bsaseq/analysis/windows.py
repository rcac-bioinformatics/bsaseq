"""Sliding window analysis for bsaseq.

This module provides functions for calculating sliding window statistics
across the genome, including tricube-smoothed delta AF, G-statistics,
and Z-scores.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as scipy_stats

from bsaseq.core.models import Variant, Window
from bsaseq.utils.logging import get_logger

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

logger = get_logger(__name__)


def tricube_weight(distance: float, max_distance: float) -> float:
    """Calculate tricube kernel weight.

    The tricube weight function is commonly used in LOESS smoothing:
        w(d) = (1 - (d/D)^3)^3 for d < D
        w(d) = 0 for d >= D

    Args:
        distance: Distance from center point.
        max_distance: Maximum distance (bandwidth).

    Returns:
        Tricube weight value (0.0 to 1.0).

    Examples:
        >>> tricube_weight(0, 100)
        1.0
        >>> tricube_weight(100, 100)
        0.0
        >>> tricube_weight(50, 100)  # ~0.67
        0.669921875
    """
    if max_distance <= 0 or distance >= max_distance:
        return 0.0
    u = distance / max_distance
    return (1 - u**3) ** 3


def calculate_g_statistic(
    obs_high_alt: int,
    obs_high_ref: int,
    obs_low_alt: int,
    obs_low_ref: int,
) -> float:
    """Calculate G-statistic for independence test.

    The G-test is a likelihood ratio test for independence between
    allele frequencies in the two bulks. The statistic is approximately
    chi-squared distributed with df=1.

    Args:
        obs_high_alt: Observed alt allele count in high bulk.
        obs_high_ref: Observed ref allele count in high bulk.
        obs_low_alt: Observed alt allele count in low bulk.
        obs_low_ref: Observed ref allele count in low bulk.

    Returns:
        G-statistic value.
    """
    total = obs_high_alt + obs_high_ref + obs_low_alt + obs_low_ref
    if total == 0:
        return 0.0

    # Row and column totals
    row_high = obs_high_alt + obs_high_ref
    row_low = obs_low_alt + obs_low_ref
    col_alt = obs_high_alt + obs_low_alt
    col_ref = obs_high_ref + obs_low_ref

    # Expected values under independence
    eps = 1e-10
    if row_high == 0 or row_low == 0 or col_alt == 0 or col_ref == 0:
        return 0.0

    exp_high_alt = (row_high * col_alt) / total
    exp_high_ref = (row_high * col_ref) / total
    exp_low_alt = (row_low * col_alt) / total
    exp_low_ref = (row_low * col_ref) / total

    # G-statistic: 2 * sum(O * ln(O/E))
    g = 0.0
    for obs, exp in [
        (obs_high_alt, exp_high_alt),
        (obs_high_ref, exp_high_ref),
        (obs_low_alt, exp_low_alt),
        (obs_low_ref, exp_low_ref),
    ]:
        if obs > 0 and exp > eps:
            g += obs * math.log(obs / exp)

    return 2 * g


def g_statistic_to_pvalue(g_stat: float) -> float:
    """Convert G-statistic to p-value using chi-squared distribution.

    Args:
        g_stat: G-statistic value.

    Returns:
        P-value (two-tailed, df=1).
    """
    if g_stat <= 0:
        return 1.0
    return float(scipy_stats.chi2.sf(g_stat, df=1))


def _calculate_window_stats(
    variants: list[Variant],
    chrom: str,
    start: int,
    end: int,
) -> Window | None:
    """Calculate statistics for a single window.

    Args:
        variants: List of variants within the window.
        chrom: Chromosome name.
        start: Window start position.
        end: Window end position.

    Returns:
        Window object with computed statistics, or None if no variants.
    """
    if not variants:
        return None

    n_variants = len(variants)

    # Extract delta_af values
    delta_afs = np.array([v.delta_af for v in variants])
    af_highs = np.array([v.af_high for v in variants])
    af_lows = np.array([v.af_low for v in variants])

    # Basic statistics
    mean_delta_af = float(np.mean(delta_afs))
    median_delta_af = float(np.median(delta_afs))
    mean_af_high = float(np.mean(af_highs))
    mean_af_low = float(np.mean(af_lows))

    # Tricube-smoothed delta AF
    midpoint = (start + end) // 2
    half_window = (end - start) // 2

    positions = np.array([v.pos for v in variants])
    distances = np.abs(positions - midpoint)

    # Calculate tricube weights
    weights = np.array([tricube_weight(d, half_window) for d in distances])
    weight_sum = np.sum(weights)

    if weight_sum > 0:
        tricube_delta_af = float(np.sum(delta_afs * weights) / weight_sum)
    else:
        tricube_delta_af = mean_delta_af

    # Pool allele counts for G-statistic
    total_high_alt = sum(v.alt_count_high for v in variants)
    total_high_ref = sum(v.ref_count_high for v in variants)
    total_low_alt = sum(v.alt_count_low for v in variants)
    total_low_ref = sum(v.ref_count_low for v in variants)

    g_stat = calculate_g_statistic(
        total_high_alt, total_high_ref, total_low_alt, total_low_ref
    )
    p_value = g_statistic_to_pvalue(g_stat)

    return Window(
        chrom=chrom,
        start=start,
        end=end,
        n_variants=n_variants,
        mean_delta_af=mean_delta_af,
        median_delta_af=median_delta_af,
        mean_af_high=mean_af_high,
        mean_af_low=mean_af_low,
        tricube_delta_af=tricube_delta_af,
        g_statistic=g_stat,
        p_value=p_value,
        z_score=None,  # Added in post-processing
    )


def _process_chromosome(
    variants: list[Variant],
    window_size: int,
    step_size: int,
    min_variants: int,
) -> list[Window]:
    """Process all windows for a single chromosome.

    Args:
        variants: List of variants for the chromosome (sorted by position).
        window_size: Window width in base pairs.
        step_size: Step between window starts in base pairs.
        min_variants: Minimum variants required per window.

    Returns:
        List of Window objects for the chromosome.
    """
    if not variants:
        return []

    chrom = variants[0].chrom
    positions = [v.pos for v in variants]

    # Determine chromosome extent from variant positions
    min_pos = min(positions)
    max_pos = max(positions)

    windows = []

    # Generate windows starting from just before the first variant
    window_start = max(1, min_pos - window_size // 2)

    while window_start <= max_pos:
        window_end = window_start + window_size - 1

        # Find variants within this window using binary search
        # (variants are already sorted by position)
        window_variants = [
            v for v in variants if window_start <= v.pos <= window_end
        ]

        if len(window_variants) >= min_variants:
            window = _calculate_window_stats(
                window_variants, chrom, window_start, window_end
            )
            if window is not None:
                windows.append(window)

        window_start += step_size

    return windows


def _add_z_scores(windows: list[Window]) -> list[Window]:
    """Add Z-scores based on genome-wide distribution of mean delta AF.

    Z-scores are calculated as:
        z = (window_mean_delta_af - genome_mean) / genome_sd

    Args:
        windows: List of Window objects.

    Returns:
        Same list with z_score field populated.
    """
    if not windows:
        return windows

    mean_values = np.array([w.mean_delta_af for w in windows])
    global_mean = float(np.mean(mean_values))
    global_sd = float(np.std(mean_values))

    if global_sd == 0 or np.isnan(global_sd):
        global_sd = 1e-10  # Avoid division by zero

    for window in windows:
        window.z_score = (window.mean_delta_af - global_mean) / global_sd

    return windows


def calculate_windows(
    variants: Iterable[Variant],
    window_size: int = 1_000_000,
    step_size: int = 250_000,
    min_variants: int = 5,
) -> Iterator[Window]:
    """Calculate sliding window statistics across the genome.

    Processes variants chromosome by chromosome, calculating various
    statistics for each window including mean/median delta AF,
    tricube-smoothed values, G-statistics, and p-values. Z-scores
    are computed in a post-processing step using genome-wide statistics.

    Args:
        variants: Iterable of Variant objects (must be sorted by chrom, pos).
        window_size: Window width in base pairs (default: 1 Mb).
        step_size: Step between window starts in base pairs (default: 250 kb).
        min_variants: Minimum variants required per window (default: 5).
            Windows with fewer variants are skipped.

    Yields:
        Window objects with computed statistics.

    Note:
        Variants are consumed into memory per-chromosome for efficiency.
        Memory usage is O(variants per chromosome), not O(total variants).

    Example:
        >>> variants = parse_vcf("data.vcf", "MUT", "WT")
        >>> for window in calculate_windows(variants, window_size=500000):
        ...     print(f"{window.chrom}:{window.start}-{window.end} "
        ...           f"delta={window.mean_delta_af:.3f}")
    """
    logger.info(
        f"Starting window analysis: size={window_size:,} bp, "
        f"step={step_size:,} bp, min_variants={min_variants}"
    )

    current_chrom: str | None = None
    chrom_variants: list[Variant] = []
    all_windows: list[Window] = []

    for variant in variants:
        if current_chrom is None:
            current_chrom = variant.chrom

        if variant.chrom != current_chrom:
            # Process completed chromosome
            logger.info(
                f"Processing chromosome {current_chrom}: "
                f"{len(chrom_variants):,} variants"
            )
            windows = _process_chromosome(
                chrom_variants, window_size, step_size, min_variants
            )
            all_windows.extend(windows)
            logger.info(f"  Generated {len(windows):,} windows")

            # Start new chromosome
            current_chrom = variant.chrom
            chrom_variants = []

        chrom_variants.append(variant)

    # Process final chromosome
    if chrom_variants:
        logger.info(
            f"Processing chromosome {current_chrom}: "
            f"{len(chrom_variants):,} variants"
        )
        windows = _process_chromosome(
            chrom_variants, window_size, step_size, min_variants
        )
        all_windows.extend(windows)
        logger.info(f"  Generated {len(windows):,} windows")

    # Compute Z-scores (requires all windows)
    logger.info(f"Computing Z-scores for {len(all_windows):,} total windows")
    all_windows = _add_z_scores(all_windows)

    yield from all_windows


def calculate_windows_to_list(
    variants: Iterable[Variant],
    window_size: int = 1_000_000,
    step_size: int = 250_000,
    min_variants: int = 5,
) -> list[Window]:
    """Calculate sliding window statistics and return as list.

    Convenience wrapper around calculate_windows() that returns a list
    instead of an iterator.

    Args:
        variants: Iterable of Variant objects (must be sorted by chrom, pos).
        window_size: Window width in base pairs.
        step_size: Step between window starts in base pairs.
        min_variants: Minimum variants required per window.

    Returns:
        List of Window objects with computed statistics.
    """
    return list(
        calculate_windows(
            variants=variants,
            window_size=window_size,
            step_size=step_size,
            min_variants=min_variants,
        )
    )
