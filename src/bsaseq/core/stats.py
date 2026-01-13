"""Statistical calculations for bsaseq.

This module provides statistical functions for BSA analysis,
including allele frequency calculations and significance testing.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

    from bsaseq.core.models import Variant


def calculate_allele_frequency(ref_depth: int, alt_depth: int) -> float:
    """Calculate allele frequency from read depths.

    Args:
        ref_depth: Number of reads supporting reference allele.
        alt_depth: Number of reads supporting alternate allele.

    Returns:
        Frequency of the alternate allele (0.0 to 1.0).
        Returns 0.0 if total depth is 0.
    """
    total = ref_depth + alt_depth
    if total == 0:
        return 0.0
    return alt_depth / total


def calculate_delta_af(variants: Sequence[Variant]) -> np.ndarray:
    """Calculate delta allele frequencies for a sequence of variants.

    Args:
        variants: Sequence of Variant objects.

    Returns:
        NumPy array of delta AF values.
    """
    return np.array([v.delta_af for v in variants], dtype=np.float64)


def calculate_mean_delta_af(variants: Sequence[Variant]) -> float:
    """Calculate mean delta AF for a set of variants.

    Args:
        variants: Sequence of Variant objects.

    Returns:
        Mean delta AF value, or 0.0 if no variants.
    """
    if not variants:
        return 0.0
    delta_afs = calculate_delta_af(variants)
    return float(np.mean(delta_afs))


def calculate_median_delta_af(variants: Sequence[Variant]) -> float:
    """Calculate median delta AF for a set of variants.

    Args:
        variants: Sequence of Variant objects.

    Returns:
        Median delta AF value, or 0.0 if no variants.
    """
    if not variants:
        return 0.0
    delta_afs = calculate_delta_af(variants)
    return float(np.median(delta_afs))


def tricube_weight(distance: float, max_distance: float) -> float:
    """Calculate tricube weight for LOESS smoothing.

    The tricube weight function is:
        w(d) = (1 - (d/D)^3)^3 for d < D
        w(d) = 0 for d >= D

    Args:
        distance: Distance from center point.
        max_distance: Maximum distance (bandwidth).

    Returns:
        Tricube weight value (0.0 to 1.0).
    """
    if max_distance <= 0 or distance >= max_distance:
        return 0.0
    u = distance / max_distance
    return (1 - u**3) ** 3


def weighted_mean(
    values: np.ndarray,
    weights: np.ndarray,
) -> float:
    """Calculate weighted mean.

    Args:
        values: Array of values.
        weights: Array of weights (must be same length as values).

    Returns:
        Weighted mean, or 0.0 if sum of weights is 0.
    """
    weight_sum = np.sum(weights)
    if weight_sum == 0:
        return 0.0
    return float(np.sum(values * weights) / weight_sum)


def calculate_g_statistic(
    af_high: float,
    af_low: float,
    dp_high: int,
    dp_low: int,
) -> float:
    """Calculate G-statistic for allele frequency divergence.

    The G-statistic measures the divergence between observed and
    expected allele frequencies under the null hypothesis of
    equal frequencies in both bulks.

    Args:
        af_high: Allele frequency in high bulk.
        af_low: Allele frequency in low bulk.
        dp_high: Read depth in high bulk.
        dp_low: Read depth in low bulk.

    Returns:
        G-statistic value.
    """
    total_dp = dp_high + dp_low
    if total_dp == 0:
        return 0.0

    # Observed counts
    obs_high_alt = af_high * dp_high
    obs_high_ref = (1 - af_high) * dp_high
    obs_low_alt = af_low * dp_low
    obs_low_ref = (1 - af_low) * dp_low

    # Expected under null
    pooled_af = (obs_high_alt + obs_low_alt) / total_dp
    exp_high_alt = pooled_af * dp_high
    exp_high_ref = (1 - pooled_af) * dp_high
    exp_low_alt = pooled_af * dp_low
    exp_low_ref = (1 - pooled_af) * dp_low

    # G-statistic
    eps = 1e-10
    g = 0.0
    for obs, exp in [
        (obs_high_alt, exp_high_alt),
        (obs_high_ref, exp_high_ref),
        (obs_low_alt, exp_low_alt),
        (obs_low_ref, exp_low_ref),
    ]:
        if obs > 0 and exp > eps:
            g += obs * np.log(obs / exp)

    return float(2 * g)
