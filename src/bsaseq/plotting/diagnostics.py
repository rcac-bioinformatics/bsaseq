"""Diagnostic plots for bsaseq.

This module provides functions for generating diagnostic plots
to assess data quality and analysis results.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend for headless environments

import matplotlib.pyplot as plt
import numpy as np

from bsaseq.plotting.style import (
    AF_CMAP,
    set_publication_style,
)
from bsaseq.utils.logging import get_logger

if TYPE_CHECKING:
    from collections.abc import Sequence

    from bsaseq.core.models import Variant

logger = get_logger(__name__)


def plot_af_distribution(
    variants: Sequence[Variant],
    output_path: str | Path | None = None,
    figsize: tuple[float, float] = (10, 8),
    dpi: int = 150,
) -> plt.Figure:
    """Plot allele frequency distributions for diagnostic purposes.

    Creates a 2x2 panel figure:
    - Top left: Histogram of af_high
    - Top right: Histogram of af_low
    - Bottom left: Histogram of delta_af
    - Bottom right: Scatter of af_high vs af_low

    Useful for QC and identifying issues with bulk composition.

    Args:
        variants: Sequence of Variant objects.
        output_path: If provided, save figure to this path.
        figsize: Figure dimensions (width, height) in inches.
        dpi: Resolution for PNG output.

    Returns:
        matplotlib Figure object.
    """
    set_publication_style()

    if not variants:
        logger.warning("No variants provided for AF distribution plot")
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        return fig

    logger.info("Generating AF distribution plot...")

    # Extract data
    af_high = np.array([v.af_high for v in variants])
    af_low = np.array([v.af_low for v in variants])
    delta_af = np.array([v.delta_af for v in variants])

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # Common histogram parameters
    bins = 50
    hist_kwargs = {
        "bins": bins,
        "edgecolor": "white",
        "linewidth": 0.5,
        "alpha": 0.8,
    }

    # Top left: af_high histogram
    ax = axes[0, 0]
    ax.hist(af_high, color="#d62728", **hist_kwargs)
    ax.set_xlabel("Allele Frequency (High Bulk)")
    ax.set_ylabel("Count")
    ax.set_xlim(0, 1)
    ax.set_title("High Bulk AF Distribution")
    _add_stats_text(ax, af_high, position="upper left")

    # Top right: af_low histogram
    ax = axes[0, 1]
    ax.hist(af_low, color="#1f77b4", **hist_kwargs)
    ax.set_xlabel("Allele Frequency (Low Bulk)")
    ax.set_ylabel("Count")
    ax.set_xlim(0, 1)
    ax.set_title("Low Bulk AF Distribution")
    _add_stats_text(ax, af_low, position="upper right")

    # Bottom left: delta_af histogram
    ax = axes[1, 0]
    ax.hist(delta_af, color="#2ca02c", **hist_kwargs)
    ax.axvline(x=0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("Delta AF (High - Low)")
    ax.set_ylabel("Count")
    ax.set_xlim(-1, 1)
    ax.set_title("Delta AF Distribution")
    _add_stats_text(ax, delta_af, position="upper left")

    # Bottom right: scatter of af_high vs af_low
    ax = axes[1, 1]

    # Use hexbin for dense data, scatter for sparse
    if len(af_high) > 5000:
        hb = ax.hexbin(
            af_low,
            af_high,
            gridsize=50,
            cmap="Blues",
            mincnt=1,
        )
        plt.colorbar(hb, ax=ax, shrink=0.8, label="Count")
    else:
        ax.scatter(
            af_low,
            af_high,
            c=delta_af,
            cmap=AF_CMAP,
            s=8,
            alpha=0.5,
            linewidths=0,
        )
        # Add colorbar for delta_af
        sm = plt.cm.ScalarMappable(
            cmap=AF_CMAP,
            norm=plt.Normalize(vmin=-1, vmax=1),
        )
        sm.set_array([])
        plt.colorbar(sm, ax=ax, shrink=0.8, label="Delta AF")

    # Add diagonal line (equal AF)
    ax.plot([0, 1], [0, 1], "k--", linewidth=1, alpha=0.5, label="Equal AF")
    ax.set_xlabel("AF (Low Bulk)")
    ax.set_ylabel("AF (High Bulk)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("AF High vs AF Low")
    ax.set_aspect("equal")

    # Overall title
    fig.suptitle(
        f"Allele Frequency Diagnostics (n={len(variants):,} variants)",
        fontsize=14,
        fontweight="bold",
    )

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        _save_figure(fig, output_path, dpi)

    return fig


def plot_depth_distribution(
    variants: Sequence[Variant],
    output_path: str | Path | None = None,
    figsize: tuple[float, float] = (10, 4),
    dpi: int = 150,
) -> plt.Figure:
    """Plot read depth distributions.

    Two-panel figure:
    - Left: Histogram of dp_high
    - Right: Histogram of dp_low

    Useful for checking coverage uniformity and filter appropriateness.

    Args:
        variants: Sequence of Variant objects.
        output_path: If provided, save figure to this path.
        figsize: Figure dimensions (width, height) in inches.
        dpi: Resolution for PNG output.

    Returns:
        matplotlib Figure object.
    """
    set_publication_style()

    if not variants:
        logger.warning("No variants provided for depth distribution plot")
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        return fig

    logger.info("Generating depth distribution plot...")

    # Extract data
    dp_high = np.array([v.dp_high for v in variants])
    dp_low = np.array([v.dp_low for v in variants])

    # Create figure
    fig, (ax_high, ax_low) = plt.subplots(1, 2, figsize=figsize)

    # Determine reasonable bin range
    max_dp = max(np.percentile(dp_high, 99), np.percentile(dp_low, 99))
    bins = np.linspace(0, max_dp * 1.1, 50)

    hist_kwargs = {
        "bins": bins,
        "edgecolor": "white",
        "linewidth": 0.5,
        "alpha": 0.8,
    }

    # Left: dp_high histogram
    ax_high.hist(dp_high, color="#d62728", **hist_kwargs)
    ax_high.set_xlabel("Read Depth")
    ax_high.set_ylabel("Count")
    ax_high.set_title("High Bulk Depth Distribution")
    _add_stats_text(ax_high, dp_high, position="upper right", decimals=0)

    # Right: dp_low histogram
    ax_low.hist(dp_low, color="#1f77b4", **hist_kwargs)
    ax_low.set_xlabel("Read Depth")
    ax_low.set_ylabel("Count")
    ax_low.set_title("Low Bulk Depth Distribution")
    _add_stats_text(ax_low, dp_low, position="upper right", decimals=0)

    # Overall title
    fig.suptitle(
        f"Read Depth Diagnostics (n={len(variants):,} variants)",
        fontsize=14,
        fontweight="bold",
    )

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        _save_figure(fig, output_path, dpi)

    return fig


def _add_stats_text(
    ax: plt.Axes,
    data: np.ndarray,
    position: str = "upper right",
    decimals: int = 3,
) -> None:
    """Add statistics text box to axes.

    Args:
        ax: matplotlib Axes object.
        data: Data array.
        position: Position for text box ('upper left', 'upper right').
        decimals: Number of decimal places for display.
    """
    mean_val = np.mean(data)
    median_val = np.median(data)
    std_val = np.std(data)

    if decimals == 0:
        text = f"Mean: {mean_val:.0f}\nMedian: {median_val:.0f}\nStd: {std_val:.0f}"
    else:
        text = (
            f"Mean: {mean_val:.{decimals}f}\n"
            f"Median: {median_val:.{decimals}f}\n"
            f"Std: {std_val:.{decimals}f}"
        )

    # Determine position
    if position == "upper left":
        x, y = 0.05, 0.95
        ha, va = "left", "top"
    else:  # upper right
        x, y = 0.95, 0.95
        ha, va = "right", "top"

    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment=va,
        horizontalalignment=ha,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )


def _save_figure(fig: plt.Figure, output_path: Path, dpi: int = 150) -> None:
    """Save figure in PNG format (and PDF for publication).

    Args:
        fig: matplotlib Figure to save.
        output_path: Base output path (without extension).
        dpi: Resolution for PNG output.
    """
    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Remove extension if present
    base_path = output_path.with_suffix("")

    # Save PNG
    png_path = base_path.with_suffix(".png")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    logger.info(f"Saved: {png_path}")
