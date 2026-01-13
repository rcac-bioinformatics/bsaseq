"""Matplotlib style configuration for bsaseq plots.

This module provides consistent styling for publication-quality figures.
"""

from __future__ import annotations

import matplotlib.pyplot as plt


def set_publication_style() -> None:
    """Configure matplotlib for publication-quality figures.

    Sets font family, sizes, and other visual parameters for
    clean, professional-looking plots suitable for publication.
    """
    plt.rcParams.update(
        {
            # Font settings
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
            "font.size": 10,
            # Title and label sizes
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.titlesize": 14,
            # Spine and line settings
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.0,
            "patch.linewidth": 0.5,
            # Grid settings
            "axes.grid": False,
            "grid.alpha": 0.3,
            # Figure settings
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "white",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.1,
        }
    )


def reset_style() -> None:
    """Reset matplotlib to default style."""
    plt.rcdefaults()


# Color palette for chromosomes (colorblind-friendly alternating colors)
CHROM_COLORS = ["#1f77b4", "#7f7f7f"]  # Blue / Gray

# Color for significant points
SIGNIFICANT_COLOR = "#d62728"  # Red

# Color for non-significant points
NONSIGNIFICANT_COLOR = "#7f7f7f"  # Gray

# Colormap for allele frequency
AF_CMAP = "RdYlBu_r"  # Red (high) - Yellow - Blue (low)

# Region highlight color
REGION_HIGHLIGHT_COLOR = "#ffcccc"  # Light red/pink
REGION_HIGHLIGHT_ALPHA = 0.3

# Peak marker color
PEAK_MARKER_COLOR = "#d62728"  # Red
PEAK_MARKER_ALPHA = 0.5

# Window line color
WINDOW_LINE_COLOR = "#2ca02c"  # Green
WINDOW_LINE_ALPHA = 0.7

# Threshold line color
THRESHOLD_LINE_COLOR = "#d62728"  # Red

# Highlighted variant marker
HIGHLIGHT_MARKER_COLOR = "#ff7f0e"  # Orange
HIGHLIGHT_MARKER_SIZE = 80
HIGHLIGHT_MARKER_EDGE = "black"

# Default marker size for scatter plots
DEFAULT_MARKER_SIZE = 8

# Small marker for dense plots
SMALL_MARKER_SIZE = 4
