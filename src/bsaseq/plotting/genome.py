"""Genome-wide and regional plotting for bsaseq.

This module provides functions for generating Manhattan-style plots
and zoomed regional plots for BSA analysis results.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend for headless environments

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from bsaseq.plotting.style import (
    AF_CMAP,
    CHROM_COLORS,
    DEFAULT_MARKER_SIZE,
    HIGHLIGHT_MARKER_COLOR,
    HIGHLIGHT_MARKER_EDGE,
    HIGHLIGHT_MARKER_SIZE,
    PEAK_MARKER_ALPHA,
    PEAK_MARKER_COLOR,
    REGION_HIGHLIGHT_ALPHA,
    REGION_HIGHLIGHT_COLOR,
    SIGNIFICANT_COLOR,
    SMALL_MARKER_SIZE,
    THRESHOLD_LINE_COLOR,
    WINDOW_LINE_ALPHA,
    WINDOW_LINE_COLOR,
    set_publication_style,
)
from bsaseq.utils.logging import get_logger
from bsaseq.utils.sorting import (
    get_chromosome_order,
    simplify_chromosome_label,
    sort_chromosomes,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

    from bsaseq.analysis.candidates import CandidateRegion, CandidateVariant
    from bsaseq.core.models import Variant, Window

logger = get_logger(__name__)


def _calculate_genome_positions(
    chroms: Sequence[str],
    positions: Sequence[int],
    chrom_offsets: dict[str, int],
) -> np.ndarray:
    """Calculate genome-wide positions for plotting.

    Args:
        chroms: Chromosome names for each point.
        positions: Genomic positions for each point.
        chrom_offsets: Mapping of chromosome to x-offset.

    Returns:
        Array of genome-wide x positions.
    """
    return np.array(
        [chrom_offsets.get(c, 0) + p for c, p in zip(chroms, positions)]
    )


def _get_chromosome_offsets(
    windows: Sequence[Window],
    gap_fraction: float = 0.02,
) -> tuple[dict[str, int], dict[str, int], int]:
    """Calculate chromosome offsets for genome-wide plotting.

    Args:
        windows: Sequence of Window objects.
        gap_fraction: Fraction of total genome to use for gaps.

    Returns:
        Tuple of (chrom_offsets, chrom_lengths, total_length).
    """
    # Get chromosome lengths from windows
    chrom_lengths: dict[str, int] = {}
    for w in windows:
        if w.chrom not in chrom_lengths:
            chrom_lengths[w.chrom] = w.end
        else:
            chrom_lengths[w.chrom] = max(chrom_lengths[w.chrom], w.end)

    # Sort chromosomes naturally
    sorted_chroms = sort_chromosomes(list(chrom_lengths.keys()))

    # Calculate total length for gap sizing
    total_length = sum(chrom_lengths.values())
    gap_size = int(total_length * gap_fraction / max(1, len(sorted_chroms) - 1))

    # Calculate offsets
    chrom_offsets: dict[str, int] = {}
    current_offset = 0
    for chrom in sorted_chroms:
        chrom_offsets[chrom] = current_offset
        current_offset += chrom_lengths[chrom] + gap_size

    total_genome_length = current_offset - gap_size if sorted_chroms else 0

    return chrom_offsets, chrom_lengths, total_genome_length


def plot_genome_wide(
    windows: Sequence[Window],
    regions: Sequence[CandidateRegion] | None = None,
    z_threshold: float = 3.0,
    output_path: str | Path | None = None,
    figsize: tuple[float, float] = (14, 8),
    dpi: int = 150,
    title: str | None = None,
) -> plt.Figure:
    """Generate genome-wide Manhattan-style plot.

    Creates a two-panel figure:
    - Top panel: Delta AF (tricube-smoothed or mean) by genomic position
    - Bottom panel: Z-score by genomic position with threshold line

    Args:
        windows: Sequence of Window objects.
        regions: Optional candidate regions to highlight.
        z_threshold: Z-score threshold line to draw.
        output_path: If provided, save figure to this path (PNG and PDF).
        figsize: Figure dimensions (width, height) in inches.
        dpi: Resolution for PNG output.
        title: Optional figure title.

    Returns:
        matplotlib Figure object.
    """
    set_publication_style()

    if not windows:
        logger.warning("No windows provided for genome-wide plot")
        fig, axes = plt.subplots(2, 1, figsize=figsize)
        return fig

    logger.info("Generating genome-wide plot...")

    # Calculate chromosome offsets
    chrom_offsets, chrom_lengths, total_length = _get_chromosome_offsets(windows)
    sorted_chroms = sort_chromosomes(list(chrom_offsets.keys()))

    # Get chromosome order for coloring
    chrom_order = get_chromosome_order(sorted_chroms)

    # Extract data from windows
    chroms = [w.chrom for w in windows]
    midpoints = [w.midpoint for w in windows]
    delta_afs = [w.tricube_delta_af for w in windows]
    z_scores = [w.z_score if w.z_score is not None else 0.0 for w in windows]

    # Calculate genome-wide positions
    x_positions = _calculate_genome_positions(chroms, midpoints, chrom_offsets)

    # Assign colors based on chromosome
    colors = [CHROM_COLORS[chrom_order[c] % 2] for c in chroms]

    # Create figure
    fig, (ax_delta, ax_zscore) = plt.subplots(
        2, 1, figsize=figsize, sharex=True, height_ratios=[1, 1]
    )

    # === Top panel: Delta AF ===
    ax_delta.scatter(
        x_positions,
        delta_afs,
        c=colors,
        s=SMALL_MARKER_SIZE,
        alpha=0.7,
        linewidths=0,
    )
    ax_delta.axhline(y=0, color="gray", linestyle="-", linewidth=0.5, alpha=0.5)
    ax_delta.set_ylabel("Delta AF (smoothed)")
    ax_delta.set_ylim(-1.05, 1.05)

    # === Bottom panel: Z-score ===
    # Color points above threshold differently
    z_array = np.array(z_scores)
    significant_mask = np.abs(z_array) >= z_threshold

    # Non-significant points
    ax_zscore.scatter(
        x_positions[~significant_mask],
        z_array[~significant_mask],
        c=[colors[i] for i in range(len(colors)) if not significant_mask[i]],
        s=SMALL_MARKER_SIZE,
        alpha=0.7,
        linewidths=0,
    )

    # Significant points
    if np.any(significant_mask):
        ax_zscore.scatter(
            x_positions[significant_mask],
            z_array[significant_mask],
            c=SIGNIFICANT_COLOR,
            s=DEFAULT_MARKER_SIZE,
            alpha=0.9,
            linewidths=0,
            label=f"|Z| >= {z_threshold}",
        )

    # Threshold lines
    ax_zscore.axhline(
        y=z_threshold,
        color=THRESHOLD_LINE_COLOR,
        linestyle="--",
        linewidth=1,
        alpha=0.7,
    )
    ax_zscore.axhline(
        y=-z_threshold,
        color=THRESHOLD_LINE_COLOR,
        linestyle="--",
        linewidth=1,
        alpha=0.7,
    )
    ax_zscore.axhline(y=0, color="gray", linestyle="-", linewidth=0.5, alpha=0.5)

    ax_zscore.set_ylabel("Z-score")

    # Highlight regions if provided
    if regions:
        for region in regions:
            if region.chrom in chrom_offsets:
                region_start = chrom_offsets[region.chrom] + region.start
                region_end = chrom_offsets[region.chrom] + region.end
                # Shade on both panels
                for ax in [ax_delta, ax_zscore]:
                    ax.axvspan(
                        region_start,
                        region_end,
                        color=REGION_HIGHLIGHT_COLOR,
                        alpha=REGION_HIGHLIGHT_ALPHA,
                        zorder=0,
                    )

    # === X-axis configuration ===
    # Set ticks at chromosome midpoints
    tick_positions = []
    tick_labels = []
    for chrom in sorted_chroms:
        mid = chrom_offsets[chrom] + chrom_lengths[chrom] // 2
        tick_positions.append(mid)
        tick_labels.append(simplify_chromosome_label(chrom))

    ax_zscore.set_xticks(tick_positions)
    ax_zscore.set_xticklabels(tick_labels, rotation=0)
    ax_zscore.set_xlabel("Chromosome")
    ax_zscore.set_xlim(0, total_length)

    # Legend for significant points
    if np.any(significant_mask):
        ax_zscore.legend(loc="upper right", framealpha=0.9)

    # Title
    if title:
        fig.suptitle(title, fontsize=14, fontweight="bold")
    else:
        fig.suptitle("Genome-wide BSA Analysis", fontsize=14, fontweight="bold")

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        _save_figure(fig, output_path, dpi)

    return fig


def plot_region(
    variants: Sequence[Variant],
    region: CandidateRegion,
    windows: Sequence[Window] | None = None,
    highlight_variants: Sequence[Variant | CandidateVariant] | None = None,
    output_path: str | Path | None = None,
    figsize: tuple[float, float] = (12, 6),
    dpi: int = 150,
    padding_fraction: float = 0.1,
) -> plt.Figure:
    """Generate zoomed plot for a single candidate region.

    Shows individual variant delta_AF values with optional window overlay.

    Args:
        variants: All variants (will be filtered to region).
        region: CandidateRegion to plot.
        windows: Optional windows to overlay as line.
        highlight_variants: Variants to highlight (e.g., candidate causal variants).
        output_path: Save path (PNG and PDF).
        figsize: Figure dimensions.
        dpi: Resolution.
        padding_fraction: Fraction of region width to add as padding on each side.

    Returns:
        matplotlib Figure object.
    """
    set_publication_style()

    logger.info(
        f"Generating region plot for {region.chrom}:{region.start:,}-{region.end:,}"
    )

    # Filter variants to region with padding
    region_width = region.end - region.start
    padding = int(region_width * padding_fraction)
    plot_start = max(1, region.start - padding)
    plot_end = region.end + padding

    region_variants = [
        v
        for v in variants
        if v.chrom == region.chrom and plot_start <= v.pos <= plot_end
    ]

    if not region_variants:
        logger.warning(f"No variants found in region {region.chrom}:{region.start}-{region.end}")
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No variants in region", ha="center", va="center")
        return fig

    # Extract variant data
    positions = np.array([v.pos for v in region_variants])
    delta_afs = np.array([v.delta_af for v in region_variants])
    af_highs = np.array([v.af_high for v in region_variants])

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Shade region boundaries
    ax.axvspan(
        region.start,
        region.end,
        color=REGION_HIGHLIGHT_COLOR,
        alpha=REGION_HIGHLIGHT_ALPHA,
        zorder=0,
        label="Candidate region",
    )

    # Mark peak position
    ax.axvline(
        region.peak_position,
        color=PEAK_MARKER_COLOR,
        linestyle="--",
        linewidth=1.5,
        alpha=PEAK_MARKER_ALPHA,
        label="Peak position",
    )

    # Scatter plot of variants colored by af_high
    scatter = ax.scatter(
        positions,
        delta_afs,
        c=af_highs,
        cmap=AF_CMAP,
        s=DEFAULT_MARKER_SIZE * 2,
        alpha=0.8,
        linewidths=0.5,
        edgecolors="white",
        vmin=0,
        vmax=1,
        zorder=2,
    )

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("AF (high bulk)")

    # Overlay window smoothed line if provided
    if windows:
        region_windows = [
            w
            for w in windows
            if w.chrom == region.chrom and w.start <= plot_end and w.end >= plot_start
        ]
        if region_windows:
            # Sort by position
            region_windows.sort(key=lambda w: w.midpoint)
            win_positions = [w.midpoint for w in region_windows]
            win_delta_afs = [w.tricube_delta_af for w in region_windows]
            ax.plot(
                win_positions,
                win_delta_afs,
                color=WINDOW_LINE_COLOR,
                linewidth=2,
                alpha=WINDOW_LINE_ALPHA,
                label="Smoothed (tricube)",
                zorder=3,
            )

    # Highlight candidate variants if provided
    if highlight_variants:
        highlight_positions = []
        highlight_deltas = []
        for hv in highlight_variants:
            # Check if variant is in this region
            if hasattr(hv, "chrom") and hasattr(hv, "pos") and hasattr(hv, "delta_af"):
                if (
                    hv.chrom == region.chrom
                    and plot_start <= hv.pos <= plot_end
                ):
                    highlight_positions.append(hv.pos)
                    highlight_deltas.append(hv.delta_af)

        if highlight_positions:
            ax.scatter(
                highlight_positions,
                highlight_deltas,
                c=HIGHLIGHT_MARKER_COLOR,
                s=HIGHLIGHT_MARKER_SIZE,
                marker="o",
                edgecolors=HIGHLIGHT_MARKER_EDGE,
                linewidths=1.5,
                zorder=4,
                label="Candidate variants",
            )

    # Reference line at 0
    ax.axhline(y=0, color="gray", linestyle="-", linewidth=0.5, alpha=0.5)

    # Format x-axis as Mb
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda x, pos: f"{x / 1_000_000:.1f}")
    )
    ax.set_xlabel(f"Position on {region.chrom} (Mb)")
    ax.set_ylabel("Delta AF")
    ax.set_xlim(plot_start, plot_end)
    ax.set_ylim(-1.05, 1.05)

    # Title
    title = (
        f"{region.chrom}:{region.start:,}-{region.end:,} "
        f"(Z={region.max_z_score:.2f})"
    )
    ax.set_title(title, fontsize=12, fontweight="bold")

    # Legend
    ax.legend(loc="upper right", framealpha=0.9)

    plt.tight_layout()

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        _save_figure(fig, output_path, dpi)

    return fig


def plot_all_regions(
    variants: Sequence[Variant],
    regions: Sequence[CandidateRegion],
    windows: Sequence[Window],
    candidate_variants: Sequence[Variant | CandidateVariant],
    output_dir: str | Path,
    max_regions: int = 10,
    **kwargs,
) -> list[Path]:
    """Generate individual plots for top candidate regions.

    Args:
        variants: All variants.
        regions: All candidate regions (will use top max_regions).
        windows: All windows.
        candidate_variants: Filtered candidate variants to highlight.
        output_dir: Directory to save plots.
        max_regions: Maximum number of region plots to generate.
        **kwargs: Passed to plot_region().

    Returns:
        List of paths to generated plot files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated_paths: list[Path] = []

    # Limit to max_regions
    regions_to_plot = list(regions)[:max_regions]

    if not regions_to_plot:
        logger.info("No regions to plot")
        return generated_paths

    logger.info(f"Generating plots for {len(regions_to_plot)} region(s)...")

    for rank, region in enumerate(regions_to_plot, 1):
        # Create filename with rank and chromosome
        chrom_label = simplify_chromosome_label(region.chrom)
        if not chrom_label:
            chrom_label = region.chrom
        filename = f"region_{rank}_{region.chrom}"
        output_path = output_dir / filename

        try:
            fig = plot_region(
                variants=variants,
                region=region,
                windows=windows,
                highlight_variants=candidate_variants,
                output_path=output_path,
                **kwargs,
            )
            plt.close(fig)

            # Track generated files
            generated_paths.append(output_path.with_suffix(".png"))
            generated_paths.append(output_path.with_suffix(".pdf"))

        except Exception as e:
            logger.error(f"Failed to generate plot for region {rank}: {e}")
            continue

    logger.info(f"Generated {len(regions_to_plot)} region plot(s)")
    return generated_paths


def _save_figure(fig: plt.Figure, output_path: Path, dpi: int = 150) -> None:
    """Save figure in both PNG and PDF formats.

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

    # Save PDF
    pdf_path = base_path.with_suffix(".pdf")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    logger.info(f"Saved: {pdf_path}")
