"""Plotting modules for bsaseq.

This package provides publication-quality visualization functions
for BSA analysis results.
"""

from bsaseq.plotting.diagnostics import (
    plot_af_distribution,
    plot_depth_distribution,
)
from bsaseq.plotting.genome import (
    plot_all_regions,
    plot_genome_wide,
    plot_region,
)
from bsaseq.plotting.style import (
    AF_CMAP,
    CHROM_COLORS,
    SIGNIFICANT_COLOR,
    reset_style,
    set_publication_style,
)

__all__ = [
    # Genome plots
    "plot_genome_wide",
    "plot_region",
    "plot_all_regions",
    # Diagnostic plots
    "plot_af_distribution",
    "plot_depth_distribution",
    # Style
    "set_publication_style",
    "reset_style",
    "CHROM_COLORS",
    "SIGNIFICANT_COLOR",
    "AF_CMAP",
]
