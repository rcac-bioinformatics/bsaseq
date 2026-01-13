"""Natural sorting utilities for chromosome names.

This module provides functions for sorting chromosome names in natural order,
where Chr1 < Chr2 < Chr10 < ChrX rather than lexicographic order.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence


def natural_sort_key(chrom: str) -> tuple:
    """Generate sort key for natural chromosome ordering.

    Splits the chromosome name into text and numeric parts, converting
    numeric parts to integers for proper numeric comparison.

    Args:
        chrom: Chromosome name (e.g., "Chr1", "chr10", "X").

    Returns:
        Tuple suitable for sorting comparison.

    Examples:
        >>> natural_sort_key("Chr1") < natural_sort_key("Chr2")
        True
        >>> natural_sort_key("Chr2") < natural_sort_key("Chr10")
        True
        >>> natural_sort_key("Chr10") < natural_sort_key("ChrX")
        True
    """
    parts = re.split(r"(\d+)", chrom)
    return tuple(int(p) if p.isdigit() else p.lower() for p in parts)


def sort_chromosomes(chroms: Sequence[str]) -> list[str]:
    """Return chromosomes in natural sort order.

    Args:
        chroms: Sequence of chromosome names.

    Returns:
        List of chromosome names sorted in natural order.

    Examples:
        >>> sort_chromosomes(['Chr10', 'Chr2', 'Chr1', 'ChrX', 'Chr3'])
        ['Chr1', 'Chr2', 'Chr3', 'Chr10', 'ChrX']
        >>> sort_chromosomes(['chr1', 'chr10', 'chr2'])
        ['chr1', 'chr2', 'chr10']
    """
    return sorted(chroms, key=natural_sort_key)


def get_chromosome_order(chroms: Sequence[str]) -> dict[str, int]:
    """Get mapping of chromosome names to their sort order.

    Args:
        chroms: Sequence of chromosome names.

    Returns:
        Dictionary mapping chromosome name to sort order index.

    Examples:
        >>> order = get_chromosome_order(['Chr10', 'Chr2', 'Chr1'])
        >>> order['Chr1']
        0
        >>> order['Chr2']
        1
        >>> order['Chr10']
        2
    """
    sorted_chroms = sort_chromosomes(chroms)
    return {chrom: i for i, chrom in enumerate(sorted_chroms)}


def simplify_chromosome_label(chrom: str) -> str:
    """Simplify chromosome label for display.

    Strips common prefixes like "chr", "Chr", "chromosome" to get
    cleaner axis labels.

    Args:
        chrom: Full chromosome name.

    Returns:
        Simplified label for display.

    Examples:
        >>> simplify_chromosome_label("Chr01")
        '1'
        >>> simplify_chromosome_label("chromosome10")
        '10'
        >>> simplify_chromosome_label("X")
        'X'
    """
    # Remove common prefixes (case insensitive)
    label = re.sub(r"^(chromosome|chrom|chr)", "", chrom, flags=re.IGNORECASE)

    # Remove leading zeros from numeric parts
    if label.isdigit():
        label = str(int(label))

    return label if label else chrom
