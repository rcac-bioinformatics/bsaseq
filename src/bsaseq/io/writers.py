"""Output writers for bsaseq.

This module provides functions for writing analysis results to various
file formats, including TSV for variant, window, region, and candidate data.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from bsaseq.utils.logging import get_logger

if TYPE_CHECKING:
    from collections.abc import Iterable

    from bsaseq.analysis.candidates import (
        AnnotatedCandidate,
        CandidateGene,
        CandidateRegion,
        CandidateVariant,
    )
    from bsaseq.core.models import Variant, Window

logger = get_logger(__name__)


def write_variants_tsv(
    variants: Iterable[Variant],
    path: str | Path,
) -> int:
    """Write per-variant table to TSV file.

    Columns:
        chrom, pos, ref, alt, af_high, af_low, delta_af,
        dp_high, dp_low, n_samples_high, n_samples_low

    Args:
        variants: Iterable of Variant objects to write.
        path: Output file path.

    Returns:
        Number of variants written.

    Example:
        >>> variants = parse_vcf("data.vcf", "MUT", "WT")
        >>> n = write_variants_tsv(variants, "output_variants.tsv")
        >>> print(f"Wrote {n} variants")
    """
    path = Path(path)
    logger.info(f"Writing variants to {path}")

    count = 0
    with open(path, "w") as f:
        # Write header
        f.write(
            "chrom\tpos\tref\talt\taf_high\taf_low\tdelta_af\t"
            "dp_high\tdp_low\tn_samples_high\tn_samples_low\n"
        )

        # Write data rows
        for v in variants:
            f.write(
                f"{v.chrom}\t{v.pos}\t{v.ref}\t{v.alt}\t"
                f"{v.af_high:.6f}\t{v.af_low:.6f}\t{v.delta_af:.6f}\t"
                f"{v.dp_high}\t{v.dp_low}\t"
                f"{v.n_samples_high}\t{v.n_samples_low}\n"
            )
            count += 1

    logger.info(f"Wrote {count:,} variants to {path}")
    return count


def write_windows_tsv(
    windows: Iterable[Window],
    path: str | Path,
) -> int:
    """Write per-window table to TSV file.

    Columns:
        chrom, start, end, midpoint, n_variants, mean_delta_af,
        median_delta_af, tricube_delta_af, mean_af_high, mean_af_low,
        g_statistic, p_value, z_score

    Args:
        windows: Iterable of Window objects to write.
        path: Output file path.

    Returns:
        Number of windows written.

    Example:
        >>> windows = calculate_windows(variants)
        >>> n = write_windows_tsv(windows, "output_windows.tsv")
        >>> print(f"Wrote {n} windows")
    """
    path = Path(path)
    logger.info(f"Writing windows to {path}")

    count = 0
    with open(path, "w") as f:
        # Write header
        f.write(
            "chrom\tstart\tend\tmidpoint\tn_variants\t"
            "mean_delta_af\tmedian_delta_af\ttricube_delta_af\t"
            "mean_af_high\tmean_af_low\t"
            "g_statistic\tp_value\tz_score\n"
        )

        # Write data rows
        for w in windows:
            z_score_str = f"{w.z_score:.6f}" if w.z_score is not None else "NA"
            f.write(
                f"{w.chrom}\t{w.start}\t{w.end}\t{w.midpoint}\t{w.n_variants}\t"
                f"{w.mean_delta_af:.6f}\t{w.median_delta_af:.6f}\t"
                f"{w.tricube_delta_af:.6f}\t"
                f"{w.mean_af_high:.6f}\t{w.mean_af_low:.6f}\t"
                f"{w.g_statistic:.6f}\t{w.p_value:.6e}\t{z_score_str}\n"
            )
            count += 1

    logger.info(f"Wrote {count:,} windows to {path}")
    return count


def write_regions_tsv(
    regions: Iterable[CandidateRegion],
    path: str | Path,
) -> int:
    """Write candidate regions to TSV file.

    Columns:
        rank, chrom, start, end, length, n_windows, n_variants,
        max_z_score, mean_z_score, max_delta_af, mean_delta_af, peak_position

    Args:
        regions: Iterable of CandidateRegion objects to write.
        path: Output file path.

    Returns:
        Number of regions written.
    """
    path = Path(path)
    logger.info(f"Writing regions to {path}")

    count = 0
    with open(path, "w") as f:
        # Write header
        f.write(
            "rank\tchrom\tstart\tend\tlength\tn_windows\tn_variants\t"
            "max_z_score\tmean_z_score\tmax_delta_af\tmean_delta_af\t"
            "peak_position\n"
        )

        # Write data rows
        for rank, r in enumerate(regions, 1):
            f.write(
                f"{rank}\t{r.chrom}\t{r.start}\t{r.end}\t{r.length}\t"
                f"{r.n_windows}\t{r.n_variants}\t"
                f"{r.max_z_score:.6f}\t{r.mean_z_score:.6f}\t"
                f"{r.max_delta_af:.6f}\t{r.mean_delta_af:.6f}\t"
                f"{r.peak_position}\n"
            )
            count += 1

    logger.info(f"Wrote {count:,} regions to {path}")
    return count


def write_regions_bed(
    regions: Iterable[CandidateRegion],
    path: str | Path,
) -> int:
    """Write candidate regions to BED format for downstream tools.

    BED columns: chrom, start, end, name (candidate_N), score (max_z_score * 100)

    Note: BED format uses 0-based, half-open coordinates.
    Start is converted from 1-based to 0-based.

    Args:
        regions: Iterable of CandidateRegion objects to write.
        path: Output file path.

    Returns:
        Number of regions written.
    """
    path = Path(path)
    logger.info(f"Writing regions to BED format: {path}")

    count = 0
    with open(path, "w") as f:
        # BED has no header by default
        for rank, r in enumerate(regions, 1):
            # Convert to 0-based coordinates
            start_0 = r.start - 1
            # Score is typically 0-1000 in BED format
            score = min(1000, int(r.max_z_score * 100))
            f.write(
                f"{r.chrom}\t{start_0}\t{r.end}\t"
                f"candidate_{rank}\t{score}\n"
            )
            count += 1

    logger.info(f"Wrote {count:,} regions to {path}")
    return count


def write_candidate_variants_tsv(
    variants: Iterable[CandidateVariant],
    path: str | Path,
) -> int:
    """Write candidate variants to TSV file.

    Columns:
        chrom, pos, ref, alt, af_high, af_low, delta_af,
        dp_high, dp_low, region_rank, distance_to_peak

    Args:
        variants: Iterable of CandidateVariant objects to write.
        path: Output file path.

    Returns:
        Number of variants written.
    """
    path = Path(path)
    logger.info(f"Writing candidate variants to {path}")

    count = 0
    with open(path, "w") as f:
        # Write header
        f.write(
            "chrom\tpos\tref\talt\taf_high\taf_low\tdelta_af\t"
            "dp_high\tdp_low\tregion_rank\tdistance_to_peak\n"
        )

        # Write data rows
        for v in variants:
            f.write(
                f"{v.chrom}\t{v.pos}\t{v.ref}\t{v.alt}\t"
                f"{v.af_high:.6f}\t{v.af_low:.6f}\t{v.delta_af:.6f}\t"
                f"{v.dp_high}\t{v.dp_low}\t"
                f"{v.region_rank}\t{v.distance_to_peak}\n"
            )
            count += 1

    logger.info(f"Wrote {count:,} candidate variants to {path}")
    return count


def write_annotated_candidates_tsv(
    candidates: Iterable[AnnotatedCandidate],
    path: str | Path,
) -> int:
    """Write annotated candidate variants to TSV file.

    Columns:
        chrom, pos, ref, alt, af_high, af_low, delta_af,
        dp_high, dp_low, region_rank, distance_to_peak,
        effect, impact, gene_name, gene_id, hgvs_c, hgvs_p, is_lof

    Args:
        candidates: Iterable of AnnotatedCandidate objects to write.
        path: Output file path.

    Returns:
        Number of candidates written.
    """
    path = Path(path)
    logger.info(f"Writing annotated candidates to {path}")

    count = 0
    with open(path, "w") as f:
        # Write header
        f.write(
            "chrom\tpos\tref\talt\taf_high\taf_low\tdelta_af\t"
            "dp_high\tdp_low\tregion_rank\tdistance_to_peak\t"
            "effect\timpact\tgene_name\tgene_id\thgvs_c\thgvs_p\tis_lof\n"
        )

        # Write data rows
        for c in candidates:
            effect = c.effect or ""
            impact = c.impact or ""
            gene_name = c.gene_name or ""
            gene_id = c.gene_id or ""
            hgvs_c = c.hgvs_c or ""
            hgvs_p = c.hgvs_p or ""
            is_lof = "TRUE" if c.is_lof else "FALSE"

            f.write(
                f"{c.chrom}\t{c.pos}\t{c.ref}\t{c.alt}\t"
                f"{c.af_high:.6f}\t{c.af_low:.6f}\t{c.delta_af:.6f}\t"
                f"{c.dp_high}\t{c.dp_low}\t"
                f"{c.region_rank}\t{c.distance_to_peak}\t"
                f"{effect}\t{impact}\t{gene_name}\t{gene_id}\t"
                f"{hgvs_c}\t{hgvs_p}\t{is_lof}\n"
            )
            count += 1

    logger.info(f"Wrote {count:,} annotated candidates to {path}")
    return count


def write_candidate_genes_tsv(
    genes: Iterable[CandidateGene],
    path: str | Path,
) -> int:
    """Write candidate genes summary to TSV file.

    Columns:
        rank, gene_name, gene_id, n_variants, n_lof,
        n_high_impact, n_moderate_impact, best_region_rank

    Args:
        genes: Iterable of CandidateGene objects to write.
        path: Output file path.

    Returns:
        Number of genes written.
    """
    path = Path(path)
    logger.info(f"Writing candidate genes to {path}")

    count = 0
    with open(path, "w") as f:
        # Write header
        f.write(
            "rank\tgene_name\tgene_id\tn_variants\tn_lof\t"
            "n_high_impact\tn_moderate_impact\tbest_region_rank\n"
        )

        # Write data rows
        for rank, g in enumerate(genes, 1):
            f.write(
                f"{rank}\t{g.gene_name}\t{g.gene_id}\t"
                f"{g.n_variants}\t{g.n_lof}\t"
                f"{g.n_high_impact}\t{g.n_moderate_impact}\t"
                f"{g.best_region_rank}\n"
            )
            count += 1

    logger.info(f"Wrote {count:,} candidate genes to {path}")
    return count
