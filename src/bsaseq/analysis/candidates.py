"""Candidate region detection and variant filtering for bsaseq.

This module provides functions for identifying genomic regions enriched
for high delta_AF signal and filtering variants that may be causal.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING

import numpy as np

from bsaseq.utils.logging import get_logger

if TYPE_CHECKING:
    from bsaseq.annotation.snpeff import VariantAnnotation
    from bsaseq.core.models import Variant, Window

logger = get_logger(__name__)


class InheritanceMode(str, Enum):
    """Expected inheritance pattern for the trait."""

    RECESSIVE = "recessive"  # AF_high ~ 1.0, AF_low ~ 0.0
    DOMINANT = "dominant"  # AF_high ~ 0.5+, AF_low ~ 0.0


@dataclass
class CandidateRegion:
    """A contiguous genomic region with elevated delta_AF signal.

    Attributes:
        chrom: Chromosome name.
        start: Start position of the region (1-based).
        end: End position of the region (1-based).
        n_windows: Number of significant windows in the region.
        n_variants: Number of variants in the region.
        max_z_score: Maximum Z-score among windows in the region.
        mean_z_score: Mean Z-score of windows in the region.
        max_delta_af: Maximum mean_delta_af among windows.
        mean_delta_af: Mean of mean_delta_af across windows.
        peak_position: Position of maximum signal within region.
    """

    chrom: str
    start: int
    end: int
    n_windows: int
    n_variants: int
    max_z_score: float
    mean_z_score: float
    max_delta_af: float
    mean_delta_af: float
    peak_position: int

    @property
    def length(self) -> int:
        """Return the length of the region in base pairs."""
        return self.end - self.start

    def contains(self, chrom: str, pos: int) -> bool:
        """Check if a position falls within this region.

        Args:
            chrom: Chromosome name.
            pos: Genomic position.

        Returns:
            True if position is within this region.
        """
        return self.chrom == chrom and self.start <= pos <= self.end

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"CandidateRegion({self.chrom}:{self.start}-{self.end}, "
            f"z={self.max_z_score:.2f}, n_windows={self.n_windows})"
        )


@dataclass
class CandidateVariant:
    """A variant within a candidate region that may be causal.

    Attributes:
        chrom: Chromosome name.
        pos: Genomic position.
        ref: Reference allele.
        alt: Alternate allele.
        af_high: Allele frequency in high bulk.
        af_low: Allele frequency in low bulk.
        delta_af: Difference in allele frequency (high - low).
        dp_high: Read depth in high bulk.
        dp_low: Read depth in low bulk.
        region_rank: Rank of containing region (1 = best).
        distance_to_peak: Distance from region peak position.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    af_high: float
    af_low: float
    delta_af: float
    dp_high: int
    dp_low: int
    region_rank: int
    distance_to_peak: int

    @classmethod
    def from_variant(
        cls,
        v: Variant,
        region: CandidateRegion,
        region_rank: int,
    ) -> CandidateVariant:
        """Create CandidateVariant from Variant and its containing region.

        Args:
            v: Source Variant object.
            region: CandidateRegion containing this variant.
            region_rank: Rank of the region (1 = best).

        Returns:
            New CandidateVariant instance.
        """
        return cls(
            chrom=v.chrom,
            pos=v.pos,
            ref=v.ref,
            alt=v.alt,
            af_high=v.af_high,
            af_low=v.af_low,
            delta_af=v.delta_af,
            dp_high=v.dp_high,
            dp_low=v.dp_low,
            region_rank=region_rank,
            distance_to_peak=abs(v.pos - region.peak_position),
        )

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"CandidateVariant({self.chrom}:{self.pos} {self.ref}>{self.alt}, "
            f"delta={self.delta_af:.3f}, rank={self.region_rank})"
        )


def get_variant_filter_thresholds(mode: InheritanceMode) -> dict[str, float]:
    """Return appropriate filtering thresholds for inheritance mode.

    Args:
        mode: Expected inheritance pattern.

    Returns:
        Dictionary with min_delta_af, min_af_high, max_af_low thresholds.
    """
    if mode == InheritanceMode.RECESSIVE:
        return {
            "min_delta_af": 0.8,
            "min_af_high": 0.9,
            "max_af_low": 0.1,
        }
    elif mode == InheritanceMode.DOMINANT:
        return {
            "min_delta_af": 0.3,
            "min_af_high": 0.4,
            "max_af_low": 0.1,
        }
    else:
        raise ValueError(f"Unknown inheritance mode: {mode}")


def _find_consecutive_runs(
    windows: list[Window],
    is_significant: list[bool],
    min_consecutive: int,
) -> list[list[Window]]:
    """Find runs of consecutive significant windows.

    Args:
        windows: List of windows for a single chromosome.
        is_significant: Boolean mask for significant windows.
        min_consecutive: Minimum consecutive windows for a run.

    Returns:
        List of runs, where each run is a list of consecutive significant windows.
    """
    runs = []
    current_run: list[Window] = []

    for window, sig in zip(windows, is_significant):
        if sig:
            current_run.append(window)
        else:
            if len(current_run) >= min_consecutive:
                runs.append(current_run)
            current_run = []

    # Don't forget the last run
    if len(current_run) >= min_consecutive:
        runs.append(current_run)

    return runs


def _create_region_from_windows(windows: list[Window]) -> CandidateRegion:
    """Create a CandidateRegion from a list of windows.

    Args:
        windows: List of consecutive significant windows.

    Returns:
        CandidateRegion summarizing the windows.
    """
    z_scores = [w.z_score for w in windows if w.z_score is not None]
    delta_afs = [w.mean_delta_af for w in windows]

    # Find peak window (maximum Z-score)
    max_z_idx = 0
    max_z = z_scores[0] if z_scores else 0
    for i, z in enumerate(z_scores):
        if z > max_z:
            max_z = z
            max_z_idx = i

    peak_window = windows[max_z_idx]

    return CandidateRegion(
        chrom=windows[0].chrom,
        start=windows[0].start,
        end=windows[-1].end,
        n_windows=len(windows),
        n_variants=sum(w.n_variants for w in windows),
        max_z_score=max(z_scores) if z_scores else 0.0,
        mean_z_score=float(np.mean(z_scores)) if z_scores else 0.0,
        max_delta_af=max(delta_afs),
        mean_delta_af=float(np.mean(delta_afs)),
        peak_position=peak_window.midpoint,
    )


def _merge_regions(
    regions: list[CandidateRegion],
    merge_distance: int,
) -> list[CandidateRegion]:
    """Merge nearby regions on the same chromosome.

    Args:
        regions: List of regions to potentially merge.
        merge_distance: Maximum gap between regions to merge.

    Returns:
        List of merged regions.
    """
    if not regions:
        return []

    # Sort by chromosome, then start position
    sorted_regions = sorted(regions, key=lambda r: (r.chrom, r.start))

    merged = []
    current = sorted_regions[0]

    for region in sorted_regions[1:]:
        # Check if regions should be merged
        if (
            region.chrom == current.chrom
            and region.start - current.end <= merge_distance
        ):
            # Merge the regions
            current = CandidateRegion(
                chrom=current.chrom,
                start=current.start,
                end=region.end,
                n_windows=current.n_windows + region.n_windows,
                n_variants=current.n_variants + region.n_variants,
                max_z_score=max(current.max_z_score, region.max_z_score),
                mean_z_score=(
                    (current.mean_z_score * current.n_windows +
                     region.mean_z_score * region.n_windows) /
                    (current.n_windows + region.n_windows)
                ),
                max_delta_af=max(current.max_delta_af, region.max_delta_af),
                mean_delta_af=(
                    (current.mean_delta_af * current.n_windows +
                     region.mean_delta_af * region.n_windows) /
                    (current.n_windows + region.n_windows)
                ),
                peak_position=(
                    current.peak_position
                    if current.max_z_score >= region.max_z_score
                    else region.peak_position
                ),
            )
        else:
            merged.append(current)
            current = region

    merged.append(current)
    return merged


def identify_candidate_regions(
    windows: list[Window],
    z_threshold: float = 3.0,
    min_consecutive: int = 2,
    merge_distance: int = 500_000,
) -> list[CandidateRegion]:
    """Identify candidate regions from sliding window analysis.

    A candidate region consists of consecutive windows exceeding the Z-score
    threshold. Adjacent regions within merge_distance are merged.

    Args:
        windows: List of Window objects (must be sorted by chrom, start).
        z_threshold: Minimum Z-score to consider window significant.
        min_consecutive: Minimum consecutive significant windows to form region.
        merge_distance: Merge regions within this distance (bp) on same chromosome.

    Returns:
        List of CandidateRegion objects, sorted by max_z_score descending.
    """
    if not windows:
        logger.warning("No windows provided for candidate region detection")
        return []

    logger.info(
        f"Identifying candidate regions: z_threshold={z_threshold}, "
        f"min_consecutive={min_consecutive}, merge_distance={merge_distance:,}"
    )

    # Group windows by chromosome
    chrom_windows: dict[str, list[Window]] = {}
    for window in windows:
        if window.chrom not in chrom_windows:
            chrom_windows[window.chrom] = []
        chrom_windows[window.chrom].append(window)

    # Find significant windows and create regions
    all_regions = []

    for _chrom, chrom_wins in chrom_windows.items():
        # Sort by start position
        chrom_wins.sort(key=lambda w: w.start)

        # Mark significant windows
        is_significant = [
            w.z_score is not None and w.z_score >= z_threshold
            for w in chrom_wins
        ]

        # Find consecutive runs
        runs = _find_consecutive_runs(chrom_wins, is_significant, min_consecutive)

        # Create regions from runs
        for run in runs:
            region = _create_region_from_windows(run)
            all_regions.append(region)

    if not all_regions:
        logger.warning(
            f"No candidate regions found with z_threshold={z_threshold}. "
            "Consider lowering the threshold."
        )
        return []

    # Merge nearby regions
    merged_regions = _merge_regions(all_regions, merge_distance)

    # Sort by max_z_score descending
    merged_regions.sort(key=lambda r: r.max_z_score, reverse=True)

    logger.info(f"Found {len(merged_regions)} candidate region(s)")
    for i, region in enumerate(merged_regions[:3], 1):
        logger.info(
            f"  Region {i}: {region.chrom}:{region.start:,}-{region.end:,} "
            f"(z={region.max_z_score:.2f})"
        )

    return merged_regions


def identify_candidate_regions_percentile(
    windows: list[Window],
    percentile: float = 99.0,
    min_consecutive: int = 2,
    merge_distance: int = 500_000,
) -> list[CandidateRegion]:
    """Identify candidates using percentile threshold instead of Z-score.

    Useful when delta_AF distribution is non-normal (e.g., multiple QTLs).

    Args:
        windows: List of Window objects.
        percentile: Windows above this percentile are considered significant.
        min_consecutive: Minimum consecutive significant windows.
        merge_distance: Merge distance in bp.

    Returns:
        List of CandidateRegion objects, sorted by max_z_score descending.
    """
    if not windows:
        logger.warning("No windows provided for candidate region detection")
        return []

    # Calculate percentile threshold based on mean_delta_af
    delta_afs = [w.mean_delta_af for w in windows]
    threshold = float(np.percentile(delta_afs, percentile))

    logger.info(
        f"Identifying candidate regions using {percentile}th percentile "
        f"(delta_af >= {threshold:.3f})"
    )

    # Group windows by chromosome
    chrom_windows: dict[str, list[Window]] = {}
    for window in windows:
        if window.chrom not in chrom_windows:
            chrom_windows[window.chrom] = []
        chrom_windows[window.chrom].append(window)

    # Find significant windows and create regions
    all_regions = []

    for _chrom, chrom_wins in chrom_windows.items():
        # Sort by start position
        chrom_wins.sort(key=lambda w: w.start)

        # Mark significant windows (using delta_af threshold)
        is_significant = [w.mean_delta_af >= threshold for w in chrom_wins]

        # Find consecutive runs
        runs = _find_consecutive_runs(chrom_wins, is_significant, min_consecutive)

        # Create regions from runs
        for run in runs:
            region = _create_region_from_windows(run)
            all_regions.append(region)

    if not all_regions:
        logger.warning(
            f"No candidate regions found at {percentile}th percentile. "
            "Consider lowering the percentile threshold."
        )
        return []

    # Merge nearby regions
    merged_regions = _merge_regions(all_regions, merge_distance)

    # Sort by max_z_score descending
    merged_regions.sort(key=lambda r: r.max_z_score, reverse=True)

    logger.info(f"Found {len(merged_regions)} candidate region(s)")

    return merged_regions


def filter_candidate_variants(
    variants: list[Variant],
    regions: list[CandidateRegion],
    min_delta_af: float = 0.8,
    min_af_high: float = 0.9,
    max_af_low: float = 0.1,
) -> list[CandidateVariant]:
    """Extract variants within candidate regions that meet causal variant criteria.

    For recessive mutations, expect:
    - AF near 1.0 in high bulk (homozygous mutant)
    - AF near 0.0 in low bulk (homozygous reference)
    - delta_AF near 1.0

    Args:
        variants: List of all Variant objects.
        regions: List of CandidateRegion objects (assumed sorted by max_z_score desc).
        min_delta_af: Minimum delta_AF threshold.
        min_af_high: Minimum AF in high bulk.
        max_af_low: Maximum AF in low bulk.

    Returns:
        List of CandidateVariant objects, sorted by region_rank then distance_to_peak.
    """
    if not regions:
        logger.info("No candidate regions provided for variant filtering")
        return []

    logger.info(
        f"Filtering candidate variants: min_delta_af={min_delta_af}, "
        f"min_af_high={min_af_high}, max_af_low={max_af_low}"
    )

    candidates = []

    for variant in variants:
        # Check if variant meets thresholds
        if variant.delta_af < min_delta_af:
            continue
        if variant.af_high < min_af_high:
            continue
        if variant.af_low > max_af_low:
            continue

        # Find which region contains this variant
        for rank, region in enumerate(regions, 1):
            if region.contains(variant.chrom, variant.pos):
                candidate = CandidateVariant.from_variant(variant, region, rank)
                candidates.append(candidate)
                break

    # Sort by region_rank, then distance_to_peak
    candidates.sort(key=lambda c: (c.region_rank, c.distance_to_peak))

    logger.info(f"Found {len(candidates)} candidate variant(s)")

    return candidates


def count_variants_in_regions(
    variants: list[Variant],
    regions: list[CandidateRegion],
) -> dict[int, int]:
    """Count variants falling within each region.

    Args:
        variants: List of Variant objects.
        regions: List of CandidateRegion objects.

    Returns:
        Dictionary mapping region index to variant count.
    """
    counts: dict[int, int] = dict.fromkeys(range(len(regions)), 0)

    for variant in variants:
        for i, region in enumerate(regions):
            if region.contains(variant.chrom, variant.pos):
                counts[i] += 1
                break

    return counts


@dataclass
class AnnotatedCandidate:
    """Candidate variant with functional annotation from snpEff.

    Attributes:
        chrom: Chromosome name.
        pos: Genomic position.
        ref: Reference allele.
        alt: Alternate allele.
        af_high: Allele frequency in high bulk.
        af_low: Allele frequency in low bulk.
        delta_af: Difference in allele frequency (high - low).
        dp_high: Read depth in high bulk.
        dp_low: Read depth in low bulk.
        region_rank: Rank of containing region (1 = best).
        distance_to_peak: Distance from region peak position.
        effect: Variant effect type (e.g., "missense_variant").
        impact: Impact severity (HIGH, MODERATE, LOW, MODIFIER).
        gene_name: Gene symbol.
        gene_id: Gene identifier.
        hgvs_c: Coding DNA change notation.
        hgvs_p: Protein change notation.
        is_lof: Whether variant is loss-of-function.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    af_high: float
    af_low: float
    delta_af: float
    dp_high: int
    dp_low: int
    region_rank: int
    distance_to_peak: int
    effect: str | None = None
    impact: str | None = None
    gene_name: str | None = None
    gene_id: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    is_lof: bool = False

    @classmethod
    def from_candidate(
        cls,
        candidate: CandidateVariant,
        annotation: VariantAnnotation | None = None,
    ) -> AnnotatedCandidate:
        """Create AnnotatedCandidate from CandidateVariant and optional annotation.

        Args:
            candidate: Source CandidateVariant object.
            annotation: Optional VariantAnnotation from snpEff.

        Returns:
            New AnnotatedCandidate instance.
        """
        if annotation is not None:
            return cls(
                chrom=candidate.chrom,
                pos=candidate.pos,
                ref=candidate.ref,
                alt=candidate.alt,
                af_high=candidate.af_high,
                af_low=candidate.af_low,
                delta_af=candidate.delta_af,
                dp_high=candidate.dp_high,
                dp_low=candidate.dp_low,
                region_rank=candidate.region_rank,
                distance_to_peak=candidate.distance_to_peak,
                effect=annotation.effect,
                impact=annotation.impact.value,
                gene_name=annotation.gene_name,
                gene_id=annotation.gene_id,
                hgvs_c=annotation.hgvs_c,
                hgvs_p=annotation.hgvs_p,
                is_lof=annotation.is_lof,
            )
        else:
            return cls(
                chrom=candidate.chrom,
                pos=candidate.pos,
                ref=candidate.ref,
                alt=candidate.alt,
                af_high=candidate.af_high,
                af_low=candidate.af_low,
                delta_af=candidate.delta_af,
                dp_high=candidate.dp_high,
                dp_low=candidate.dp_low,
                region_rank=candidate.region_rank,
                distance_to_peak=candidate.distance_to_peak,
            )

    def __repr__(self) -> str:
        """Return string representation."""
        gene = self.gene_name or "unknown"
        effect = self.effect or "unannotated"
        return (
            f"AnnotatedCandidate({self.chrom}:{self.pos} {self.ref}>{self.alt}, "
            f"{gene}, {effect}, rank={self.region_rank})"
        )


@dataclass
class CandidateGene:
    """Summary of a gene containing candidate variants.

    Attributes:
        gene_name: Gene symbol.
        gene_id: Gene identifier.
        n_variants: Number of candidate variants in gene.
        n_lof: Number of loss-of-function variants.
        n_high_impact: Number of HIGH impact variants.
        n_moderate_impact: Number of MODERATE impact variants.
        best_region_rank: Best (lowest) region rank among variants.
        variants: List of AnnotatedCandidate in this gene.
    """

    gene_name: str
    gene_id: str
    n_variants: int
    n_lof: int
    n_high_impact: int
    n_moderate_impact: int
    best_region_rank: int
    variants: list[AnnotatedCandidate] = field(default_factory=list)

    @property
    def priority_score(self) -> int:
        """Calculate priority score for ranking genes.

        Lower score = higher priority.
        """
        # Prioritize: region rank, LOF count, high impact, moderate impact
        return (
            self.best_region_rank * 1000
            - self.n_lof * 100
            - self.n_high_impact * 10
            - self.n_moderate_impact
        )

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"CandidateGene({self.gene_name}, n_variants={self.n_variants}, "
            f"n_lof={self.n_lof}, rank={self.best_region_rank})"
        )


def annotate_candidates(
    candidates: list[CandidateVariant],
    annotations: dict[tuple[str, int, str, str], VariantAnnotation],
) -> list[AnnotatedCandidate]:
    """Merge candidate variants with snpEff annotations.

    Args:
        candidates: List of CandidateVariant objects.
        annotations: Dictionary mapping (chrom, pos, ref, alt) to best annotation.

    Returns:
        List of AnnotatedCandidate objects.
    """
    annotated = []
    matched = 0
    unmatched = 0

    for candidate in candidates:
        key = (candidate.chrom, candidate.pos, candidate.ref, candidate.alt)
        annotation = annotations.get(key)
        if annotation is not None:
            matched += 1
        else:
            unmatched += 1
        annotated.append(AnnotatedCandidate.from_candidate(candidate, annotation))

    logger.info(f"Annotated {matched} candidates, {unmatched} without annotation")
    return annotated


def summarize_candidate_genes(
    annotated: list[AnnotatedCandidate],
) -> list[CandidateGene]:
    """Summarize candidates by gene.

    Args:
        annotated: List of AnnotatedCandidate objects.

    Returns:
        List of CandidateGene objects, sorted by priority score.
    """
    # Group by gene
    gene_variants: dict[str, list[AnnotatedCandidate]] = {}
    gene_ids: dict[str, str] = {}

    for candidate in annotated:
        if candidate.gene_name:
            gene = candidate.gene_name
            if gene not in gene_variants:
                gene_variants[gene] = []
                gene_ids[gene] = candidate.gene_id or ""
            gene_variants[gene].append(candidate)

    # Create gene summaries
    genes = []
    for gene_name, variants in gene_variants.items():
        n_lof = sum(1 for v in variants if v.is_lof)
        n_high = sum(1 for v in variants if v.impact == "HIGH")
        n_moderate = sum(1 for v in variants if v.impact == "MODERATE")
        best_rank = min(v.region_rank for v in variants)

        genes.append(
            CandidateGene(
                gene_name=gene_name,
                gene_id=gene_ids[gene_name],
                n_variants=len(variants),
                n_lof=n_lof,
                n_high_impact=n_high,
                n_moderate_impact=n_moderate,
                best_region_rank=best_rank,
                variants=variants,
            )
        )

    # Sort by priority score
    genes.sort(key=lambda g: g.priority_score)

    logger.info(f"Summarized {len(genes)} candidate genes")
    if genes:
        top_genes = genes[:5]
        for i, gene in enumerate(top_genes, 1):
            logger.info(
                f"  Gene {i}: {gene.gene_name} "
                f"(n={gene.n_variants}, lof={gene.n_lof}, rank={gene.best_region_rank})"
            )

    return genes
