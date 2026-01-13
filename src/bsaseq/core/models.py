"""Data models for bsaseq.

This module defines the core data structures used throughout the package
for representing variants and genomic windows.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Variant:
    """Single variant with allele frequency data from both bulks.

    Attributes:
        chrom: Chromosome name.
        pos: 1-based genomic position.
        ref: Reference allele.
        alt: Alternate allele.
        af_high: Allele frequency in high-bulk (e.g., mutant pool).
        af_low: Allele frequency in low-bulk (e.g., wild-type pool).
        dp_high: Read depth in high-bulk (pooled across samples).
        dp_low: Read depth in low-bulk (pooled across samples).
        n_samples_high: Number of samples in high bulk.
        n_samples_low: Number of samples in low bulk.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    af_high: float
    af_low: float
    dp_high: int
    dp_low: int
    n_samples_high: int = 1
    n_samples_low: int = 1

    @property
    def delta_af(self) -> float:
        """Allele frequency difference (high - low).

        A positive value indicates higher alternate allele frequency
        in the high bulk (mutant) compared to low bulk (wild-type).

        Returns:
            Difference in allele frequencies between bulks.
        """
        return self.af_high - self.af_low

    @property
    def snp_index_high(self) -> float:
        """SNP-index for high bulk.

        The SNP-index represents the proportion of reads supporting
        the alternate allele. For the alternate allele, this equals
        the allele frequency.

        Returns:
            SNP-index value for high bulk.
        """
        return self.af_high

    @property
    def snp_index_low(self) -> float:
        """SNP-index for low bulk.

        Returns:
            SNP-index value for low bulk.
        """
        return self.af_low

    @property
    def ref_count_high(self) -> int:
        """Approximate reference allele count in high bulk."""
        return int(round((1 - self.af_high) * self.dp_high))

    @property
    def alt_count_high(self) -> int:
        """Approximate alternate allele count in high bulk."""
        return int(round(self.af_high * self.dp_high))

    @property
    def ref_count_low(self) -> int:
        """Approximate reference allele count in low bulk."""
        return int(round((1 - self.af_low) * self.dp_low))

    @property
    def alt_count_low(self) -> int:
        """Approximate alternate allele count in low bulk."""
        return int(round(self.af_low * self.dp_low))

    @property
    def g_statistic(self) -> float:
        """Calculate G-statistic for this variant.

        The G-statistic is a measure of allele frequency divergence
        between the two bulks. Higher values indicate stronger
        association with the phenotype.

        Returns:
            G-statistic value.
        """
        import math

        # Avoid log(0) by adding small epsilon
        eps = 1e-10

        # Calculate expected frequencies (pooled)
        total_dp = self.dp_high + self.dp_low
        if total_dp == 0:
            return 0.0

        # Observed counts (approximate from AF and DP)
        obs_high_alt = self.af_high * self.dp_high
        obs_high_ref = (1 - self.af_high) * self.dp_high
        obs_low_alt = self.af_low * self.dp_low
        obs_low_ref = (1 - self.af_low) * self.dp_low

        # Expected under null (equal AF in both bulks)
        pooled_af = (obs_high_alt + obs_low_alt) / total_dp
        exp_high_alt = pooled_af * self.dp_high
        exp_high_ref = (1 - pooled_af) * self.dp_high
        exp_low_alt = pooled_af * self.dp_low
        exp_low_ref = (1 - pooled_af) * self.dp_low

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

    def __repr__(self) -> str:
        """Return string representation of variant."""
        return (
            f"Variant({self.chrom}:{self.pos} {self.ref}>{self.alt}, "
            f"AF_high={self.af_high:.3f}, AF_low={self.af_low:.3f}, "
            f"delta={self.delta_af:.3f})"
        )


@dataclass
class Window:
    """Genomic window with aggregated statistics.

    Represents a sliding or fixed window across the genome with
    summary statistics from variants within the window.

    Note:
        This class does not use slots=True to allow modification of
        z_score after initial creation during post-processing.

    Attributes:
        chrom: Chromosome name.
        start: Window start position (1-based, inclusive).
        end: Window end position (1-based, inclusive).
        n_variants: Number of variants in the window.
        mean_delta_af: Mean allele frequency difference.
        median_delta_af: Median allele frequency difference.
        mean_af_high: Mean AF in high bulk.
        mean_af_low: Mean AF in low bulk.
        tricube_delta_af: Tricube-weighted delta AF (for LOESS smoothing).
        g_statistic: G-statistic for the window (pooled across variants).
        p_value: P-value from G-statistic (chi-squared, df=1).
        z_score: Z-score for statistical significance (added in post-processing).
    """

    chrom: str
    start: int
    end: int
    n_variants: int
    mean_delta_af: float
    median_delta_af: float
    mean_af_high: float
    mean_af_low: float
    tricube_delta_af: float
    g_statistic: float
    p_value: float
    z_score: float | None = None

    @property
    def midpoint(self) -> int:
        """Calculate the midpoint of the window.

        Returns:
            Midpoint position of the window.
        """
        return (self.start + self.end) // 2

    @property
    def size(self) -> int:
        """Calculate the size of the window in base pairs.

        Returns:
            Window size in bp.
        """
        return self.end - self.start + 1

    def __repr__(self) -> str:
        """Return string representation of window."""
        return (
            f"Window({self.chrom}:{self.start}-{self.end}, "
            f"n={self.n_variants}, mean_delta={self.mean_delta_af:.3f})"
        )
