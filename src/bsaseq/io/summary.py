"""Summary report generation for bsaseq.

This module provides functions for generating analysis summary reports.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

from bsaseq import __version__
from bsaseq.utils.logging import get_logger

if TYPE_CHECKING:
    from bsaseq.analysis.candidates import CandidateRegion

logger = get_logger(__name__)


@dataclass
class AnalysisSummary:
    """Summary statistics for the analysis.

    Attributes:
        vcf_path: Path to the input VCF file.
        high_bulk_samples: List of sample names in high bulk.
        low_bulk_samples: List of sample names in low bulk.
        total_variants: Total variants in VCF before filtering.
        passing_variants: Variants passing quality filters.
        total_windows: Total windows analyzed.
        significant_windows: Windows above Z-score threshold.
        candidate_regions: Number of candidate regions identified.
        candidate_variants: Number of candidate causal variants.
        top_region: The top-ranked candidate region (or None).
        parameters: Dictionary of all analysis parameters.
        annotated_candidates: Number of candidates with annotations (if annotated).
        candidate_genes: Number of unique candidate genes (if annotated).
        lof_variants: Number of loss-of-function variants (if annotated).
        top_gene: Top candidate gene name (if annotated).
    """

    vcf_path: str
    high_bulk_samples: list[str]
    low_bulk_samples: list[str]
    total_variants: int
    passing_variants: int
    total_windows: int
    significant_windows: int
    candidate_regions: int
    candidate_variants: int
    top_region: CandidateRegion | None
    parameters: dict[str, Any] = field(default_factory=dict)
    annotated_candidates: int = 0
    candidate_genes: int = 0
    lof_variants: int = 0
    top_gene: str | None = None

    def to_text(self) -> str:
        """Format summary as human-readable text.

        Returns:
            Multi-line string with summary information.
        """
        lines = [
            "=" * 70,
            "BSAseq Analysis Summary",
            "=" * 70,
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Version: {__version__}",
            "",
            "INPUT",
            "-" * 70,
            f"VCF file: {self.vcf_path}",
            f"High bulk samples: {', '.join(self.high_bulk_samples)}",
            f"Low bulk samples: {', '.join(self.low_bulk_samples)}",
            "",
            "VARIANT FILTERING",
            "-" * 70,
            f"Total variants in VCF: {self.total_variants:,}",
            f"Variants passing filters: {self.passing_variants:,}",
            f"Filter rate: {100 * (1 - self.passing_variants / max(self.total_variants, 1)):.1f}%",
            "",
            "WINDOW ANALYSIS",
            "-" * 70,
            f"Total windows: {self.total_windows:,}",
            f"Significant windows (z >= {self.parameters.get('z_threshold', 'N/A')}): {self.significant_windows:,}",
            f"Window size: {self.parameters.get('window_size', 'N/A'):,} bp",
            f"Step size: {self.parameters.get('step_size', 'N/A'):,} bp",
            "",
            "CANDIDATE REGIONS",
            "-" * 70,
            f"Candidate regions identified: {self.candidate_regions}",
        ]

        if self.top_region:
            lines.extend([
                "",
                "Top candidate region:",
                f"  Location: {self.top_region.chrom}:{self.top_region.start:,}-{self.top_region.end:,}",
                f"  Length: {self.top_region.length:,} bp",
                f"  Max Z-score: {self.top_region.max_z_score:.2f}",
                f"  Mean delta-AF: {self.top_region.mean_delta_af:.3f}",
                f"  Peak position: {self.top_region.peak_position:,}",
                f"  Windows: {self.top_region.n_windows}",
                f"  Variants: {self.top_region.n_variants}",
            ])

        lines.extend([
            "",
            "CANDIDATE VARIANTS",
            "-" * 70,
            f"Candidate causal variants: {self.candidate_variants}",
            f"Inheritance mode: {self.parameters.get('mode', 'N/A')}",
        ])

        if self.parameters.get("mode") == "recessive":
            lines.append("  (Expected: AF_high ~ 1.0, AF_low ~ 0.0)")
        elif self.parameters.get("mode") == "dominant":
            lines.append("  (Expected: AF_high ~ 0.5+, AF_low ~ 0.0)")

        # Add annotation section if candidates were annotated
        if self.annotated_candidates > 0 or self.candidate_genes > 0:
            lines.extend([
                "",
                "ANNOTATION (snpEff)",
                "-" * 70,
                f"Annotated candidates: {self.annotated_candidates}",
                f"Candidate genes: {self.candidate_genes}",
                f"Loss-of-function variants: {self.lof_variants}",
            ])
            if self.top_gene:
                lines.append(f"Top candidate gene: {self.top_gene}")

        lines.extend([
            "",
            "PARAMETERS",
            "-" * 70,
        ])

        for key, value in sorted(self.parameters.items()):
            if isinstance(value, float):
                lines.append(f"  {key}: {value:.4f}")
            elif isinstance(value, int):
                lines.append(f"  {key}: {value:,}")
            else:
                lines.append(f"  {key}: {value}")

        lines.extend([
            "",
            "=" * 70,
        ])

        return "\n".join(lines)

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for JSON serialization.

        Returns:
            Dictionary representation of the summary.
        """
        result = {
            "version": __version__,
            "timestamp": datetime.now().isoformat(),
            "input": {
                "vcf_path": self.vcf_path,
                "high_bulk_samples": self.high_bulk_samples,
                "low_bulk_samples": self.low_bulk_samples,
            },
            "variants": {
                "total": self.total_variants,
                "passing": self.passing_variants,
            },
            "windows": {
                "total": self.total_windows,
                "significant": self.significant_windows,
            },
            "candidates": {
                "regions": self.candidate_regions,
                "variants": self.candidate_variants,
            },
            "annotation": {
                "annotated_candidates": self.annotated_candidates,
                "candidate_genes": self.candidate_genes,
                "lof_variants": self.lof_variants,
                "top_gene": self.top_gene,
            },
            "parameters": self.parameters,
        }

        if self.top_region:
            result["top_region"] = {
                "chrom": self.top_region.chrom,
                "start": self.top_region.start,
                "end": self.top_region.end,
                "length": self.top_region.length,
                "max_z_score": self.top_region.max_z_score,
                "mean_delta_af": self.top_region.mean_delta_af,
                "peak_position": self.top_region.peak_position,
            }

        return result


def write_summary(summary: AnalysisSummary, path: str | Path) -> None:
    """Write summary report to text file.

    Args:
        summary: AnalysisSummary object to write.
        path: Output file path.
    """
    path = Path(path)
    logger.info(f"Writing summary to {path}")

    with open(path, "w") as f:
        f.write(summary.to_text())
        f.write("\n")

    logger.info(f"Summary written to {path}")
