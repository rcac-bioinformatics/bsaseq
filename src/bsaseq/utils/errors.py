"""User-friendly error messages for bsaseq.

This module provides formatted error messages with helpful suggestions
for common error conditions encountered during BSA analysis.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from rich.console import Console
from rich.panel import Panel

console = Console(stderr=True)


class BSASeqError(Exception):
    """Base exception for bsaseq errors with user-friendly formatting."""

    def __init__(self, message: str, suggestion: str | None = None):
        """Initialize error with message and optional suggestion.

        Args:
            message: Main error message.
            suggestion: Optional suggestion for how to fix the error.
        """
        self.message = message
        self.suggestion = suggestion
        super().__init__(message)

    def display(self) -> None:
        """Display the error in a formatted panel."""
        content = f"[red bold]Error:[/red bold] {self.message}"
        if self.suggestion:
            content += f"\n\n[yellow]Suggestion:[/yellow] {self.suggestion}"
        console.print(Panel(content, title="bsaseq Error", border_style="red"))


def format_file_not_found(path: str | Path, file_type: str = "File") -> str:
    """Format a file not found error message.

    Args:
        path: Path to the missing file.
        file_type: Type of file (e.g., "VCF file", "Output directory").

    Returns:
        Formatted error message.
    """
    path = Path(path)
    msg = f"{file_type} not found: {path}"

    # Check if parent directory exists
    if not path.parent.exists():
        msg += f"\n\nThe parent directory does not exist: {path.parent}"
        msg += "\nCreate the directory first or check the path."

    return msg


def format_sample_not_found(
    sample_names: Sequence[str],
    available_samples: Sequence[str],
    vcf_path: str | Path,
) -> str:
    """Format a sample not found error message.

    Args:
        sample_names: Sample names that were not found.
        available_samples: List of available sample names.
        vcf_path: Path to the VCF file.

    Returns:
        Formatted error message with suggestions.
    """
    missing = ", ".join(sample_names)
    available = ", ".join(sorted(available_samples))

    msg = f"Sample(s) not found in VCF: {missing}\n\n"
    msg += f"Available samples:\n  {available}\n\n"
    msg += "Hints:\n"
    msg += "  - Sample names are case-sensitive\n"
    msg += f"  - Use 'bsaseq samples --vcf {vcf_path}' to list all samples\n"
    msg += "  - Check for typos or extra whitespace in sample names"

    return msg


def format_sample_overlap(overlapping: Sequence[str]) -> str:
    """Format an error for samples appearing in both bulks.

    Args:
        overlapping: Sample names that appear in both bulks.

    Returns:
        Formatted error message.
    """
    samples = ", ".join(overlapping)
    msg = f"Sample(s) cannot be in both high and low bulk: {samples}\n\n"
    msg += "Each sample must be assigned to exactly one bulk.\n"
    msg += "Check your --high-bulk and --low-bulk arguments."

    return msg


def format_missing_format_field(field: str, description: str) -> str:
    """Format an error for missing VCF FORMAT field.

    Args:
        field: The missing field name (e.g., "AD", "GT").
        description: Description of what the field is used for.

    Returns:
        Formatted error message with suggestions.
    """
    msg = f"VCF does not contain {field} FORMAT field.\n\n"
    msg += f"{description}\n\n"

    if field == "AD":
        msg += "The AD (allelic depth) field is required to calculate allele frequencies.\n"
        msg += "Ensure your VCF was generated with a variant caller that outputs per-allele\n"
        msg += "read depths, such as:\n"
        msg += "  - GATK HaplotypeCaller\n"
        msg += "  - bcftools mpileup with -a AD option\n"
        msg += "  - FreeBayes\n"
    elif field == "GT":
        msg += "The GT (genotype) field is a standard VCF field required for variant analysis.\n"
        msg += "This may indicate a malformed VCF file."

    return msg


def format_insufficient_samples(found: int, required: int = 2) -> str:
    """Format an error for insufficient samples in VCF.

    Args:
        found: Number of samples found.
        required: Minimum number of samples required.

    Returns:
        Formatted error message.
    """
    msg = f"VCF must contain at least {required} samples for BSA analysis.\n"
    msg += f"Found {found} sample(s).\n\n"
    msg += "BSA requires comparing allele frequencies between two bulked DNA pools.\n"
    msg += "Ensure your VCF contains samples from both the mutant and wild-type pools."

    return msg


def format_invalid_parameter(
    param_name: str,
    value: int | float | str,
    reason: str,
    suggestion: str | None = None,
) -> str:
    """Format an error for an invalid parameter value.

    Args:
        param_name: Name of the parameter.
        value: Invalid value provided.
        reason: Why the value is invalid.
        suggestion: Optional suggestion for valid values.

    Returns:
        Formatted error message.
    """
    msg = f"Invalid value for {param_name}: {value}\n\n"
    msg += f"Reason: {reason}"

    if suggestion:
        msg += f"\n\nSuggestion: {suggestion}"

    return msg


def format_output_error(path: str | Path, reason: str) -> str:
    """Format an error for output path issues.

    Args:
        path: Output path that caused the error.
        reason: Reason for the error.

    Returns:
        Formatted error message.
    """
    msg = f"Cannot write output to: {path}\n\n"
    msg += f"Reason: {reason}\n\n"
    msg += "Suggestions:\n"
    msg += "  - Check that you have write permissions to the directory\n"
    msg += "  - Ensure the parent directory exists\n"
    msg += "  - Check available disk space"

    return msg


def format_no_variants_error(reason: str) -> str:
    """Format an error for no variants passing filters.

    Args:
        reason: Specific reason why no variants passed.

    Returns:
        Formatted error message with suggestions.
    """
    msg = "No variants remaining after filtering.\n\n"
    msg += f"Reason: {reason}\n\n"
    msg += "Suggestions:\n"
    msg += "  - Lower the --min-dp threshold (current filter may be too strict)\n"
    msg += "  - Raise the --max-dp threshold (may be filtering high-coverage regions)\n"
    msg += "  - Lower the --min-qual threshold\n"
    msg += "  - Check that your VCF contains biallelic SNPs (indels are skipped)\n"
    msg += "  - Verify that samples have sufficient coverage"

    return msg


def format_no_candidates_error() -> str:
    """Format an error for no candidate regions found.

    Returns:
        Formatted error message with suggestions.
    """
    msg = "No candidate regions identified.\n\n"
    msg += "This may indicate:\n"
    msg += "  - No significant linkage between markers and the trait\n"
    msg += "  - Insufficient variant density in the region of interest\n"
    msg += "  - Filter thresholds may be too strict\n\n"
    msg += "Suggestions:\n"
    msg += "  - Lower the --z-threshold (default: 3.0) to detect weaker signals\n"
    msg += "  - Increase --window-size to smooth over sparse regions\n"
    msg += "  - Check the diagnostic plots for delta-AF distribution\n"
    msg += "  - Verify that bulk pools are phenotypically distinct"

    return msg


def format_snpeff_error(reason: str, command: str | None = None) -> str:
    """Format an error for snpEff issues.

    Args:
        reason: Why snpEff failed.
        command: Optional failed command for debugging.

    Returns:
        Formatted error message.
    """
    msg = "snpEff annotation failed.\n\n"
    msg += f"Reason: {reason}\n\n"
    msg += "Suggestions:\n"
    msg += "  - Run 'bsaseq check-snpeff' to verify installation\n"
    msg += "  - Ensure the --snpeff-db matches your reference genome\n"
    msg += "  - Check that Java is installed and in PATH\n"
    msg += "  - Try running snpEff manually to diagnose issues"

    if command:
        msg += f"\n\nFailed command:\n  {command}"

    return msg


def display_error(message: str, suggestion: str | None = None) -> None:
    """Display an error message in a formatted panel.

    Args:
        message: Main error message.
        suggestion: Optional suggestion for fixing the error.
    """
    content = f"[red bold]Error:[/red bold] {message}"
    if suggestion:
        content += f"\n\n[yellow]Suggestion:[/yellow] {suggestion}"
    console.print(Panel(content, title="bsaseq Error", border_style="red"))


def display_warning(message: str) -> None:
    """Display a warning message.

    Args:
        message: Warning message to display.
    """
    console.print(f"[yellow]Warning:[/yellow] {message}")
