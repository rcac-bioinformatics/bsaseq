"""snpEff annotation wrapper for bsaseq.

This module provides functions for running snpEff variant annotation
and parsing the results.
"""

from __future__ import annotations

import shutil
import subprocess
from collections.abc import Iterator
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING

from bsaseq.utils.logging import get_logger

if TYPE_CHECKING:
    pass

logger = get_logger(__name__)


class VariantImpact(str, Enum):
    """snpEff impact categories."""

    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    MODIFIER = "MODIFIER"

    @property
    def rank(self) -> int:
        """Numeric rank for sorting (lower = more severe)."""
        return {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}[self.value]


@dataclass
class VariantAnnotation:
    """Parsed snpEff annotation for a single variant effect.

    Attributes:
        chrom: Chromosome name.
        pos: Genomic position.
        ref: Reference allele.
        alt: Alternate allele.
        effect: Effect type (e.g., "missense_variant", "stop_gained").
        impact: Impact severity category.
        gene_name: Gene symbol.
        gene_id: Gene identifier.
        feature_type: Feature type (e.g., "transcript").
        feature_id: Feature/transcript identifier.
        biotype: Gene biotype (e.g., "protein_coding").
        hgvs_c: Coding DNA change (e.g., "c.123A>G").
        hgvs_p: Protein change (e.g., "p.Lys41Arg").
        cdna_pos: Position in cDNA.
        cds_pos: Position in CDS.
        protein_pos: Position in protein.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    effect: str
    impact: VariantImpact
    gene_name: str
    gene_id: str
    feature_type: str
    feature_id: str
    biotype: str
    hgvs_c: str
    hgvs_p: str
    cdna_pos: str
    cds_pos: str
    protein_pos: str

    @property
    def is_coding(self) -> bool:
        """Check if variant affects coding sequence."""
        return self.biotype == "protein_coding" and self.cds_pos != ""

    @property
    def is_lof(self) -> bool:
        """Check if variant is likely loss-of-function."""
        lof_effects = {
            "stop_gained",
            "stop_lost",
            "start_lost",
            "frameshift_variant",
            "splice_acceptor_variant",
            "splice_donor_variant",
            "transcript_ablation",
        }
        return self.effect in lof_effects

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"VariantAnnotation({self.chrom}:{self.pos} {self.ref}>{self.alt}, "
            f"{self.effect}, {self.impact.value}, {self.gene_name})"
        )


def check_snpeff_available() -> bool:
    """Check if snpEff is available in PATH.

    Returns:
        True if snpEff is found, False otherwise.
    """
    return shutil.which("snpEff") is not None


def get_snpeff_version() -> str | None:
    """Get snpEff version string.

    Returns:
        Version string if snpEff is available, None otherwise.
    """
    if not check_snpeff_available():
        return None
    try:
        result = subprocess.run(
            ["snpEff", "-version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        # Parse version from output (usually in stderr)
        output = result.stderr or result.stdout
        for line in output.split("\n"):
            if "snpEff" in line:
                return line.strip()
        # Return first non-empty line
        for line in output.split("\n"):
            if line.strip():
                return line.strip()
        return None
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        return None


def list_snpeff_databases() -> list[str]:
    """List available snpEff databases.

    Returns:
        List of database names.
    """
    if not check_snpeff_available():
        return []
    try:
        result = subprocess.run(
            ["snpEff", "databases"],
            capture_output=True,
            text=True,
            timeout=60,
        )
        databases = []
        for line in result.stdout.split("\n"):
            if line.strip() and not line.startswith("#"):
                parts = line.split("\t")
                if parts:
                    db_name = parts[0].strip()
                    if db_name:
                        databases.append(db_name)
        return databases
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        return []


def run_snpeff(
    vcf_path: str | Path,
    database: str,
    output_path: str | Path | None = None,
    extra_args: list[str] | None = None,
    java_mem: str = "4g",
    timeout: int = 3600,
) -> Path:
    """Run snpEff annotation on a VCF file.

    Args:
        vcf_path: Input VCF file path.
        database: snpEff database name (e.g., "Sorghum_bicolor").
        output_path: Output VCF path (default: input with .annotated.vcf suffix).
        extra_args: Additional arguments to pass to snpEff.
        java_mem: Java heap memory (e.g., "4g", "8g").
        timeout: Maximum runtime in seconds.

    Returns:
        Path to annotated VCF file.

    Raises:
        FileNotFoundError: If snpEff not available.
        RuntimeError: If snpEff fails.
        subprocess.TimeoutExpired: If timeout exceeded.
    """
    if not check_snpeff_available():
        raise FileNotFoundError(
            "snpEff not found in PATH. Install with: conda install -c bioconda snpeff"
        )

    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"Input VCF file not found: {vcf_path}")

    if output_path is None:
        # Handle .vcf.gz extension
        if vcf_path.suffix == ".gz":
            base = vcf_path.with_suffix("")  # Remove .gz
            output_path = base.with_suffix(".annotated.vcf")
        else:
            output_path = vcf_path.with_suffix(".annotated.vcf")
    output_path = Path(output_path)

    logger.info(f"Running snpEff with database: {database}")
    logger.info(f"Input: {vcf_path}")
    logger.info(f"Output: {output_path}")

    cmd = [
        "snpEff",
        f"-Xmx{java_mem}",
        "-v",  # Verbose
        "-noStats",  # Skip HTML stats (faster)
        database,
        str(vcf_path),
    ]

    if extra_args:
        cmd.extend(extra_args)

    logger.debug(f"Command: {' '.join(cmd)}")

    try:
        with open(output_path, "w") as outfile:
            result = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout,
            )

        if result.returncode != 0:
            raise RuntimeError(f"snpEff failed with code {result.returncode}: {result.stderr}")

        logger.info(f"snpEff annotation complete: {output_path}")
        return output_path

    except subprocess.TimeoutExpired:
        logger.error(f"snpEff timed out after {timeout} seconds")
        raise


def parse_snpeff_vcf(vcf_path: str | Path) -> Iterator[VariantAnnotation]:
    """Parse snpEff-annotated VCF and yield VariantAnnotation objects.

    Parses the ANN field in INFO column according to snpEff format:
    ANN=Allele|Effect|Impact|GeneName|GeneID|FeatureType|FeatureID|
        Biotype|Rank|HGVS.c|HGVS.p|cDNA_pos|CDS_pos|Protein_pos|...

    Args:
        vcf_path: Path to snpEff-annotated VCF.

    Yields:
        VariantAnnotation objects (one per effect, multiple per variant possible).
    """
    from cyvcf2 import VCF

    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"Annotated VCF file not found: {vcf_path}")

    logger.info(f"Parsing snpEff annotations from {vcf_path}")

    vcf = VCF(str(vcf_path))

    annotation_count = 0
    for variant in vcf:
        ann_str = variant.INFO.get("ANN")
        if not ann_str:
            continue

        # Multiple annotations separated by comma
        for ann in ann_str.split(","):
            fields = ann.split("|")
            if len(fields) < 11:
                continue  # Malformed annotation

            try:
                impact = VariantImpact(fields[2])
            except ValueError:
                impact = VariantImpact.MODIFIER

            yield VariantAnnotation(
                chrom=variant.CHROM,
                pos=variant.POS,
                ref=variant.REF,
                alt=fields[0] or (variant.ALT[0] if variant.ALT else ""),
                effect=fields[1],
                impact=impact,
                gene_name=fields[3],
                gene_id=fields[4],
                feature_type=fields[5],
                feature_id=fields[6],
                biotype=fields[7] if len(fields) > 7 else "",
                hgvs_c=fields[9] if len(fields) > 9 else "",
                hgvs_p=fields[10] if len(fields) > 10 else "",
                cdna_pos=fields[11] if len(fields) > 11 else "",
                cds_pos=fields[12] if len(fields) > 12 else "",
                protein_pos=fields[13] if len(fields) > 13 else "",
            )
            annotation_count += 1

    vcf.close()
    logger.info(f"Parsed {annotation_count:,} annotations")


def get_best_annotation(
    annotations: list[VariantAnnotation],
    prioritize_lof: bool = True,
) -> VariantAnnotation | None:
    """Select the best annotation for a variant.

    When multiple annotations exist (multiple transcripts), selects
    the most severe effect.

    Args:
        annotations: List of annotations for a single variant.
        prioritize_lof: If True, prefer LOF annotations over others.

    Returns:
        Best annotation or None if list is empty.
    """
    if not annotations:
        return None

    def sort_key(ann: VariantAnnotation) -> tuple:
        # Lower is better: (lof_priority, impact_rank, gene_name)
        lof_priority = 0 if (prioritize_lof and ann.is_lof) else 1
        return (lof_priority, ann.impact.rank, ann.gene_name or "")

    return min(annotations, key=sort_key)
