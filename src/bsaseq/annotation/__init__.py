"""Annotation modules for bsaseq.

This package provides functions for annotating variants using
external tools like snpEff.
"""

from bsaseq.annotation.snpeff import (
    VariantAnnotation,
    VariantImpact,
    check_snpeff_available,
    get_best_annotation,
    get_snpeff_version,
    list_snpeff_databases,
    parse_snpeff_vcf,
    run_snpeff,
)

__all__ = [
    "VariantAnnotation",
    "VariantImpact",
    "check_snpeff_available",
    "get_best_annotation",
    "get_snpeff_version",
    "list_snpeff_databases",
    "parse_snpeff_vcf",
    "run_snpeff",
]
