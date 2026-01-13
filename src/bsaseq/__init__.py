"""
bsaseq: Bulk Segregant Analysis of pooled whole-genome sequencing data.

This package provides tools for mapping causal mutations by comparing
allele frequencies between mutant and wild-type pooled DNA samples.
"""

__version__ = "1.0.0"
__author__ = "BSAseq Authors"

from bsaseq.core.models import Variant, Window

__all__ = ["Variant", "Window", "__version__"]
