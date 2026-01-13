"""I/O utilities for bsaseq."""

from bsaseq.io.readers import (
    read_candidate_variants_tsv,
    read_regions_tsv,
    read_variants_tsv,
    read_windows_tsv,
)
from bsaseq.io.summary import AnalysisSummary, write_summary
from bsaseq.io.vcf import (
    get_sample_names,
    parse_vcf,
    parse_vcf_to_list,
    write_candidates_vcf,
)
from bsaseq.io.writers import (
    write_annotated_candidates_tsv,
    write_candidate_genes_tsv,
    write_candidate_variants_tsv,
    write_regions_bed,
    write_regions_tsv,
    write_variants_tsv,
    write_windows_tsv,
)

__all__ = [
    # VCF
    "get_sample_names",
    "parse_vcf",
    "parse_vcf_to_list",
    "write_candidates_vcf",
    # Readers
    "read_candidate_variants_tsv",
    "read_regions_tsv",
    "read_variants_tsv",
    "read_windows_tsv",
    # Writers
    "write_annotated_candidates_tsv",
    "write_candidate_genes_tsv",
    "write_candidate_variants_tsv",
    "write_regions_bed",
    "write_regions_tsv",
    "write_variants_tsv",
    "write_windows_tsv",
    # Summary
    "AnalysisSummary",
    "write_summary",
]
