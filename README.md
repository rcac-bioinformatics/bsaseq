# bsaseq

[![PyPI version](https://badge.fury.io/py/bsaseq.svg)](https://badge.fury.io/py/bsaseq)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Bulk Segregant Analysis for QTL mapping from pooled whole-genome sequencing data.**

bsaseq identifies genomic loci controlling traits by comparing allele frequencies
between phenotypically distinct bulked DNA pools. It provides a complete workflow
from VCF input to annotated candidate genes.

## Features

- Multi-sample bulk support (pool technical replicates)
- Sliding window analysis with tricube smoothing
- Z-score and G-statistic for candidate region detection
- Publication-quality genome-wide and regional plots
- Optional snpEff integration for variant annotation
- Gene-level candidate ranking

## Installation

### From PyPI (recommended)

```bash
pip install bsaseq
```

### From source

```bash
git clone https://github.com/username/bsaseq.git
cd bsaseq
pip install -e .
```

### Dependencies

**Required:**
- Python 3.9-3.13
- cyvcf2, numpy, scipy, matplotlib, click, rich

**Optional:**
- snpEff (for variant annotation)

## Quick start

```bash
# Basic analysis
bsaseq run \
    --vcf joint_calls.vcf.gz \
    --high-bulk mutant_pool \
    --low-bulk wildtype_pool \
    --out results/my_analysis

# With multiple samples per bulk
bsaseq run \
    --vcf joint_calls.vcf.gz \
    --high-bulk "mut_rep1,mut_rep2" \
    --low-bulk "wt_rep1,wt_rep2" \
    --out results/my_analysis

# With annotation
bsaseq run \
    --vcf joint_calls.vcf.gz \
    --high-bulk mutant_pool \
    --low-bulk wildtype_pool \
    --out results/my_analysis \
    --annotate \
    --snpeff-db Sorghum_bicolor
```

## Commands

| Command | Description |
|---------|-------------|
| `bsaseq run` | Complete BSA analysis pipeline |
| `bsaseq samples` | List sample names in VCF |
| `bsaseq plot` | Regenerate plots from existing output |
| `bsaseq annotate` | Annotate candidates with snpEff |
| `bsaseq check-snpeff` | Verify snpEff installation |

## Output files

| File | Description |
|------|-------------|
| `*_variants.tsv` | Per-variant allele frequencies |
| `*_windows.tsv` | Sliding window statistics |
| `*_regions.tsv` | Candidate genomic regions |
| `*_regions.bed` | Candidate regions in BED format |
| `*_candidates.tsv` | Filtered candidate variants |
| `*_annotated_candidates.tsv` | Candidates with snpEff annotation |
| `*_candidate_genes.tsv` | Gene-level summary |
| `*_summary.txt` | Analysis summary report |
| `*_genome_wide.png/pdf` | Genome-wide Manhattan plot |
| `*_region_*.png/pdf` | Zoomed candidate region plots |
| `*_af_distribution.png` | Allele frequency diagnostics |
| `*_depth_distribution.png` | Read depth diagnostics |

## Parameters

### Filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-dp` | 10 | Minimum read depth per sample |
| `--max-dp` | 200 | Maximum read depth (excludes repeats) |
| `--min-gq` | 20 | Minimum genotype quality |
| `--min-qual` | 30 | Minimum variant QUAL score |

### Window analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--window-size` | 1000000 | Window width in bp |
| `--step-size` | 250000 | Step between windows in bp |
| `--min-variants` | 5 | Minimum variants per window |

### Candidate detection

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--z-threshold` | 3.0 | Z-score cutoff for significance |
| `--mode` | recessive | Inheritance mode (recessive/dominant) |

## Methodology

### Delta allele frequency (delta-AF)

For each SNP, bsaseq calculates the difference in alternate allele frequency
between the high bulk (mutant pool) and low bulk (wild-type pool):

```
delta_AF = AF_high - AF_low
```

For a recessive causal mutation:
- In the mutant bulk: individuals are homozygous mutant, so AF ~ 1.0
- In the wild-type bulk: individuals are homozygous reference, so AF ~ 0.0
- Expected delta_AF ~ 1.0 at the causal locus

### Sliding window analysis

Variants are analyzed in overlapping genomic windows to reduce noise:

1. **Tricube-weighted smoothing**: Variants near window center contribute more
2. **G-statistic**: Tests for significant allele frequency differences
3. **Z-score normalization**: Genome-wide standardization for peak calling

### Candidate region detection

Regions are identified where consecutive windows exceed the Z-score threshold.
Adjacent significant regions within 500 kb are merged.

### Candidate variant filtering

Within candidate regions, variants are filtered based on inheritance mode:

| Mode | min_delta_AF | min_AF_high | max_AF_low |
|------|--------------|-------------|------------|
| Recessive | 0.8 | 0.9 | 0.1 |
| Dominant | 0.3 | 0.4 | 0.1 |

## Input requirements

The VCF file must contain:
- Biallelic SNPs (multiallelic and indels are skipped)
- AD (allelic depth) FORMAT field for allele frequency calculation
- GQ (genotype quality) FORMAT field for filtering (optional but recommended)

Recommended variant calling:
```bash
# GATK HaplotypeCaller (recommended)
gatk HaplotypeCaller -R ref.fa -I mutant_pool.bam -I wildtype_pool.bam -O calls.vcf.gz

# bcftools mpileup (alternative)
bcftools mpileup -f ref.fa mutant_pool.bam wildtype_pool.bam | bcftools call -mv -O z -o calls.vcf.gz
```

## Examples

### Example 1: Basic recessive mutation mapping

```bash
bsaseq run \
    --vcf pooled_calls.vcf.gz \
    --high-bulk mutant \
    --low-bulk wildtype \
    --out results/analysis \
    --mode recessive
```

### Example 2: Dominant trait with relaxed thresholds

```bash
bsaseq run \
    --vcf pooled_calls.vcf.gz \
    --high-bulk affected \
    --low-bulk unaffected \
    --out results/dominant \
    --mode dominant \
    --z-threshold 2.5 \
    --min-dp 5
```

### Example 3: High-coverage data with annotation

```bash
bsaseq run \
    --vcf deep_seq.vcf.gz \
    --high-bulk mut1,mut2 \
    --low-bulk wt1,wt2 \
    --out results/annotated \
    --min-dp 30 \
    --max-dp 500 \
    --annotate \
    --snpeff-db Arabidopsis_thaliana
```

### Example 4: Regenerate plots with different format

```bash
bsaseq plot \
    --windows results/analysis_windows.tsv \
    --regions results/analysis_regions.tsv \
    --variants results/analysis_variants.tsv \
    --out new_plots \
    --format pdf
```

## Development

```bash
# Install with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=bsaseq --cov-report=term-missing

# Type checking
mypy src/bsaseq
```

## Citation

If you use bsaseq in your research, please cite:

> bsaseq: A tool for Bulk Segregant Analysis of pooled sequencing data (in preparation)

## License

MIT License - see [LICENSE](LICENSE) for details.
