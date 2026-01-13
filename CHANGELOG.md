# Changelog

All notable changes to bsaseq will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2025-01-XX

### Added

#### Core Analysis
- Sliding window analysis with configurable window and step sizes
- Tricube-weighted smoothing for delta allele frequency
- G-statistic calculation for each window
- Z-score normalization for genome-wide significance testing
- Candidate region detection with automatic merging of adjacent peaks
- Candidate variant filtering for recessive and dominant inheritance modes

#### Multi-sample Support
- Pool multiple samples per bulk (e.g., technical replicates)
- Weighted allele frequency calculation across pooled samples
- Sample name validation and overlap detection

#### Visualization
- Publication-quality genome-wide Manhattan plots
- Regional zoom plots for candidate intervals
- Diagnostic plots for allele frequency and depth distributions
- Support for PNG and PDF output formats
- Customizable z-threshold highlighting

#### Annotation
- Optional snpEff integration for variant effect prediction
- Gene-level candidate ranking based on impact and position
- Annotated candidates output with effect predictions
- `check-snpeff` command for installation verification

#### Input/Output
- VCF input with automatic format validation
- TSV output for variants, windows, regions, and candidates
- BED output for candidate regions
- Comprehensive analysis summary reports

#### CLI Features
- Intuitive command-line interface with Click
- Rich progress bars and formatted output
- Comprehensive help text with examples
- Input validation with helpful error messages
- `samples` command for listing VCF sample names
- `plot` command for regenerating plots
- `annotate` command for standalone annotation

#### Packaging & Deployment
- PyPI package with hatchling build system
- Docker container with optional snpEff
- Bioconda recipe for conda installation
- Nextflow pipeline wrapper
- Snakemake pipeline wrapper
- GitHub Actions CI/CD pipeline

### Dependencies
- Python 3.9-3.12
- cyvcf2 >= 0.30.0
- numpy >= 1.21.0
- pandas >= 1.4.0
- scipy >= 1.7.0
- matplotlib >= 3.5.0
- click >= 8.0.0
- rich >= 12.0.0

### Documentation
- Comprehensive README with methodology and examples
- CLI help text with usage examples
- Nextflow and Snakemake workflow documentation
- Contributing guidelines

---

## Version History

### Pre-release Development

#### v0.5.0 - Annotation Integration
- Added snpEff annotation support
- Gene-level candidate ranking
- Annotated output files

#### v0.4.0 - Visualization
- Publication-quality plotting
- Genome-wide Manhattan plots
- Regional zoom plots
- Diagnostic plots

#### v0.3.0 - Candidate Detection
- Candidate region identification
- Variant filtering by inheritance mode
- Region merging algorithm

#### v0.2.0 - Multi-sample Support
- Multiple samples per bulk
- Pooled allele frequency calculation
- Sliding window analysis

#### v0.1.0 - Initial Development
- Basic VCF parsing
- Variant data models
- CLI skeleton

---

[Unreleased]: https://github.com/username/bsaseq/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/username/bsaseq/releases/tag/v1.0.0
