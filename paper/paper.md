---
title: 'bsaseq: A Python package for bulk segregant analysis and QTL mapping from pooled sequencing data'
tags:
  - Python
  - bioinformatics
  - genomics
  - QTL mapping
  - bulk segregant analysis
  - variant analysis
authors:
  - name: Arun Seetharam
    orcid: 0000-0002-6789-9298
    affiliation: 1
affiliations:
  - name: Rosen Center for Advanced Computing, Purdue University, West Lafayette, IN, USA
    index: 1
date: January 2025
bibliography: paper.bib
---

# Summary

Bulk segregant analysis (BSA) is a powerful genetic mapping approach for rapidly identifying genomic loci controlling phenotypic traits. By comparing allele frequencies between pooled DNA samples from individuals with contrasting phenotypes, BSA can localize causal mutations without requiring individual genotyping of large populations.

`bsaseq` is a Python package that provides a complete workflow for BSA from variant call format (VCF) input to ranked candidate genes. The package implements sliding window analysis with tricube kernel smoothing, calculates both Z-scores and G-statistics for significance testing, generates publication-quality genome-wide Manhattan plots, and optionally integrates with snpEff for functional variant annotation. Designed for ease of use, `bsaseq` offers a comprehensive command-line interface with progress reporting and helpful error messages, making it accessible to researchers with varying levels of computational expertise.

# Statement of need

BSA combined with whole-genome sequencing has become a standard approach for QTL mapping in plant genetics, fungal genetics, and model organism research [@takagi2013qtlseq; @magwene2011statistics]. Despite its widespread adoption, researchers face significant barriers when implementing BSA analyses due to limitations in existing software tools.

QTLseqr [@mansfeld2018qtlseqr] provides a robust R implementation but lacks a command-line interface, requires complex R environment setup, and has not been actively maintained since 2021. SHOREmap offers comprehensive functionality but suffers from complex installation requirements involving Perl and C++ dependencies, and sparse documentation. Many researchers resort to custom scripts adapted from published MutMap pipelines, which are difficult to reproduce and maintain.

`bsaseq` addresses these limitations by providing:

- **Accessibility**: Pure Python implementation with minimal dependencies, installable via pip, conda, or Docker
- **Usability**: Intuitive CLI with comprehensive help text, input validation, and informative error messages
- **Reproducibility**: Integration with workflow managers (Nextflow, Snakemake) for scalable, reproducible analyses
- **HPC-ready**: Tested on SLURM clusters with appropriate resource management
- **Active maintenance**: Open-source development with continuous integration and automated testing

# Key features

- Multi-sample bulk support for pooling technical replicates
- Configurable sliding window analysis with tricube kernel smoothing
- Dual statistical framework: Z-scores for relative signal strength and G-statistics for significance testing
- Automatic candidate region detection with adjacent peak merging
- Variant filtering for recessive and dominant inheritance modes
- Publication-quality genome-wide Manhattan plots and regional zoom plots
- Optional snpEff integration for functional variant annotation [@cingolani2012snpeff]
- Gene-level candidate ranking prioritizing high-impact variants
- Comprehensive output formats including TSV, BED, and VCF

# Implementation

`bsaseq` is implemented in Python 3.9+ and leverages established scientific Python libraries: cyvcf2 for efficient VCF parsing, NumPy and SciPy for statistical calculations, pandas for data manipulation, and Matplotlib for visualization. The command-line interface uses Click with Rich for enhanced terminal output including progress bars and formatted tables.

The statistical methodology follows established approaches [@magwene2011statistics; @mansfeld2018qtlseqr]. For each genomic window, the delta allele frequency (high bulk AF minus low bulk AF) is calculated and smoothed using a tricube kernel. Windows are scored using Z-score normalization, and significant regions are identified by consecutive windows exceeding a user-defined threshold (default: Z > 3.0).

# Availability

- **Source code**: https://github.com/rcac-bioinformatics/bsaseq
- **PyPI**: `pip install bsaseq`
- **Bioconda**: `conda install -c bioconda bsaseq`
- **Docker**: `docker pull username/bsaseq`
- **Documentation**: https://github.com/rcac-bioinformatics/bsaseq#readme

# Acknowledgements

The author thanks the Purdue University Research Computing team for computational resources and feedback during development. This work was supported by the Rosen Center for Advanced Computing at Purdue University.

# References
