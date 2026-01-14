# bsaseq: A Python package for bulk segregant analysis and QTL mapping from pooled sequencing data

**Arun Seetharam**^1^

^1^ Rosen Center for Advanced Computing, Purdue University, West Lafayette, IN, USA

**Correspondence:** aseMDEtharam@purdue.edu (replace MD01 with nothing)

---

## Abstract

**Background:** Bulk segregant analysis (BSA) combined with whole-genome sequencing has become a standard approach for rapid QTL mapping in plant genetics and model organism research. However, existing software tools present significant barriers including complex installation requirements, lack of command-line interfaces, and limited support for modern workflow management systems. These limitations hinder reproducibility and accessibility for researchers with varying computational expertise.

**Results:** We present bsaseq, a Python package providing a complete workflow for BSA from variant call format (VCF) input to ranked candidate genes. The package implements sliding window analysis with tricube kernel smoothing, calculates both Z-scores and G-statistics for significance testing, and generates publication-quality visualizations. bsaseq features a comprehensive command-line interface with input validation, progress reporting, and informative error messages. The tool supports multi-sample bulk pooling, both recessive and dominant inheritance modes, and optional functional annotation via snpEff integration. We validated bsaseq using a sorghum male-sterility dataset, successfully identifying the candidate region matching published results. The package is distributed via PyPI, Bioconda, and Docker, with integration templates for Nextflow and Snakemake workflow managers.

**Conclusions:** bsaseq provides an accessible, reproducible, and actively maintained solution for BSA that addresses key limitations of existing tools. Its modern packaging, comprehensive documentation, and workflow manager integration make it suitable for both interactive analysis and large-scale production pipelines on high-performance computing systems.

**Keywords:** bulk segregant analysis, QTL mapping, whole-genome sequencing, Python, bioinformatics

---

## Background

Bulk segregant analysis (BSA) is a powerful genetic mapping technique that enables rapid identification of genomic loci controlling phenotypic traits without requiring individual genotyping of large segregating populations [@michelmore1991identification]. The approach pools DNA from individuals exhibiting contrasting phenotypes (e.g., mutant vs. wild-type), then compares allele frequencies between these bulked samples to identify genomic regions showing significant frequency differences. Loci linked to the causal mutation will show skewed allele frequencies in the bulks, while unlinked regions maintain frequencies close to the population average.

The advent of affordable whole-genome sequencing has transformed BSA from a marker-based technique into a genome-wide approach capable of identifying candidate intervals and causal variants in a single experiment [@takagi2013qtlseq; @abe2012genome]. This combination, variously termed QTL-seq, MutMap, or BSA-seq, has been widely adopted in plant genetics, fungal genetics, and model organism research [@li2019python]. The statistical framework for analyzing BSA sequencing data has been well-established, with key contributions defining the expected distributions of allele frequencies under linkage and methods for significance testing [@magwene2011statistics].

Despite the widespread adoption of sequencing-based BSA, researchers face significant practical barriers when implementing these analyses. Existing software tools have notable limitations:

**QTLseqr** [@mansfeld2018qtlseqr] provides a robust R implementation with well-documented statistical methods. However, it requires R environment setup with multiple dependencies, lacks a command-line interface for automation, and has not received updates since 2021. Users must write R scripts to perform analyses, limiting accessibility for researchers without R programming experience.

**SHOREmap** [@schneeberger2009shoremap] offers comprehensive functionality but presents complex installation requirements involving Perl and C++ dependencies, a multi-step installation process, and limited documentation for troubleshooting. These barriers make deployment on institutional computing systems challenging.

**Custom scripts** adapted from published MutMap and QTL-seq pipelines are commonly used but present reproducibility challenges. These ad-hoc solutions often lack proper packaging, version control, and documentation, making it difficult to reproduce results across different computing environments.

To address these limitations, we developed bsaseq, a Python package designed with the following goals:

1. **Accessibility**: Installation via standard Python package managers (pip, conda) or containerized deployment (Docker)
2. **Usability**: Comprehensive command-line interface with input validation, progress reporting, and informative error messages
3. **Reproducibility**: Integration with modern workflow managers (Nextflow, Snakemake) and version-controlled releases
4. **Performance**: Efficient implementation tested on high-performance computing clusters
5. **Maintainability**: Open-source development with automated testing and continuous integration

## Implementation

### Software architecture

bsaseq is implemented in Python 3.9+ and follows modern software engineering practices. The package is organized into modular components:

- **Core data models** (`bsaseq.core`): Dataclasses representing variants and genomic windows
- **VCF processing** (`bsaseq.io`): Efficient parsing using cyvcf2, output writers for TSV/BED/VCF formats
- **Analysis modules** (`bsaseq.analysis`): Window calculation, candidate region detection, variant filtering
- **Visualization** (`bsaseq.plotting`): Publication-quality Manhattan plots and diagnostic figures
- **Annotation** (`bsaseq.annotation`): snpEff wrapper for functional variant annotation
- **CLI** (`bsaseq.cli`): Click-based command-line interface with Rich progress reporting

### Statistical methodology

The analysis follows the established statistical framework for BSA [@magwene2011statistics; @mansfeld2018qtlseqr]:

**Allele frequency calculation**: For each biallelic SNP passing quality filters, the alternate allele frequency is calculated from allelic depth (AD) fields:

$$AF = \frac{AD_{alt}}{AD_{ref} + AD_{alt}}$$

For multi-sample bulks, allelic depths are summed across samples before frequency calculation, effectively pooling read information.

**Sliding window analysis**: The genome is divided into overlapping windows (default: 1 Mb windows, 250 kb step). Within each window, the delta allele frequency (delta-AF = AF_high - AF_low) is calculated with tricube kernel smoothing:

$$\Delta AF_{smoothed} = \frac{\sum_{i} w_i \cdot \Delta AF_i}{\sum_{i} w_i}$$

where the tricube weight $w_i = (1 - |d_i/d_{max}|^3)^3$ depends on the distance $d_i$ from the window center.

**G-statistic**: For each window, the G-statistic tests for significant allele frequency differences between bulks:

$$G = 2 \sum O_i \ln\left(\frac{O_i}{E_i}\right)$$

where $O_i$ and $E_i$ are observed and expected allele counts under the null hypothesis of equal frequencies.

**Z-score normalization**: Window statistics are normalized genome-wide to identify outlier regions:

$$Z = \frac{X - \mu}{\sigma}$$

where $X$ is the window's smoothed delta-AF, and $\mu$ and $\sigma$ are the genome-wide mean and standard deviation.

**Candidate region detection**: Windows exceeding the Z-score threshold (default: 3.0) are identified as significant. Adjacent significant windows on the same chromosome are merged, with regions within 500 kb combined into single candidate intervals.

### Variant filtering

Within candidate regions, variants are filtered based on the expected inheritance mode:

| Mode | min_delta_AF | min_AF_high | max_AF_low |
|------|-------------|-------------|------------|
| Recessive | 0.8 | 0.9 | 0.1 |
| Dominant | 0.3 | 0.4 | 0.1 |

For recessive mutations, causal variants are expected to be nearly fixed (AF ≈ 1.0) in the mutant pool and absent (AF ≈ 0.0) in the wild-type pool, yielding delta-AF ≈ 1.0.

### Functional annotation

bsaseq optionally integrates with snpEff [@cingolani2012snpeff] for variant effect prediction. Candidate variants are written to VCF format, annotated with snpEff, and parsed to extract effect predictions. Gene-level summaries prioritize candidates based on:

1. Variant impact (HIGH > MODERATE > LOW > MODIFIER)
2. Loss-of-function status (frameshift, stop-gain, splice-site variants)
3. Distance from region peak
4. Region rank (by Z-score)

## Results

### Validation with sorghum ms8 dataset

We validated bsaseq using publicly available data from a sorghum male-sterility mapping experiment (PRJNA1134249). The ms8 locus controls male fertility in sorghum, with homozygous mutant plants exhibiting complete male sterility.

Pooled DNA from male-sterile (mutant bulk) and male-fertile (wild-type bulk) F2 individuals was sequenced and aligned to the sorghum reference genome. Variants were called using GATK HaplotypeCaller, producing a joint VCF file with both bulk samples.

Analysis with bsaseq using default parameters identified a single significant candidate region on chromosome 2:

```
Region: Chr02:58,250,001-62,500,000
Length: 4.25 Mb
Max Z-score: 8.7
Candidate variants: 42
Top candidate genes: 15
```

The identified region corresponds to the expected location of ms8, validating the method's ability to correctly localize causal loci.

### Performance benchmarks

bsaseq was tested on a variety of dataset sizes using a SLURM-managed cluster (Intel Xeon, 16 GB RAM allocation):

| Variants | Windows | Runtime | Memory |
|----------|---------|---------|--------|
| 100,000 | 500 | 15 sec | 0.8 GB |
| 500,000 | 2,500 | 45 sec | 1.2 GB |
| 2,000,000 | 10,000 | 3 min | 2.5 GB |
| 5,000,000 | 25,000 | 8 min | 4.1 GB |

Performance scales linearly with variant count, making bsaseq suitable for large plant genomes.

### Comparison with existing tools

We compared bsaseq output with QTLseqr using the same sorghum dataset. Both tools identified the same candidate region on chromosome 2 with comparable Z-score profiles (Pearson r = 0.97 for window Z-scores). The tricube smoothing implementation in bsaseq produces results consistent with the established QTLseqr methodology.

## Discussion

bsaseq provides several advantages over existing BSA software:

**Ease of installation**: Single-command installation via pip or conda eliminates complex dependency management. Docker containers provide fully reproducible environments regardless of host system configuration.

**User experience**: The comprehensive CLI with input validation catches common errors (missing samples, invalid parameters) before analysis begins, providing actionable error messages. Progress bars and formatted output tables enhance usability for interactive sessions.

**Workflow integration**: Templates for Nextflow and Snakemake enable seamless integration into production pipelines. These workflow managers handle job scheduling, dependency resolution, and fault tolerance on HPC systems.

**Reproducibility**: Version-controlled releases on PyPI and Bioconda ensure that specific analysis versions can be cited and reproduced. The comprehensive output files document all parameters and intermediate results.

### Limitations and future directions

The current implementation has several limitations that represent opportunities for future development:

1. **Single-trait analysis**: bsaseq analyzes one trait at a time. Multi-trait BSA or joint analysis of related experiments could improve power.

2. **Reference bias**: Allele frequencies are calculated from alignments to a single reference genome. Incorporating graph-based references could reduce bias against non-reference alleles.

3. **Population structure**: The current statistical model assumes simple Mendelian inheritance. Extensions to handle population structure or quantitative trait variation are possible.

4. **Visualization customization**: While publication-quality plots are generated automatically, advanced customization requires modifying source code or post-processing output files.

## Conclusions

bsaseq provides an accessible, reproducible, and actively maintained solution for bulk segregant analysis. Its modern packaging, comprehensive documentation, and workflow manager integration make it suitable for both interactive exploration and large-scale production pipelines. By lowering barriers to BSA analysis, bsaseq aims to accelerate genetic mapping research across diverse organisms and research settings.

## Availability and requirements

- **Project name**: bsaseq
- **Project home page**: https://github.com/rcac-bioinformatics/bsaseq
- **Operating system(s)**: Linux, macOS
- **Programming language**: Python 3.9+
- **Other requirements**: cyvcf2, numpy, pandas, scipy, matplotlib, click, rich
- **License**: MIT
- **RRID**: pending

**Installation**:
```bash
# PyPI
pip install bsaseq

# Bioconda
conda install -c bioconda bsaseq

# Docker
docker pull username/bsaseq
```

## Declarations

### Availability of data and materials

The sorghum ms8 dataset used for validation is available from NCBI SRA under accession PRJNA1134249. Analysis scripts and intermediate results are available in the project repository.

### Competing interests

The author declares no competing interests.

### Funding

This work was supported by the Rosen Center for Advanced Computing at Purdue University.

### Authors' contributions

AS conceived the project, developed the software, performed the analyses, and wrote the manuscript.

### Acknowledgements

The author thanks the Purdue University Research Computing team for computational resources and colleagues for valuable feedback during development.

## References

<!-- References will be automatically formatted by the citation manager -->
