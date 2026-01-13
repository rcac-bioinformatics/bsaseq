#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * bsaseq Nextflow Pipeline
 * Bulk Segregant Analysis for QTL mapping
 *
 * Usage:
 *   nextflow run main.nf --vcf variants.vcf.gz --high_bulk MUT --low_bulk WT
 */

// Pipeline version
def pipelineVersion = '1.0.0'

// Default parameters
params.vcf = null
params.high_bulk = null
params.low_bulk = null
params.outdir = "results"
params.window_size = 1000000
params.step_size = 250000
params.min_dp = 10
params.max_dp = 200
params.min_gq = 20
params.min_qual = 30
params.min_variants = 5
params.z_threshold = 3.0
params.mode = "recessive"
params.annotate = false
params.snpeff_db = null
params.snpeff_mem = "4g"
params.plot_format = "both"
params.help = false

// Help message
def helpMessage() {
    log.info """
    =========================================
      bsaseq Nextflow Pipeline v${pipelineVersion}
    =========================================

    Usage:
      nextflow run main.nf --vcf <file> --high_bulk <samples> --low_bulk <samples>

    Required:
      --vcf           Input VCF file (bgzipped and indexed)
      --high_bulk     High bulk sample name(s), comma-separated
      --low_bulk      Low bulk sample name(s), comma-separated

    Output:
      --outdir        Output directory [default: ${params.outdir}]

    Analysis parameters:
      --window_size   Window size in bp [default: ${params.window_size}]
      --step_size     Step size in bp [default: ${params.step_size}]
      --min_dp        Minimum read depth [default: ${params.min_dp}]
      --max_dp        Maximum read depth [default: ${params.max_dp}]
      --min_gq        Minimum genotype quality [default: ${params.min_gq}]
      --min_qual      Minimum variant quality [default: ${params.min_qual}]
      --min_variants  Minimum variants per window [default: ${params.min_variants}]
      --z_threshold   Z-score threshold [default: ${params.z_threshold}]
      --mode          Inheritance mode (recessive|dominant) [default: ${params.mode}]

    Annotation:
      --annotate      Enable snpEff annotation [default: ${params.annotate}]
      --snpeff_db     snpEff database name [default: ${params.snpeff_db}]
      --snpeff_mem    Memory for snpEff [default: ${params.snpeff_mem}]

    Plotting:
      --plot_format   Plot format (png|pdf|both) [default: ${params.plot_format}]

    Profiles:
      -profile standard     Docker (default)
      -profile singularity  Singularity/Apptainer
      -profile conda        Conda environment
      -profile slurm        SLURM cluster
      -profile test         Test with example data

    Examples:
      # Basic analysis
      nextflow run main.nf \\
          --vcf variants.vcf.gz \\
          --high_bulk mutant \\
          --low_bulk wildtype

      # With multiple samples per bulk
      nextflow run main.nf \\
          --vcf variants.vcf.gz \\
          --high_bulk "mut1,mut2" \\
          --low_bulk "wt1,wt2"

      # With annotation on HPC
      nextflow run main.nf \\
          -profile slurm,singularity \\
          --vcf variants.vcf.gz \\
          --high_bulk MUT \\
          --low_bulk WT \\
          --annotate \\
          --snpeff_db Arabidopsis_thaliana
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.vcf) {
    log.error "Parameter 'vcf' is required. Use --help for usage information."
    exit 1
}
if (!params.high_bulk) {
    log.error "Parameter 'high_bulk' is required. Use --help for usage information."
    exit 1
}
if (!params.low_bulk) {
    log.error "Parameter 'low_bulk' is required. Use --help for usage information."
    exit 1
}

// Print pipeline info
log.info """
=========================================
  bsaseq Nextflow Pipeline v${pipelineVersion}
=========================================
VCF:          ${params.vcf}
High bulk:    ${params.high_bulk}
Low bulk:     ${params.low_bulk}
Output:       ${params.outdir}
Window size:  ${params.window_size}
Step size:    ${params.step_size}
Z-threshold:  ${params.z_threshold}
Mode:         ${params.mode}
Annotate:     ${params.annotate}
snpEff DB:    ${params.snpeff_db ?: 'N/A'}
=========================================
"""


/*
 * Main analysis process
 */
process BSASEQ_RUN {
    tag "BSA analysis"
    publishDir "${params.outdir}", mode: 'copy'

    container 'username/bsaseq:latest'

    input:
    path vcf
    path vcf_index

    output:
    path "bsa_*", emit: all_outputs
    path "bsa_summary.txt", emit: summary
    path "bsa_variants.tsv", emit: variants
    path "bsa_windows.tsv", emit: windows
    path "bsa_regions.tsv", emit: regions, optional: true
    path "bsa_regions.bed", emit: regions_bed, optional: true
    path "bsa_candidates.tsv", emit: candidates, optional: true
    path "bsa_genome_wide.png", emit: plot_png, optional: true
    path "bsa_genome_wide.pdf", emit: plot_pdf, optional: true

    script:
    def annotate_args = params.annotate && params.snpeff_db ? "--annotate --snpeff-db ${params.snpeff_db} --snpeff-mem ${params.snpeff_mem}" : "--no-annotate"
    """
    bsaseq run \\
        --vcf ${vcf} \\
        --high-bulk "${params.high_bulk}" \\
        --low-bulk "${params.low_bulk}" \\
        --out bsa \\
        --window-size ${params.window_size} \\
        --step-size ${params.step_size} \\
        --min-dp ${params.min_dp} \\
        --max-dp ${params.max_dp} \\
        --min-gq ${params.min_gq} \\
        --min-qual ${params.min_qual} \\
        --min-variants ${params.min_variants} \\
        --z-threshold ${params.z_threshold} \\
        --mode ${params.mode} \\
        ${annotate_args} \\
        --plot \\
        --plot-format ${params.plot_format}
    """
}


/*
 * Workflow
 */
workflow {
    // Create channels for input files
    vcf_ch = Channel.fromPath(params.vcf, checkIfExists: true)

    // Try to find index file (either .tbi or .csi)
    vcf_index_ch = Channel.fromPath("${params.vcf}.tbi", checkIfExists: false)
        .ifEmpty { Channel.fromPath("${params.vcf}.csi", checkIfExists: false) }
        .ifEmpty {
            log.warn "No index file found for ${params.vcf}. Analysis may be slower."
            Channel.empty()
        }

    // Run analysis
    BSASEQ_RUN(vcf_ch, vcf_index_ch.ifEmpty(file('NO_INDEX')))
}


/*
 * Completion handler
 */
workflow.onComplete {
    log.info """
    =========================================
      Pipeline completed!
    =========================================
    Status:     ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:   ${workflow.duration}
    Output:     ${params.outdir}
    =========================================
    """.stripIndent()
}
