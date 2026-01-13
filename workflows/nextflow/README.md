# bsaseq Nextflow Pipeline

Run bsaseq as a Nextflow workflow for reproducible, scalable analysis.

## Quick start

```bash
nextflow run main.nf \
    --vcf /path/to/variants.vcf.gz \
    --high_bulk "mut1,mut2" \
    --low_bulk "wt1,wt2" \
    --outdir results
```

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04.0
- One of: Docker, Singularity/Apptainer, or Conda

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--vcf` | Input VCF file (bgzipped recommended) |
| `--high_bulk` | High bulk sample name(s), comma-separated |
| `--low_bulk` | Low bulk sample name(s), comma-separated |

### Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | results | Output directory |

### Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--window_size` | 1000000 | Window size (bp) |
| `--step_size` | 250000 | Step between windows (bp) |
| `--min_dp` | 10 | Minimum read depth |
| `--max_dp` | 200 | Maximum read depth |
| `--min_gq` | 20 | Minimum genotype quality |
| `--min_qual` | 30 | Minimum variant quality |
| `--min_variants` | 5 | Minimum variants per window |
| `--z_threshold` | 3.0 | Z-score threshold |
| `--mode` | recessive | Inheritance mode (recessive/dominant) |

### Annotation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--annotate` | false | Enable snpEff annotation |
| `--snpeff_db` | null | snpEff database name |
| `--snpeff_mem` | 4g | Memory for snpEff |

### Plotting

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--plot_format` | both | Plot format (png/pdf/both) |

## Profiles

| Profile | Description |
|---------|-------------|
| `standard` | Docker (default) |
| `singularity` | Singularity container |
| `apptainer` | Apptainer container |
| `conda` | Conda environment |
| `mamba` | Mamba (faster conda) |
| `slurm` | SLURM cluster |
| `sge` | Sun Grid Engine |
| `lsf` | IBM LSF |
| `test` | Test with example data |

## Examples

### Basic local analysis

```bash
nextflow run main.nf \
    --vcf variants.vcf.gz \
    --high_bulk MUT \
    --low_bulk WT
```

### HPC with SLURM and Singularity

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --vcf variants.vcf.gz \
    --high_bulk MUT \
    --low_bulk WT \
    --outdir /scratch/user/bsa_results
```

### With annotation

```bash
nextflow run main.nf \
    --vcf variants.vcf.gz \
    --high_bulk MUT \
    --low_bulk WT \
    --annotate \
    --snpeff_db Sorghum_bicolor
```

### Multiple samples per bulk

```bash
nextflow run main.nf \
    --vcf pooled_calls.vcf.gz \
    --high_bulk "mut_rep1,mut_rep2,mut_rep3" \
    --low_bulk "wt_rep1,wt_rep2,wt_rep3"
```

### Dominant trait with relaxed thresholds

```bash
nextflow run main.nf \
    --vcf variants.vcf.gz \
    --high_bulk affected \
    --low_bulk unaffected \
    --mode dominant \
    --z_threshold 2.5 \
    --min_dp 5
```

## Output

The pipeline produces the following files in the output directory:

```
results/
├── bsa_variants.tsv        # Per-variant statistics
├── bsa_windows.tsv         # Sliding window statistics
├── bsa_regions.tsv         # Candidate regions
├── bsa_regions.bed         # Candidate regions (BED format)
├── bsa_candidates.tsv      # Candidate variants
├── bsa_summary.txt         # Analysis summary
├── bsa_genome_wide.png/pdf # Genome-wide plot
├── bsa_region_*.png/pdf    # Regional plots
└── pipeline_info/
    ├── timeline.html       # Execution timeline
    ├── report.html         # Pipeline report
    ├── trace.txt           # Process trace
    └── dag.svg             # Workflow DAG
```

## Resume failed runs

If a run fails, you can resume from where it left off:

```bash
nextflow run main.nf -resume [other options]
```

## Resource requirements

- Memory: 8 GB minimum (16 GB recommended)
- Time: Depends on VCF size, typically 30 min - 2 hours
- Disk: Output is typically < 100 MB

## Troubleshooting

### Container issues

```bash
# Pull container manually
docker pull username/bsaseq:latest

# Or for Singularity
singularity pull docker://username/bsaseq:latest
```

### Permission errors

For Docker, ensure you're in the docker group:
```bash
sudo usermod -aG docker $USER
# Then log out and back in
```

### Memory errors

Increase memory allocation:
```bash
nextflow run main.nf --memory 32.GB [other options]
```
