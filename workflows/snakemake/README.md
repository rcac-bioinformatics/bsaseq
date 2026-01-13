# bsaseq Snakemake Pipeline

Run bsaseq as a Snakemake workflow for reproducible analysis.

## Quick start

```bash
# Edit config.yaml with your parameters
cp config.yaml my_analysis.yaml
vim my_analysis.yaml

# Run with conda
snakemake --configfile my_analysis.yaml --cores 1 --use-conda

# Or run with existing bsaseq installation
snakemake --configfile my_analysis.yaml --cores 1
```

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) >= 7.0
- One of:
  - bsaseq installed in current environment
  - Conda/Mamba for automatic environment creation

## Configuration

Copy and edit `config.yaml`:

```yaml
# Required
vcf: "path/to/variants.vcf.gz"
high_bulk: "mutant"
low_bulk: "wildtype"

# Optional (defaults shown)
outdir: "results"
window_size: 1000000
z_threshold: 3.0
mode: "recessive"
annotate: false
```

See `config.yaml` for all available parameters.

## Usage

### Basic analysis

```bash
snakemake --configfile config.yaml --cores 1 --use-conda
```

### Dry run (preview)

```bash
snakemake --configfile config.yaml -n
```

### With specific rule

```bash
# Just generate plots
snakemake regenerate_plots --configfile config.yaml --cores 1

# Just list samples
snakemake list_samples --configfile config.yaml --cores 1
```

### On HPC cluster

```bash
# SLURM
snakemake --configfile config.yaml \
    --use-conda \
    --cluster "sbatch -p normal -t 4:00:00 --mem={resources.mem_mb}M" \
    --jobs 10

# SGE
snakemake --configfile config.yaml \
    --use-conda \
    --cluster "qsub -q all.q -l h_rt=4:00:00 -l h_vmem={resources.mem_mb}M" \
    --jobs 10
```

### With Snakemake profiles

Create a profile for your cluster (e.g., `~/.config/snakemake/slurm/config.yaml`):

```yaml
executor: slurm
default-resources:
  - mem_mb=8000
  - runtime=120
jobs: 10
use-conda: true
```

Then run:

```bash
snakemake --configfile config.yaml --profile slurm
```

## Rules

| Rule | Description |
|------|-------------|
| `all` | Default target, runs full analysis |
| `bsaseq_run` | Complete BSA analysis |
| `regenerate_plots` | Regenerate plots from existing output |
| `annotate_candidates` | Annotate variants with snpEff |
| `list_samples` | List sample names in VCF |
| `clean` | Remove all output files |

## Output

```
results/
├── bsa_summary.txt          # Analysis summary
├── bsa_variants.tsv         # Per-variant statistics
├── bsa_windows.tsv          # Sliding window statistics
├── bsa_regions.tsv          # Candidate regions
├── bsa_regions.bed          # Candidate regions (BED)
├── bsa_candidates.tsv       # Candidate variants
├── bsa_genome_wide.png/pdf  # Genome-wide plot
├── bsa_region_*.png/pdf     # Regional plots
├── logs/
│   └── bsaseq.log           # Analysis log
├── benchmarks/
│   └── bsaseq.txt           # Runtime benchmarks
└── samples.txt              # Sample list (if generated)
```

## Examples

### Multiple samples per bulk

```yaml
vcf: "pooled_calls.vcf.gz"
high_bulk: "mut1,mut2,mut3"
low_bulk: "wt1,wt2,wt3"
```

### With annotation

```yaml
vcf: "variants.vcf.gz"
high_bulk: "mut"
low_bulk: "wt"
annotate: true
snpeff_db: "Arabidopsis_thaliana"
```

### Dominant trait

```yaml
vcf: "variants.vcf.gz"
high_bulk: "affected"
low_bulk: "unaffected"
mode: "dominant"
z_threshold: 2.5
```

## Troubleshooting

### Conda environment creation fails

```bash
# Create environment manually
conda env create -f environment.yaml
conda activate bsaseq

# Then run without --use-conda
snakemake --configfile config.yaml --cores 1
```

### Out of memory

Increase memory in config or command line:

```bash
snakemake --configfile config.yaml \
    --resources mem_mb=16000 \
    --cores 1
```

### Resume after failure

Snakemake automatically resumes from where it left off:

```bash
snakemake --configfile config.yaml --cores 1
```
