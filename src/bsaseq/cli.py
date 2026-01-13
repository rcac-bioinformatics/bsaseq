"""Command-line interface for bsaseq.

This module defines the Click-based CLI for the bsaseq package,
providing commands for bulk segregant analysis of VCF data.
"""

from __future__ import annotations

from pathlib import Path

import click
from rich.console import Console
from rich.table import Table

from bsaseq import __version__
from bsaseq.analysis.candidates import (
    AnnotatedCandidate,
    CandidateGene,
    InheritanceMode,
    annotate_candidates,
    filter_candidate_variants,
    get_variant_filter_thresholds,
    identify_candidate_regions,
    summarize_candidate_genes,
)
from bsaseq.analysis.windows import calculate_windows
from bsaseq.io.summary import AnalysisSummary, write_summary
from bsaseq.io.vcf import get_sample_names, parse_vcf, write_candidates_vcf
from bsaseq.io.writers import (
    write_annotated_candidates_tsv,
    write_candidate_genes_tsv,
    write_candidate_variants_tsv,
    write_regions_bed,
    write_regions_tsv,
    write_variants_tsv,
    write_windows_tsv,
)
from bsaseq.utils.logging import (
    print_error,
    print_info,
    print_stats,
    print_success,
    print_warning,
    setup_logging,
)

console = Console(stderr=True)

# Context settings for all commands
CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__, prog_name="bsaseq")
@click.option(
    "--verbose", "-v", is_flag=True, default=False, help="Enable verbose output"
)
@click.pass_context
def cli(ctx: click.Context, verbose: bool) -> None:
    """bsaseq: Bulk Segregant Analysis for QTL mapping.

    Identify genomic loci controlling traits by comparing allele frequencies
    between phenotypically distinct pooled DNA samples.

    \b
    Quick start:
        bsaseq run --vcf calls.vcf.gz --high-bulk mutant --low-bulk wildtype -o results

    \b
    Common workflows:
        bsaseq samples --vcf data.vcf.gz          # List sample names
        bsaseq run --vcf data.vcf.gz ...          # Full analysis
        bsaseq plot --windows windows.tsv ...      # Regenerate plots
        bsaseq annotate --candidates cand.tsv ...  # Add annotations

    \b
    Documentation:
        https://github.com/username/bsaseq
    """
    import logging

    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose

    level = logging.DEBUG if verbose else logging.INFO
    setup_logging(level=level)


@cli.command()
@click.option(
    "--vcf",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    metavar="FILE",
    help="Input VCF file (bgzipped recommended).",
)
def samples(vcf: Path) -> None:
    """List sample names in a VCF file.

    Displays all sample names present in a VCF file. Use this to identify
    the correct sample names before running the analysis.

    \b
    Example:
        bsaseq samples --vcf joint_calls.vcf.gz

    Sample names are case-sensitive and must match exactly when used
    with the --high-bulk and --low-bulk options.
    """
    try:
        sample_names = get_sample_names(vcf)

        if not sample_names:
            print_error("No samples found in VCF file")
            raise SystemExit(1)

        console.print(f"\n[bold]Samples in {vcf.name}:[/bold]\n")

        table = Table(show_header=True, header_style="bold cyan")
        table.add_column("Index", justify="right", style="dim")
        table.add_column("Sample Name")

        for i, name in enumerate(sample_names):
            table.add_row(str(i), name)

        console.print(table)
        console.print(f"\n[dim]Total: {len(sample_names)} sample(s)[/dim]\n")

    except Exception as e:
        print_error(f"Failed to read VCF: {e}")
        raise SystemExit(1) from e


@cli.command()
@click.option(
    "--vcf",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    metavar="FILE",
    help="Input VCF file (bgzipped recommended).",
)
@click.option(
    "--high-bulk",
    required=True,
    metavar="SAMPLES",
    help='Sample name(s) for high/mutant bulk. Comma-separated for multiple: "mut1,mut2".',
)
@click.option(
    "--low-bulk",
    required=True,
    metavar="SAMPLES",
    help='Sample name(s) for low/wild-type bulk. Comma-separated for multiple: "wt1,wt2".',
)
@click.option(
    "--out",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    metavar="PREFIX",
    help="Output prefix. Files will be named {PREFIX}_variants.tsv, etc.",
)
@click.option(
    "--window-size",
    default=1_000_000,
    show_default=True,
    type=int,
    metavar="BP",
    help="Sliding window size in base pairs.",
)
@click.option(
    "--step-size",
    default=250_000,
    show_default=True,
    type=int,
    metavar="BP",
    help="Step size between windows in base pairs.",
)
@click.option(
    "--min-variants",
    default=5,
    show_default=True,
    type=int,
    metavar="INT",
    help="Minimum variants per window. Windows with fewer are skipped.",
)
@click.option(
    "--min-dp",
    default=10,
    show_default=True,
    type=int,
    metavar="INT",
    help="Minimum read depth. Variants below this in either bulk are excluded.",
)
@click.option(
    "--max-dp",
    default=200,
    show_default=True,
    type=int,
    metavar="INT",
    help="Maximum read depth. Excludes variants in repetitive regions.",
)
@click.option(
    "--min-gq",
    default=20,
    show_default=True,
    type=int,
    metavar="INT",
    help="Minimum genotype quality (GQ field).",
)
@click.option(
    "--min-qual",
    default=30.0,
    show_default=True,
    type=float,
    metavar="FLOAT",
    help="Minimum variant quality (QUAL field).",
)
@click.option(
    "--z-threshold",
    default=3.0,
    show_default=True,
    type=float,
    metavar="FLOAT",
    help="Z-score threshold for candidate regions. Higher = more stringent.",
)
@click.option(
    "--mode",
    type=click.Choice(["recessive", "dominant"]),
    default="recessive",
    show_default=True,
    help="Expected inheritance mode. Affects variant filtering thresholds.",
)
@click.option(
    "--plot/--no-plot",
    default=True,
    show_default=True,
    help="Generate visualization plots.",
)
@click.option(
    "--plot-format",
    type=click.Choice(["png", "pdf", "both"]),
    default="both",
    show_default=True,
    help="Output format for plots.",
)
@click.option(
    "--max-region-plots",
    default=5,
    show_default=True,
    type=int,
    metavar="INT",
    help="Maximum number of regional zoom plots to generate.",
)
@click.option(
    "--annotate/--no-annotate",
    default=False,
    show_default=True,
    help="Run snpEff annotation on candidate variants.",
)
@click.option(
    "--snpeff-db",
    default=None,
    metavar="DATABASE",
    help="snpEff database name (e.g., Sorghum_bicolor). Required if --annotate.",
)
@click.option(
    "--snpeff-mem",
    default="4g",
    show_default=True,
    metavar="MEM",
    help="Java heap memory for snpEff (e.g., 4g, 8g).",
)
@click.pass_context
def run(
    ctx: click.Context,
    vcf: Path,
    high_bulk: str,
    low_bulk: str,
    out: Path,
    window_size: int,
    step_size: int,
    min_variants: int,
    min_dp: int,
    max_dp: int,
    min_gq: int,
    min_qual: float,
    z_threshold: float,
    mode: str,
    plot: bool,
    plot_format: str,
    max_region_plots: int,
    annotate: bool,
    snpeff_db: str | None,
    snpeff_mem: str,
) -> None:
    """Run complete BSA analysis pipeline.

    Parses the input VCF file, calculates allele frequency differences
    between the high (mutant) and low (wild-type) bulk samples, performs
    sliding window analysis, identifies candidate regions, and filters
    potential causal variants.

    \b
    Examples:
      Basic analysis:
        bsaseq run --vcf calls.vcf.gz --high-bulk MUT --low-bulk WT -o results

      Multiple samples per bulk:
        bsaseq run --vcf calls.vcf.gz --high-bulk "MUT1,MUT2" --low-bulk "WT1,WT2" -o results

      Dominant trait with annotation:
        bsaseq run --vcf calls.vcf.gz --high-bulk MUT --low-bulk WT -o results \\
            --mode dominant --annotate --snpeff-db Sorghum_bicolor

      Stringent analysis (fewer false positives):
        bsaseq run --vcf calls.vcf.gz --high-bulk MUT --low-bulk WT -o results \\
            --z-threshold 4.0 --min-dp 20 --min-variants 10

    \b
    Output files:
      {PREFIX}_variants.tsv           Per-variant allele frequencies
      {PREFIX}_windows.tsv            Sliding window statistics
      {PREFIX}_regions.tsv            Candidate genomic regions
      {PREFIX}_regions.bed            Regions in BED format
      {PREFIX}_candidates.tsv         Filtered candidate variants
      {PREFIX}_summary.txt            Analysis summary
      {PREFIX}_genome_wide.png/pdf    Manhattan plot (if --plot)
      {PREFIX}_region_*.png/pdf       Regional zoom plots (if --plot)
      {PREFIX}_af_distribution.png    AF diagnostics (if --plot)
      {PREFIX}_depth_distribution.png Depth diagnostics (if --plot)
      {PREFIX}_annotated_candidates.tsv  Annotated variants (if --annotate)
      {PREFIX}_candidate_genes.tsv    Gene summary (if --annotate)
    """
    import numpy as np

    # Parse sample names
    high_samples = [s.strip() for s in high_bulk.split(",")]
    low_samples = [s.strip() for s in low_bulk.split(",")]

    print_info("Starting BSA analysis")
    print_info(f"Input VCF: {vcf}")
    print_info(f"High bulk: {', '.join(high_samples)}")
    print_info(f"Low bulk: {', '.join(low_samples)}")
    print_info(f"Output prefix: {out}")
    print_info(f"Window: {window_size:,} bp, step: {step_size:,} bp")
    print_info(f"Z-threshold: {z_threshold}, Mode: {mode}")

    # Store parameters for summary
    parameters = {
        "window_size": window_size,
        "step_size": step_size,
        "min_variants": min_variants,
        "min_dp": min_dp,
        "max_dp": max_dp,
        "min_gq": min_gq,
        "min_qual": min_qual,
        "z_threshold": z_threshold,
        "mode": mode,
    }

    try:
        # Parse VCF and collect variants
        variants = list(
            parse_vcf(
                vcf_path=vcf,
                high_bulk=high_bulk,
                low_bulk=low_bulk,
                min_dp=min_dp,
                max_dp=max_dp,
                min_gq=min_gq,
                min_qual=min_qual,
            )
        )

        if not variants:
            print_error("No variants passed filters")
            raise SystemExit(1)

        # Calculate summary statistics
        delta_afs = np.array([v.delta_af for v in variants])
        af_highs = np.array([v.af_high for v in variants])
        af_lows = np.array([v.af_low for v in variants])

        # Collect chromosome counts
        chrom_counts: dict[str, int] = {}
        for v in variants:
            chrom_counts[v.chrom] = chrom_counts.get(v.chrom, 0) + 1

        # Print variant summary statistics
        variant_stats = {
            "Total variants": len(variants),
            "Chromosomes": len(chrom_counts),
            "Mean delta-AF": float(np.mean(delta_afs)),
            "Median delta-AF": float(np.median(delta_afs)),
            "Std delta-AF": float(np.std(delta_afs)),
            "Min delta-AF": float(np.min(delta_afs)),
            "Max delta-AF": float(np.max(delta_afs)),
            "Mean AF (high bulk)": float(np.mean(af_highs)),
            "Mean AF (low bulk)": float(np.mean(af_lows)),
        }

        print_stats(variant_stats, title="Variant Summary")

        # Print per-chromosome counts
        console.print("\n[bold]Variants per chromosome:[/bold]")
        chrom_table = Table(show_header=True, header_style="bold cyan")
        chrom_table.add_column("Chromosome")
        chrom_table.add_column("Variants", justify="right")

        for chrom in sorted(chrom_counts.keys()):
            chrom_table.add_row(chrom, f"{chrom_counts[chrom]:,}")

        console.print(chrom_table)

        # Create output directory if needed
        out_dir = out.parent
        if out_dir and str(out_dir) != "." and not out_dir.exists():
            out_dir.mkdir(parents=True, exist_ok=True)

        # Write variant data to TSV file
        out_variants = Path(f"{out}_variants.tsv")
        write_variants_tsv(variants, out_variants)

        # Calculate sliding windows
        print_info("Calculating sliding windows...")
        windows = list(
            calculate_windows(
                variants=variants,
                window_size=window_size,
                step_size=step_size,
                min_variants=min_variants,
            )
        )

        if not windows:
            print_error(
                f"No windows generated. Try reducing --min-variants "
                f"(currently {min_variants}) or increasing --window-size."
            )
            raise SystemExit(1)

        # Print window summary statistics
        window_delta_afs = np.array([w.mean_delta_af for w in windows])
        window_z_scores = np.array(
            [w.z_score for w in windows if w.z_score is not None]
        )

        window_stats = {
            "Total windows": len(windows),
            "Mean window delta-AF": float(np.mean(window_delta_afs)),
            "Max window delta-AF": float(np.max(window_delta_afs)),
            "Min window delta-AF": float(np.min(window_delta_afs)),
        }

        significant_windows = sum(
            1 for w in windows if w.z_score is not None and w.z_score >= z_threshold
        )

        if len(window_z_scores) > 0:
            window_stats["Max Z-score"] = float(np.max(window_z_scores))
            window_stats["Min Z-score"] = float(np.min(window_z_scores))
            window_stats[f"Windows >= z={z_threshold}"] = significant_windows

        print_stats(window_stats, title="Window Summary")

        # Write window data to TSV file
        out_windows = Path(f"{out}_windows.tsv")
        write_windows_tsv(windows, out_windows)

        # Identify candidate regions
        print_info("Identifying candidate regions...")
        regions = identify_candidate_regions(
            windows=windows,
            z_threshold=z_threshold,
            min_consecutive=2,
            merge_distance=500_000,
        )

        # Write regions
        out_regions_tsv = Path(f"{out}_regions.tsv")
        out_regions_bed = Path(f"{out}_regions.bed")
        write_regions_tsv(regions, out_regions_tsv)
        write_regions_bed(regions, out_regions_bed)

        if regions:
            region_stats = {
                "Candidate regions": len(regions),
                "Top region chromosome": regions[0].chrom,
                "Top region start": regions[0].start,
                "Top region end": regions[0].end,
                "Top region Z-score": regions[0].max_z_score,
            }
            print_stats(region_stats, title="Candidate Regions")

        # Filter candidate variants
        print_info("Filtering candidate variants...")
        inheritance_mode = InheritanceMode(mode)
        thresholds = get_variant_filter_thresholds(inheritance_mode)

        candidate_variants = filter_candidate_variants(
            variants=variants,
            regions=regions,
            min_delta_af=thresholds["min_delta_af"],
            min_af_high=thresholds["min_af_high"],
            max_af_low=thresholds["max_af_low"],
        )

        # Write candidate variants
        out_candidates = Path(f"{out}_candidates.tsv")
        write_candidate_variants_tsv(candidate_variants, out_candidates)

        if candidate_variants:
            top_candidate = candidate_variants[0]
            candidate_stats = {
                "Candidate variants": len(candidate_variants),
                "Top candidate": f"{top_candidate.chrom}:{top_candidate.pos}",
                "Top candidate delta-AF": top_candidate.delta_af,
                "Top candidate AF (high)": top_candidate.af_high,
                "Top candidate AF (low)": top_candidate.af_low,
            }
            print_stats(candidate_stats, title="Candidate Variants")
        else:
            print_info("No candidate causal variants found with current thresholds")

        # Annotation step
        annotated_candidates: list[AnnotatedCandidate] = []
        candidate_genes: list[CandidateGene] = []
        out_annotated: Path | None = None
        out_genes: Path | None = None

        if annotate:
            if not snpeff_db:
                print_error("--snpeff-db is required when using --annotate")
                raise SystemExit(1)

            # Import annotation functions
            from bsaseq.annotation import (
                check_snpeff_available,
                get_best_annotation,
                parse_snpeff_vcf,
                run_snpeff,
            )

            if not check_snpeff_available():
                print_warning("snpEff not found in PATH, skipping annotation")
                print_warning("Install with: conda install -c bioconda snpeff")
            elif candidate_variants:
                print_info(f"Annotating candidates with snpEff database: {snpeff_db}")

                try:
                    # Write candidates to VCF for snpEff
                    out_candidates_vcf = Path(f"{out}_candidates.vcf")
                    write_candidates_vcf(candidate_variants, out_candidates_vcf)

                    # Run snpEff
                    annotated_vcf = run_snpeff(
                        vcf_path=out_candidates_vcf,
                        database=snpeff_db,
                        java_mem=snpeff_mem,
                    )

                    # Parse annotations and build lookup
                    annotations: dict[tuple[str, int, str, str], object] = {}
                    for ann in parse_snpeff_vcf(annotated_vcf):
                        key = (ann.chrom, ann.pos, ann.ref, ann.alt)
                        existing = annotations.get(key)
                        if existing is None:
                            annotations[key] = ann
                        else:
                            # Keep the best annotation
                            best = get_best_annotation([existing, ann])
                            if best:
                                annotations[key] = best

                    # Annotate candidates
                    annotated_candidates = annotate_candidates(
                        candidates=candidate_variants,
                        annotations=annotations,
                    )

                    # Summarize by gene
                    candidate_genes = summarize_candidate_genes(annotated_candidates)

                    # Write output files
                    out_annotated = Path(f"{out}_annotated_candidates.tsv")
                    out_genes = Path(f"{out}_candidate_genes.tsv")
                    write_annotated_candidates_tsv(annotated_candidates, out_annotated)
                    write_candidate_genes_tsv(candidate_genes, out_genes)

                    # Print annotation summary
                    n_annotated = sum(1 for c in annotated_candidates if c.gene_name)
                    n_lof = sum(1 for c in annotated_candidates if c.is_lof)
                    n_high = sum(1 for c in annotated_candidates if c.impact == "HIGH")

                    annotation_stats = {
                        "Annotated variants": n_annotated,
                        "Candidate genes": len(candidate_genes),
                        "Loss-of-function": n_lof,
                        "HIGH impact": n_high,
                    }
                    if candidate_genes:
                        annotation_stats["Top gene"] = candidate_genes[0].gene_name
                    print_stats(annotation_stats, title="Annotation Summary")

                except Exception as e:
                    print_warning(f"Annotation failed: {e}")
                    if ctx.obj.get("verbose"):
                        console.print_exception()
            else:
                print_info("No candidate variants to annotate")

        # Create and write summary
        # Calculate annotation stats
        n_annotated = sum(1 for c in annotated_candidates if c.gene_name) if annotated_candidates else 0
        n_lof = sum(1 for c in annotated_candidates if c.is_lof) if annotated_candidates else 0
        top_gene_name = candidate_genes[0].gene_name if candidate_genes else None

        summary = AnalysisSummary(
            vcf_path=str(vcf),
            high_bulk_samples=high_samples,
            low_bulk_samples=low_samples,
            total_variants=len(variants),  # We don't have pre-filter count easily
            passing_variants=len(variants),
            total_windows=len(windows),
            significant_windows=significant_windows,
            candidate_regions=len(regions),
            candidate_variants=len(candidate_variants),
            top_region=regions[0] if regions else None,
            parameters=parameters,
            annotated_candidates=n_annotated,
            candidate_genes=len(candidate_genes),
            lof_variants=n_lof,
            top_gene=top_gene_name,
        )

        out_summary = Path(f"{out}_summary.txt")
        write_summary(summary, out_summary)

        # Generate plots if requested
        plot_files: list[Path] = []
        if plot:
            print_info("Generating plots...")
            try:
                import matplotlib.pyplot as plt

                from bsaseq.plotting.diagnostics import (
                    plot_af_distribution,
                    plot_depth_distribution,
                )
                from bsaseq.plotting.genome import (
                    plot_all_regions,
                    plot_genome_wide,
                )

                # Genome-wide Manhattan plot
                out_genome = Path(f"{out}_genome_wide")
                fig = plot_genome_wide(
                    windows=windows,
                    regions=regions,
                    z_threshold=z_threshold,
                    output_path=out_genome,
                    title="BSA Genome-wide Analysis",
                )
                plt.close(fig)
                plot_files.extend([out_genome.with_suffix(".png"), out_genome.with_suffix(".pdf")])

                # Regional zoom plots
                if regions and max_region_plots > 0:
                    region_paths = plot_all_regions(
                        variants=variants,
                        regions=regions,
                        windows=windows,
                        candidate_variants=candidate_variants,
                        output_dir=out.parent if out.parent.name else Path("."),
                        max_regions=max_region_plots,
                    )
                    # Rename files to include output prefix
                    for p in region_paths:
                        if p.exists():
                            # Extract region info from filename
                            new_name = f"{out.name}_{p.name}"
                            new_path = p.parent / new_name
                            p.rename(new_path)
                            plot_files.append(new_path)

                # AF distribution plot
                out_af_dist = Path(f"{out}_af_distribution")
                fig = plot_af_distribution(variants=variants, output_path=out_af_dist)
                plt.close(fig)
                plot_files.append(out_af_dist.with_suffix(".png"))

                # Depth distribution plot
                out_depth_dist = Path(f"{out}_depth_distribution")
                fig = plot_depth_distribution(variants=variants, output_path=out_depth_dist)
                plt.close(fig)
                plot_files.append(out_depth_dist.with_suffix(".png"))

                print_success(f"Generated {len(plot_files)} plot file(s)")

            except ImportError as e:
                print_warning(f"Plotting skipped: {e}")
            except Exception as e:
                print_warning(f"Plotting failed: {e}")
                if ctx.obj.get("verbose"):
                    console.print_exception()

        print_success("Analysis complete!")
        print_info("Output files:")
        print_info(f"  Variants: {out_variants}")
        print_info(f"  Windows: {out_windows}")
        print_info(f"  Regions (TSV): {out_regions_tsv}")
        print_info(f"  Regions (BED): {out_regions_bed}")
        print_info(f"  Candidates: {out_candidates}")
        print_info(f"  Summary: {out_summary}")
        if out_annotated:
            print_info(f"  Annotated: {out_annotated}")
        if out_genes:
            print_info(f"  Genes: {out_genes}")
        if plot_files:
            print_info("  Plots:")
            for pf in plot_files[:10]:  # Show first 10
                if pf.exists():
                    print_info(f"    {pf}")
            if len(plot_files) > 10:
                print_info(f"    ... and {len(plot_files) - 10} more")

    except ValueError as e:
        print_error(str(e))
        raise SystemExit(1) from e
    except Exception as e:
        print_error(f"Analysis failed: {e}")
        if ctx.obj.get("verbose"):
            console.print_exception()
        raise SystemExit(1) from e


@cli.command()
@click.option(
    "--windows",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Windows TSV file from previous run",
)
@click.option(
    "--regions",
    type=click.Path(exists=True, path_type=Path),
    help="Regions TSV file from previous run",
)
@click.option(
    "--variants",
    type=click.Path(exists=True, path_type=Path),
    help="Variants TSV file from previous run (for diagnostic plots)",
)
@click.option(
    "--out",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix for plot files",
)
@click.option(
    "--z-threshold",
    default=3.0,
    show_default=True,
    type=float,
    help="Z-score threshold line for genome plot",
)
@click.option(
    "--format",
    "fmt",
    type=click.Choice(["png", "pdf", "both"]),
    default="both",
    show_default=True,
    help="Output format for plots",
)
@click.pass_context
def plot(
    ctx: click.Context,
    windows: Path,
    regions: Path | None,
    variants: Path | None,
    out: Path,
    z_threshold: float,
    fmt: str,
) -> None:
    """Generate plots from existing analysis output.

    Re-generate plots without re-running the full analysis pipeline.
    Useful for adjusting plot parameters or generating additional plots.

    Example:

        bsaseq plot --windows analysis_windows.tsv --regions analysis_regions.tsv -o new_plots
    """
    import matplotlib.pyplot as plt

    from bsaseq.io.readers import read_regions_tsv, read_variants_tsv, read_windows_tsv
    from bsaseq.plotting.diagnostics import plot_af_distribution, plot_depth_distribution
    from bsaseq.plotting.genome import plot_genome_wide

    print_info("Generating plots from existing data...")

    try:
        # Read windows
        window_list = read_windows_tsv(windows)
        if not window_list:
            print_error("No windows found in input file")
            raise SystemExit(1)

        # Read regions if provided
        region_list = []
        if regions:
            region_list = read_regions_tsv(regions)

        # Read variants if provided
        variant_list = []
        if variants:
            variant_list = read_variants_tsv(variants)

        # Create output directory if needed
        out_dir = out.parent
        if out_dir and str(out_dir) != "." and not out_dir.exists():
            out_dir.mkdir(parents=True, exist_ok=True)

        plot_count = 0

        # Genome-wide plot
        print_info("Generating genome-wide plot...")
        out_genome = Path(f"{out}_genome_wide")
        fig = plot_genome_wide(
            windows=window_list,
            regions=region_list if region_list else None,
            z_threshold=z_threshold,
            output_path=out_genome,
        )
        plt.close(fig)
        plot_count += 2  # PNG + PDF

        # Diagnostic plots (if variants provided)
        if variant_list:
            print_info("Generating diagnostic plots...")
            out_af = Path(f"{out}_af_distribution")
            fig = plot_af_distribution(variants=variant_list, output_path=out_af)
            plt.close(fig)
            plot_count += 1

            out_depth = Path(f"{out}_depth_distribution")
            fig = plot_depth_distribution(variants=variant_list, output_path=out_depth)
            plt.close(fig)
            plot_count += 1

        print_success(f"Generated {plot_count} plot file(s)")
        print_info(f"Output: {out}_*.png/pdf")

    except FileNotFoundError as e:
        print_error(str(e))
        raise SystemExit(1) from e
    except Exception as e:
        print_error(f"Plot generation failed: {e}")
        if ctx.obj.get("verbose"):
            console.print_exception()
        raise SystemExit(1) from e


@cli.command("check-snpeff")
@click.pass_context
def check_snpeff(ctx: click.Context) -> None:
    """Check if snpEff is installed and show version.

    Verifies that snpEff is available in PATH and displays version
    information. This is useful for troubleshooting annotation issues.

    Example:

        bsaseq check-snpeff
    """
    from bsaseq.annotation import (
        check_snpeff_available,
        get_snpeff_version,
        list_snpeff_databases,
    )

    if not check_snpeff_available():
        print_error("snpEff not found in PATH")
        print_info("Install with: conda install -c bioconda snpeff")
        raise SystemExit(1)

    version = get_snpeff_version()
    print_success("snpEff is installed")
    if version:
        print_info(f"Version: {version}")

    print_info("Checking available databases (this may take a moment)...")
    databases = list_snpeff_databases()
    if databases:
        console.print(f"\n[bold]Available databases ({len(databases)} total):[/bold]")
        # Show first 10 as examples
        for db in databases[:10]:
            console.print(f"  {db}")
        if len(databases) > 10:
            console.print(f"  ... and {len(databases) - 10} more")
        console.print("\n[dim]Use 'snpEff databases' for full list[/dim]")


@cli.command()
@click.option(
    "--candidates",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Candidate variants TSV file from previous run",
)
@click.option(
    "--snpeff-db",
    required=True,
    help="snpEff database name (e.g., 'Sorghum_bicolor')",
)
@click.option(
    "--out",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output prefix for annotated files",
)
@click.option(
    "--snpeff-mem",
    default="4g",
    show_default=True,
    help="Java heap memory for snpEff (e.g., '4g', '8g')",
)
@click.pass_context
def annotate(
    ctx: click.Context,
    candidates: Path,
    snpeff_db: str,
    out: Path,
    snpeff_mem: str,
) -> None:
    """Annotate candidate variants using snpEff.

    Takes candidate variants from a previous BSA analysis run and
    annotates them using snpEff. This is useful for re-annotating
    with a different database or re-running annotation after updating
    snpEff.

    Example:

        bsaseq annotate --candidates analysis_candidates.tsv \\
            --snpeff-db Sorghum_bicolor -o annotated

    Output files:

        {out}_annotated_candidates.tsv - Annotated variant details

        {out}_candidate_genes.tsv - Gene-level summary
    """
    from bsaseq.annotation import (
        check_snpeff_available,
        get_best_annotation,
        parse_snpeff_vcf,
        run_snpeff,
    )
    from bsaseq.io.readers import read_candidate_variants_tsv

    print_info("Annotating candidates with snpEff...")
    print_info(f"Input: {candidates}")
    print_info(f"Database: {snpeff_db}")

    try:
        # Check snpEff availability
        if not check_snpeff_available():
            print_error("snpEff not found in PATH")
            print_info("Install with: conda install -c bioconda snpeff")
            raise SystemExit(1)

        # Read candidate variants
        candidate_list = read_candidate_variants_tsv(candidates)
        if not candidate_list:
            print_error("No candidates found in input file")
            raise SystemExit(1)

        print_info(f"Loaded {len(candidate_list)} candidate variants")

        # Create output directory if needed
        out_dir = out.parent
        if out_dir and str(out_dir) != "." and not out_dir.exists():
            out_dir.mkdir(parents=True, exist_ok=True)

        # Write candidates to VCF
        out_vcf = Path(f"{out}_candidates.vcf")
        write_candidates_vcf(candidate_list, out_vcf)

        # Run snpEff
        annotated_vcf = run_snpeff(
            vcf_path=out_vcf,
            database=snpeff_db,
            java_mem=snpeff_mem,
        )

        # Parse annotations
        annotations: dict[tuple[str, int, str, str], object] = {}
        for ann in parse_snpeff_vcf(annotated_vcf):
            key = (ann.chrom, ann.pos, ann.ref, ann.alt)
            existing = annotations.get(key)
            if existing is None:
                annotations[key] = ann
            else:
                best = get_best_annotation([existing, ann])
                if best:
                    annotations[key] = best

        # Annotate candidates
        annotated_list = annotate_candidates(
            candidates=candidate_list,
            annotations=annotations,
        )

        # Summarize by gene
        gene_list = summarize_candidate_genes(annotated_list)

        # Write outputs
        out_annotated = Path(f"{out}_annotated_candidates.tsv")
        out_genes = Path(f"{out}_candidate_genes.tsv")
        write_annotated_candidates_tsv(annotated_list, out_annotated)
        write_candidate_genes_tsv(gene_list, out_genes)

        # Print summary
        n_annotated = sum(1 for c in annotated_list if c.gene_name)
        n_lof = sum(1 for c in annotated_list if c.is_lof)

        annotation_stats = {
            "Total candidates": len(candidate_list),
            "Annotated": n_annotated,
            "Candidate genes": len(gene_list),
            "Loss-of-function": n_lof,
        }
        if gene_list:
            annotation_stats["Top gene"] = gene_list[0].gene_name
        print_stats(annotation_stats, title="Annotation Results")

        print_success("Annotation complete!")
        print_info("Output files:")
        print_info(f"  Annotated candidates: {out_annotated}")
        print_info(f"  Candidate genes: {out_genes}")

    except FileNotFoundError as e:
        print_error(str(e))
        raise SystemExit(1) from e
    except Exception as e:
        print_error(f"Annotation failed: {e}")
        if ctx.obj.get("verbose"):
            console.print_exception()
        raise SystemExit(1) from e


def main() -> None:
    """Main entry point for CLI."""
    cli()


if __name__ == "__main__":
    main()
