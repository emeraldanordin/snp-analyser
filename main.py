#!/usr/bin/env python3
"""
main.py
-------
Command-line interface for the SNP Annotation Analyser.

SINGLE-VCF MODE (one strain/isolate, or multiple isolates within one
comparison):
    python main.py --vcf data/your_file_name.vcf --output output/
    python main.py --vcf your_file_name1.vcf your_file_name2.vcf your_file_name3.vcf --output output/
    python main.py --vcf data/your_file_name.vcf --output output/ --min-qual 30 --min-dp 10 --no-pdf (no PDF output -- OPTIONAL)


TWO-STRAIN COMPARISON MODE:
    python main.py --vcf1 data/your_file_name1.vcf --vcf2 data/your_file_name2.vcf --name1 strain1 --name2 strain2 --output output/

Produces (per strain, under output/<strain_name>/):
    <strain>_SNPs.csv / .xlsx
    <strain>_gene_summary.csv
    <strain>_position_summary.csv
    <strain>_report.html / .pdf

Plus (under output/comparison/):
    <strain1>_<strain2>_shared_SNPs.csv
    <strain1>_<strain2>_<strain1>_specific.csv
    <strain1>_<strain2>_<strain2>_specific.csv
    <strain1>_<strain2>_gene_comparison.csv
    <strain1>_<strain2>_position_comparison.csv
    <strain1>_<strain2>_comparison.html / .pdf

MULTI-STRAIN COMPARISON MODE (2 to 10+ strains):
    python main.py --vcfs data/your_file_name1.vcf data/your_file_name2.vcf data/your_file_name3.vcf --strain-names strain1 strain2 strain3 --output output/

If --strain-names is omitted, names are derived from filenames. Produces
(per strain, under output/<strain_name>/): the same per-strain files as
two-strain mode. Plus (under output/multi_strain/):
    <prefix>_snp_matrix.csv               -- full N-way position matrix
    <prefix>_core_genome_positions.csv    -- positions shared by ALL strains
    <prefix>_<strain>_unique.csv          -- one file per strain
    <prefix>_gene_comparison.csv
    <prefix>_pairwise_sharing_matrix.csv
    <prefix>_comparison.html / .pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from comparison import (
    build_gene_comparison,
    build_gene_comparison_n_way,
    build_position_comparison,
    build_snp_matrix,
    classify_snp_matrix,
    comparison_summary,
    comparison_summary_n_way,
    pairwise_sharing_matrix_n_way,
    snp_sharing_matrix,
)
from export import export_all, export_comparison_all, export_multi_strain_comparison, export_strain_labeled
from grouping import group_by_gene, group_by_position
from pipeline import (
    apply_quality_filters,
    process_multiple_strains,
    process_multiple_vcfs,
    process_two_strains,
    process_vcf,
)
from report import generate_comparison_reports, generate_multi_strain_reports, generate_reports
from statistics import compute_filter_impact, compute_summary_statistics
from utils import ensure_output_dir, ensure_strain_output_dirs, get_logger
from validation import check_reference_compatibility
from visualisations import build_all_figures, build_comparison_figures, build_multi_strain_figures
from vcf_parser import VCFValidationError

logger = get_logger(__name__)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="SNP Annotation Analyser -- local, annotated-VCF SNP analysis.",
    )
    parser.add_argument(
        "--vcf", nargs="+", default=None,
        help="Single-VCF mode: path(s) to one or more .vcf / .vcf.gz files. "
             "Multiple files here are treated as isolates of the SAME strain "
             "and compared as isolates, not as separate strains -- use "
             "--vcf1/--vcf2 or --vcfs for strain comparison instead.",
    )
    parser.add_argument("--vcf1", default=None, help="Two-strain mode: path to strain 1's VCF")
    parser.add_argument("--vcf2", default=None, help="Two-strain mode: path to strain 2's VCF")
    parser.add_argument("--name1", default=None, help="Two-strain mode: display name for strain 1 (default: derived from filename)")
    parser.add_argument("--name2", default=None, help="Two-strain mode: display name for strain 2 (default: derived from filename)")
    parser.add_argument(
        "--vcfs", nargs="+", default=None,
        help="Multi-strain mode: paths to 2 or more VCF files (any number, "
             "including 10+). Use with --strain-names.",
    )
    parser.add_argument(
        "--strain-names", nargs="+", default=None,
        help="Multi-strain mode: one name per file in --vcfs, same order. "
             "If omitted, names are derived from filenames. Must be unique.",
    )
    parser.add_argument(
        "--output", default="output", help="Output directory (default: output/)",
    )
    parser.add_argument("--min-qual", type=float, default=None, help="Minimum QUAL filter (optional)")
    parser.add_argument("--min-dp", type=int, default=None, help="Minimum DP filter (optional)")
    parser.add_argument("--no-pdf", action="store_true", help="Skip PDF report generation")
    return parser


def _strain_name_from_path(path: str) -> str:
    """Derive a default strain name from a VCF filename if --name1/--name2
    aren't given, e.g. 'data/your_file_name.vcf' -> 'strain1'."""
    stem = Path(path).name
    for suffix in (".vcf.gz", ".vcf"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem


def run_single_vcf_mode(args: argparse.Namespace) -> int:
    """Original single-VCF (or same-strain multi-isolate) analysis mode."""
    output_dir = ensure_output_dir(args.output)
    logger.info("Starting SNP Annotation Analyser CLI run (single-VCF mode)")
    logger.info("Input VCF(s): %s", args.vcf)

    try:
        if len(args.vcf) == 1:
            df, parse_report = process_vcf(args.vcf[0])
            parse_reports = [parse_report]
        else:
            df, parse_reports = process_multiple_vcfs(args.vcf)
            parse_report = parse_reports[0]
    except VCFValidationError as exc:
        logger.error("VCF validation failed: %s", exc)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        logger.exception("Unexpected error while parsing VCF(s)")
        print(f"ERROR: Unexpected failure while parsing VCF(s): {exc}", file=sys.stderr)
        return 1

    if df.empty:
        print("ERROR: No SNP records were extracted from the given VCF(s).", file=sys.stderr)
        return 1

    records_before = len(df)
    if args.min_qual is not None or args.min_dp is not None:
        df = apply_quality_filters(df, min_qual=args.min_qual, min_dp=args.min_dp)
        impact = compute_filter_impact(records_before, len(df))
        logger.info(
            "Applied QC filters (min_qual=%s, min_dp=%s): %d -> %d records (%.1f%% retained)",
            args.min_qual, args.min_dp, impact["records_before"], impact["records_after"],
            impact["pct_retained"],
        )

    gene_df = group_by_gene(df)
    position_df = group_by_position(df)
    stats = compute_summary_statistics(df)
    sharing_matrix = snp_sharing_matrix(df) if df["Sample"].nunique() > 1 else None
    figures = build_all_figures(df, gene_df, sharing_matrix)

    export_paths = export_all(df, gene_df, position_df, output_dir)

    vcf_display_name = args.vcf[0] if len(args.vcf) == 1 else f"{len(args.vcf)} VCF files"
    report_paths = generate_reports(
        vcf_display_name, df, gene_df, position_df, stats, parse_report, figures,
        output_dir, make_pdf=not args.no_pdf,
    )

    print("\n=== SNP Annotation Analyser -- Run Complete ===")
    print(f"Total SNPs processed:   {stats['total_snps']}")
    print(f"Annotated:              {stats['annotated_snps']} ({stats['annotated_pct']}%)")
    print(f"Unannotated:            {stats['unannotated_snps']} ({stats['unannotated_pct']}%)")
    print(f"Affected genes:         {stats['affected_genes']}")
    print(f"Ti/Tv ratio:            {stats['ti_tv_ratio']}")
    print("\nOutput files:")
    for label, path in {**export_paths, **report_paths}.items():
        print(f"  {label:12s} -> {path}")
    print()

    return 0


def run_two_strain_mode(args: argparse.Namespace) -> int:
    """Two-strain independent analysis + coordinate-based comparison."""
    strain1_name = args.name1 or _strain_name_from_path(args.vcf1)
    strain2_name = args.name2 or _strain_name_from_path(args.vcf2)

    logger.info("Starting SNP Annotation Analyser CLI run (two-strain comparison mode)")
    logger.info("Strain 1: %s (%s)", strain1_name, args.vcf1)
    logger.info("Strain 2: %s (%s)", strain2_name, args.vcf2)

    try:
        df1, df2, report1, report2 = process_two_strains(args.vcf1, args.vcf2, strain1_name, strain2_name)
    except VCFValidationError as exc:
        logger.error("VCF validation failed: %s", exc)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        logger.exception("Unexpected error while parsing VCF(s)")
        print(f"ERROR: Unexpected failure while parsing VCF(s): {exc}", file=sys.stderr)
        return 1

    if df1.empty or df2.empty:
        empty_strain = strain1_name if df1.empty else strain2_name
        print(f"ERROR: No SNP records were extracted from {empty_strain}'s VCF.", file=sys.stderr)
        return 1

    if args.min_qual is not None or args.min_dp is not None:
        df1 = apply_quality_filters(df1, min_qual=args.min_qual, min_dp=args.min_dp)
        df2 = apply_quality_filters(df2, min_qual=args.min_qual, min_dp=args.min_dp)

    # Reference compatibility check (spec §26) -- warns, never blocks.
    compat = check_reference_compatibility(report1, report2, strain1_name, strain2_name)
    if compat.warnings:
        for w in compat.warnings:
            print(f"WARNING: {w}", file=sys.stderr)

    dirs = ensure_strain_output_dirs(args.output, strain1_name, strain2_name)

    #  Individual analysis, unchanged from single-VCF mode, run twice 
    gene_df1 = group_by_gene(df1)
    position_df1 = group_by_position(df1)
    stats1 = compute_summary_statistics(df1)
    figures1 = build_all_figures(df1, gene_df1)

    gene_df2 = group_by_gene(df2)
    position_df2 = group_by_position(df2)
    stats2 = compute_summary_statistics(df2)
    figures2 = build_all_figures(df2, gene_df2)

    export_strain_labeled(df1, gene_df1, position_df1, dirs["strain1"], strain1_name)
    export_strain_labeled(df2, gene_df2, position_df2, dirs["strain2"], strain2_name)

    generate_reports(Path(args.vcf1).name, df1, gene_df1, position_df1, stats1, report1,
                      figures1, dirs["strain1"], make_pdf=not args.no_pdf)
    generate_reports(Path(args.vcf2).name, df2, gene_df2, position_df2, stats2, report2,
                      figures2, dirs["strain2"], make_pdf=not args.no_pdf)

    # --- Comparison layer: reads df1/df2, never modifies them ---
    position_comparison = build_position_comparison(df1, df2, strain1_name, strain2_name)
    gene_comparison = build_gene_comparison(df1, df2, strain1_name, strain2_name)
    summary = comparison_summary(position_comparison, gene_comparison, strain1_name, strain2_name)

    export_comparison_all(position_comparison, gene_comparison, dirs["comparison"], strain1_name, strain2_name)

    comparison_figures = build_comparison_figures(df1, df2, gene_comparison, summary, strain1_name, strain2_name)
    generate_comparison_reports(
        strain1_name, strain2_name, Path(args.vcf1).name, Path(args.vcf2).name,
        stats1, stats2, position_comparison, gene_comparison, summary,
        comparison_figures, dirs["comparison"], make_pdf=not args.no_pdf,
        compatibility_warnings=compat.warnings,
    )

    print(f"\n=== SNP Annotation Analyser -- {strain1_name} vs {strain2_name} Comparison Complete ===")
    print(f"\n{strain1_name}: {stats1['total_snps']} SNPs, {stats1['affected_genes']} genes affected")
    print(f"{strain2_name}: {stats2['total_snps']} SNPs, {stats2['affected_genes']} genes affected")
    print(f"\nTotal unique positions:              {summary['total_unique_positions']}")
    print(f"Shared positions:                     {summary['shared_positions']}")
    print(f"  Identical mutations:                {summary['identical_mutations']}")
    print(f"  Same position, different mutation:  {summary['same_position_different_mutation']}")
    print(f"  Reference discrepancies:             {summary['reference_discrepancies']}")
    print(f"{strain1_name}-specific positions:{'':<10}{summary[f'{strain1_name}_specific']}")
    print(f"{strain2_name}-specific positions:{'':<10}{summary[f'{strain2_name}_specific']}")
    print(f"Shared genes:                          {summary['shared_genes']}")
    print(f"\nOutput written to: {ensure_output_dir(args.output)}")
    print(f"  {strain1_name}/   -> individual {strain1_name} analysis")
    print(f"  {strain2_name}/   -> individual {strain2_name} analysis")
    print(f"  comparison/  -> comparison tables and report")
    print()

    return 0


def run_multi_strain_mode(args: argparse.Namespace) -> int:
    """N-strain (2 to 10+) independent analysis + N-way coordinate-based
    comparison, generalising run_two_strain_mode above."""
    vcf_paths = args.vcfs
    strain_names = args.strain_names or [_strain_name_from_path(p) for p in vcf_paths]

    if len(strain_names) != len(vcf_paths):
        print(
            f"ERROR: Got {len(vcf_paths)} VCF file(s) but {len(strain_names)} "
            "--strain-names -- these must match 1:1.", file=sys.stderr,
        )
        return 1
    if len(set(strain_names)) != len(strain_names):
        print("ERROR: --strain-names must all be unique.", file=sys.stderr)
        return 1

    logger.info("Starting SNP Annotation Analyser CLI run (multi-strain, N=%d)", len(strain_names))
    for name, path in zip(strain_names, vcf_paths):
        logger.info("  %s -> %s", name, path)

    try:
        strains = process_multiple_strains(vcf_paths, strain_names)
    except VCFValidationError as exc:
        logger.error("VCF validation failed: %s", exc)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        logger.exception("Unexpected error while parsing VCF(s)")
        print(f"ERROR: Unexpected failure while parsing VCF(s): {exc}", file=sys.stderr)
        return 1

    empty_strains = [name for name, data in strains.items() if data["df"].empty]
    if empty_strains:
        print(f"ERROR: No SNP records extracted for: {', '.join(empty_strains)}", file=sys.stderr)
        return 1

    if args.min_qual is not None or args.min_dp is not None:
        for name in strains:
            strains[name]["df"] = apply_quality_filters(
                strains[name]["df"], min_qual=args.min_qual, min_dp=args.min_dp
            )

    # Pairwise reference compatibility checks -- warn, never block.
    all_warnings = []
    for i, name_a in enumerate(strain_names):
        for name_b in strain_names[i + 1:]:
            compat = check_reference_compatibility(
                strains[name_a]["parse_report"], strains[name_b]["parse_report"], name_a, name_b,
            )
            if compat.warnings:
                all_warnings.extend(compat.warnings)
                for w in compat.warnings:
                    print(f"WARNING ({name_a} vs {name_b}): {w}", file=sys.stderr)

    output_root = ensure_output_dir(args.output)

    # --- Individual analysis, unchanged pipeline, run once per strain ---
    per_strain_stats = {}
    vcf_filenames = {}
    for name, path in zip(strain_names, vcf_paths):
        df = strains[name]["df"]
        gene_df = group_by_gene(df)
        position_df = group_by_position(df)
        stats = compute_summary_statistics(df)
        figures = build_all_figures(df, gene_df)
        per_strain_stats[name] = stats
        vcf_filenames[name] = Path(path).name

        strain_dir = ensure_output_dir(output_root / name)
        export_strain_labeled(df, gene_df, position_df, strain_dir, name)
        generate_reports(Path(path).name, df, gene_df, position_df, stats,
                          strains[name]["parse_report"], figures, strain_dir, make_pdf=not args.no_pdf)

    # --- N-way comparison layer: reads all strain DataFrames, modifies none ---
    matrix = build_snp_matrix(strains)
    classified_matrix = classify_snp_matrix(matrix, strain_names)
    gene_comparison = build_gene_comparison_n_way(strains)
    summary = comparison_summary_n_way(classified_matrix, gene_comparison, strain_names)
    pairwise_matrix = pairwise_sharing_matrix_n_way(strains)

    multi_strain_dir = ensure_output_dir(output_root / "multi_strain")
    export_multi_strain_comparison(classified_matrix, gene_comparison, pairwise_matrix, multi_strain_dir, strain_names)

    figures = build_multi_strain_figures(strains, classified_matrix, gene_comparison, summary, pairwise_matrix)
    generate_multi_strain_reports(
        strain_names, vcf_filenames, per_strain_stats, classified_matrix, gene_comparison,
        summary, figures, multi_strain_dir, make_pdf=not args.no_pdf,
        compatibility_warnings=all_warnings,
    )

    print(f"\n=== SNP Annotation Analyser -- {len(strain_names)}-Strain Comparison Complete ===\n")
    for name in strain_names:
        s = per_strain_stats[name]
        print(f"{name}: {s['total_snps']} SNPs, {s['affected_genes']} genes affected")
    print(f"\nTotal unique positions:     {summary['total_unique_positions']}")
    print(f"Core positions (all {summary['num_strains']}):    {summary['core_positions_total']}")
    print(f"  Core, identical:          {summary['core_identical']}")
    print(f"  Core, variable:           {summary['core_variable']}")
    print(f"Partial sharing:            {summary['partial_positions']}")
    print(f"Reference discrepancies:    {summary['reference_discrepancies']}")
    print(f"Core genes:                 {summary['core_genes']} / {summary['total_genes']}")
    print("\nStrain-specific position counts:")
    for name, count in summary["strain_unique_counts"].items():
        print(f"  {name}: {count}")
    print(f"\nOutput written to: {output_root}")
    for name in strain_names:
        print(f"  {name}/  -> individual {name} analysis")
    print(f"  multi_strain/  -> N-way comparison tables and report")
    print()

    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)

    two_strain_mode = args.vcf1 is not None or args.vcf2 is not None
    multi_strain_mode = args.vcfs is not None
    single_vcf_mode = args.vcf is not None

    modes_selected = sum([two_strain_mode, multi_strain_mode, single_vcf_mode])
    if modes_selected > 1:
        print(
            "ERROR: Use exactly one of --vcf (single-VCF), --vcf1/--vcf2 "
            "(two-strain), or --vcfs (multi-strain).", file=sys.stderr,
        )
        return 1

    if multi_strain_mode:
        if len(args.vcfs) < 2:
            print("ERROR: --vcfs requires at least 2 VCF files. Use --vcf for a single VCF.", file=sys.stderr)
            return 1
        return run_multi_strain_mode(args)
    if two_strain_mode:
        if not (args.vcf1 and args.vcf2):
            print("ERROR: Two-strain comparison mode requires both --vcf1 and --vcf2.", file=sys.stderr)
            return 1
        return run_two_strain_mode(args)
    if single_vcf_mode:
        return run_single_vcf_mode(args)

    print(
        "ERROR: Provide --vcf (single-VCF), --vcf1/--vcf2 (two-strain), "
        "or --vcfs (multi-strain).", file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
