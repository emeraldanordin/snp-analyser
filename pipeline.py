"""
pipeline.py
-----------
Single orchestration entry point that runs a VCF (or several) through
parse -> annotate -> classify and returns one clean pandas DataFrame,
plus the ParseReport(s). Both app.py (Streamlit) and main.py (CLI) call
this instead of re-implementing the wiring themselves, guaranteeing they
never diverge in behavior.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from annotation import annotate_rows
from classifier import classify_rows
from utils import get_logger
from vcf_parser import ParseReport, parse_vcf

logger = get_logger(__name__)


def process_vcf(path: str | Path) -> tuple[pd.DataFrame, ParseReport]:
    """Run one VCF file through the full parse -> annotate -> classify
    pipeline. Returns (snp_dataframe, parse_report)."""
    rows, report = parse_vcf(path)
    rows = annotate_rows(rows)
    rows = classify_rows(rows)
    df = pd.DataFrame(rows)
    return df, report


def process_vcf_with_strain_name(path: str | Path, strain_name: str) -> tuple[pd.DataFrame, ParseReport]:
    """
    Run one VCF through the full pipeline and tag every row with a
    user-assigned 'Strain' label, independent of whatever the VCF's own
    sample column happens to say.

    This exists because two different strains can share the exact same
    FORMAT/sample column name in their VCFs (e.g. both files call their
    sample 'isolate1') -- the VCF alone cannot distinguish them. The
    original 'Sample' column from the VCF header is preserved unchanged
    alongside the new 'Strain' column, so neither piece of information
    is lost.
    """
    df, report = process_vcf(path)
    if df.empty:
        df["Strain"] = pd.Series(dtype=str)
        return df, report
    df = df.copy()
    df["Strain"] = strain_name
    logger.info(
        "Tagged %d rows from %s as strain '%s' (VCF sample column: %s)",
        len(df), Path(path).name, strain_name, report.samples,
    )
    return df, report


def process_multiple_vcfs(paths: list[str | Path]) -> tuple[pd.DataFrame, list[ParseReport]]:
    """
    Run multiple VCF files (e.g. one per isolate/strain) through the pipeline and
    concatenate into one long-format DataFrame. Used for multi-isolate
    comparison when isolates arrive as separate files rather than one
    multi-sample VCF.
    """
    all_dfs = []
    reports = []
    for path in paths:
        df, report = process_vcf(path)
        all_dfs.append(df)
        reports.append(report)

    if not all_dfs:
        return pd.DataFrame(), reports

    combined = pd.concat(all_dfs, ignore_index=True)
    logger.info(
        "Combined %d VCF files into one table: %d total SNP rows across %d samples",
        len(paths), len(combined), combined["Sample"].nunique() if "Sample" in combined else 0,
    )
    return combined, reports


def apply_quality_filters(
    df: pd.DataFrame,
    min_qual: float | None = None,
    min_dp: int | None = None,
) -> pd.DataFrame:
    """
    Apply optional QC filters (§31) to the processed DataFrame. Never
    modifies the original VCF -- this only filters the in-memory analysis
    table, and the caller is responsible for reporting before/after counts
    (see statistics.compute_filter_impact).
    """
    filtered = df.copy()
    if min_qual is not None:
        filtered = filtered[filtered["QUAL"].fillna(0) >= min_qual]
    if min_dp is not None:
        filtered = filtered[filtered["INFO_DP"].fillna(0) >= min_dp]
    return filtered


def process_two_strains(
    path1: str | Path,
    path2: str | Path,
    name1: str,
    name2: str,
) -> tuple[pd.DataFrame, pd.DataFrame, ParseReport, ParseReport]:
    """
    Run two VCF files through the pipeline independently, each tagged
    with its own strain name. Returns two fully independent DataFrames
    (df1, df2) plus their ParseReports -- nothing here compares or links
    them; that is comparison.py's job, operating strictly read-only on
    these two results.
    """
    df1, report1 = process_vcf_with_strain_name(path1, name1)
    df2, report2 = process_vcf_with_strain_name(path2, name2)
    return df1, df2, report1, report2


def process_multiple_strains(
    paths: list[str | Path],
    names: list[str],
) -> dict[str, dict]:
    """
    Run any number (2, 10, 50...) of VCF files through the pipeline,
    each independently and each tagged with its own strain name.

    Returns a dict keyed by strain name:
        {
            "STB7A":  {"df": <DataFrame>, "parse_report": <ParseReport>},
            "STB36A": {"df": <DataFrame>, "parse_report": <ParseReport>},
            ...
        }

    This is the N-strain generalization of process_two_strains(). Each
    strain's dataset is fully independent of every other's -- nothing in
    this function compares or links them; that happens in
    comparison.build_snp_matrix() and friends, operating strictly
    read-only on the dict this function returns.

    Raises ValueError if paths and names have different lengths, or if
    any two strain names collide (comparison logic depends on names
    being unique to distinguish strains whose VCFs may share the same
    internal sample name).
    """
    if len(paths) != len(names):
        raise ValueError(
            f"Got {len(paths)} VCF path(s) but {len(names)} strain name(s) -- "
            "these lists must be the same length and in the same order."
        )
    if len(set(names)) != len(names):
        duplicates = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(
            f"Strain names must be unique; duplicate name(s) found: {', '.join(duplicates)}"
        )
    if len(paths) < 2:
        raise ValueError("process_multiple_strains requires at least 2 VCF files.")

    strains = {}
    for path, name in zip(paths, names):
        df, report = process_vcf_with_strain_name(path, name)
        strains[name] = {"df": df, "parse_report": report}
        logger.info("Processed strain '%s': %d SNPs", name, len(df))

    return strains
