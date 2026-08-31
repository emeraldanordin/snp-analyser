"""
export.py
---------
CSV and Excel export helpers for the SNP table, gene-level summary, and
position-level summary. List/dict-valued columns (positions, allele
counts, etc.) are stringified for flat-file formats since CSV/Excel
can't natively hold nested structures.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from utils import ensure_output_dir, get_logger

logger = get_logger(__name__)


def _flatten_for_export(df: pd.DataFrame) -> pd.DataFrame:
    """Convert list/dict-valued cells to readable strings so CSV/Excel
    writers don't choke on nested Python objects."""
    flat = df.copy()
    for col in flat.columns:
        if flat[col].apply(lambda v: isinstance(v, (list, dict))).any():
            flat[col] = flat[col].apply(
                lambda v: "; ".join(map(str, v)) if isinstance(v, list)
                else "; ".join(f"{k}={val}" for k, val in v.items()) if isinstance(v, dict)
                else v
            )
    return flat


def export_csv(df: pd.DataFrame, path: str | Path) -> Path:
    """Export a DataFrame to CSV, flattening nested columns first."""
    path = Path(path)
    ensure_output_dir(path.parent)
    flat = _flatten_for_export(df)
    flat.to_csv(path, index=False)
    logger.info("Exported CSV: %s (%d rows)", path, len(flat))
    return path


def export_excel(df: pd.DataFrame, path: str | Path, sheet_name: str = "Data") -> Path:
    """Export a DataFrame to a single-sheet Excel file, flattening nested columns first."""
    path = Path(path)
    ensure_output_dir(path.parent)
    flat = _flatten_for_export(df)
    flat.to_excel(path, index=False, sheet_name=sheet_name, engine="openpyxl")
    logger.info("Exported Excel: %s (%d rows)", path, len(flat))
    return path


def export_multi_sheet_excel(sheets: dict[str, pd.DataFrame], path: str | Path) -> Path:
    """Export multiple DataFrames to one Excel workbook, one sheet each.
    Used for a combined 'full analysis' workbook (SNP table + gene summary
    + position summary in one file)."""
    path = Path(path)
    ensure_output_dir(path.parent)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            flat = _flatten_for_export(df)
            # Excel sheet names are capped at 31 characters
            safe_name = sheet_name[:31]
            flat.to_excel(writer, index=False, sheet_name=safe_name)
    logger.info("Exported multi-sheet Excel: %s (%d sheets)", path, len(sheets))
    return path


def export_all(
    snp_df: pd.DataFrame,
    gene_df: pd.DataFrame,
    position_df: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, Path]:
    """
    Convenience wrapper used by both the CLI and the Streamlit download
    buttons: writes snp_table.csv/.xlsx, gene_summary.csv, and
    position_summary.csv into output_dir. Returns a dict of the paths
    written, matching the CLI's expected output structure (§26).
    """
    out = ensure_output_dir(output_dir)
    paths = {}
    paths["snp_csv"] = export_csv(snp_df, out / "snp_table.csv")
    paths["snp_xlsx"] = export_excel(snp_df, out / "snp_table.xlsx", sheet_name="SNPs")
    paths["gene_csv"] = export_csv(gene_df, out / "gene_summary.csv")
    paths["position_csv"] = export_csv(position_df, out / "position_summary.csv")
    return paths


def export_strain_labeled(
    snp_df: pd.DataFrame,
    gene_df: pd.DataFrame,
    position_df: pd.DataFrame,
    output_dir: str | Path,
    strain_name: str,
) -> dict[str, Path]:
    """
    Same as export_all, but with strain-prefixed filenames (spec §29):
    <strain>_SNPs.csv/.xlsx, <strain>_gene_summary.csv,
    <strain>_position_summary.csv -- used for two-strain comparison runs
    where each strain's individual exports must stay clearly labeled and
    never overwrite the other strain's files.
    """
    out = ensure_output_dir(output_dir)
    paths = {}
    paths["snp_csv"] = export_csv(snp_df, out / f"{strain_name}_SNPs.csv")
    paths["snp_xlsx"] = export_excel(snp_df, out / f"{strain_name}_SNPs.xlsx", sheet_name="SNPs")
    paths["gene_csv"] = export_csv(gene_df, out / f"{strain_name}_gene_summary.csv")
    paths["position_csv"] = export_csv(position_df, out / f"{strain_name}_position_summary.csv")
    return paths


def export_comparison_all(
    position_comparison_df: pd.DataFrame,
    gene_comparison_df: pd.DataFrame,
    output_dir: str | Path,
    strain1_name: str,
    strain2_name: str,
) -> dict[str, Path]:
    """
    Export the comparison-layer outputs (spec §29): shared SNPs,
    strain-specific SNPs (split per strain), gene comparison, and the
    full position comparison table. Filenames use both strain names so
    they're self-documenting.
    """
    out = ensure_output_dir(output_dir)
    prefix = f"{strain1_name}_{strain2_name}"
    paths = {}

    shared = position_comparison_df[
        position_comparison_df["Comparison_Status"].isin(
            ["Identical mutation", "Same position, different mutation", "Reference allele discrepancy"]
        )
    ]
    strain1_specific = position_comparison_df[
        position_comparison_df["Comparison_Status"] == f"{strain1_name}-specific"
    ]
    strain2_specific = position_comparison_df[
        position_comparison_df["Comparison_Status"] == f"{strain2_name}-specific"
    ]

    paths["shared_csv"] = export_csv(shared, out / f"{prefix}_shared_SNPs.csv")
    paths["strain1_specific_csv"] = export_csv(strain1_specific, out / f"{prefix}_{strain1_name}_specific.csv")
    paths["strain2_specific_csv"] = export_csv(strain2_specific, out / f"{prefix}_{strain2_name}_specific.csv")
    paths["gene_comparison_csv"] = export_csv(gene_comparison_df, out / f"{prefix}_gene_comparison.csv")
    paths["position_comparison_csv"] = export_csv(position_comparison_df, out / f"{prefix}_position_comparison.csv")

    return paths


def export_multi_strain_comparison(
    classified_matrix: pd.DataFrame,
    gene_comparison_df: pd.DataFrame,
    pairwise_matrix: pd.DataFrame,
    output_dir: str | Path,
    strain_names: list[str],
) -> dict[str, Path]:
    """
    Export N-strain comparison outputs: the full SNP matrix (all
    strains' positions, comparison status), the core-genome-only subset
    (useful as a starting point for phylogenetic SNP alignment, per the
    original spec's forward-looking §32), per-strain-unique subsets, the
    gene comparison table, and the pairwise sharing matrix.
    """
    out = ensure_output_dir(output_dir)
    prefix = "_".join(strain_names) if len(strain_names) <= 4 else f"{len(strain_names)}strains"
    paths = {}

    paths["snp_matrix_csv"] = export_csv(classified_matrix, out / f"{prefix}_snp_matrix.csv")

    core_only = classified_matrix[
        classified_matrix["Comparison_Status"].isin(["Core, identical", "Core, variable"])
    ]
    paths["core_genome_csv"] = export_csv(core_only, out / f"{prefix}_core_genome_positions.csv")

    for name in strain_names:
        unique_subset = classified_matrix[classified_matrix["Comparison_Status"] == f"{name}-unique"]
        paths[f"{name}_unique_csv"] = export_csv(unique_subset, out / f"{prefix}_{name}_unique.csv")

    paths["gene_comparison_csv"] = export_csv(gene_comparison_df, out / f"{prefix}_gene_comparison.csv")
    paths["pairwise_sharing_csv"] = export_csv(
        pairwise_matrix.reset_index().rename(columns={"index": "Strain"}),
        out / f"{prefix}_pairwise_sharing_matrix.csv",
    )

    return paths
