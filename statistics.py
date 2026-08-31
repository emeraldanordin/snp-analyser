"""
statistics.py
-------------
Computes overall summary statistics for the SNP table: counts and
percentages for mutation categories, annotation coverage, transitions/
transversions, and affected gene/position counts.

Pure function: DataFrame in, dict out. Never mutates the input.
"""

from __future__ import annotations

import pandas as pd

from utils import get_logger

logger = get_logger(__name__)


def _pct(part: int, whole: int) -> float:
    return round(100 * part / whole, 1) if whole else 0.0


def compute_summary_statistics(df: pd.DataFrame) -> dict:
    """
    Compute the full summary statistics dict described in the spec:
    total SNPs, annotated/unannotated counts, affected genes, unique
    positions, per-category counts + percentages, and Ti/Tv counts/ratio.
    """
    if df.empty:
        return {
            "total_snps": 0,
            "annotated_snps": 0,
            "unannotated_snps": 0,
            "affected_genes": 0,
            "unique_positions": 0,
            "categories": {},
            "transitions": 0,
            "transversions": 0,
            "indel_complex": 0,
            "ti_tv_ratio": 0.0,
        }

    total = len(df)
    annotated = int((df["Effect_Raw"] != "Not annotated").sum())
    unannotated = total - annotated

    affected_genes = int(
        df.loc[df["Gene_Name"] != "Not annotated", "Gene_Name"].nunique()
    )
    unique_positions = int(df["POS"].nunique())

    category_counts = df["Mutation_Category"].value_counts().to_dict()
    categories = {
        cat: {"count": int(count), "percentage": _pct(count, total)}
        for cat, count in category_counts.items()
    }

    transitions = int((df["SNP_Type"] == "Transition").sum())
    transversions = int((df["SNP_Type"] == "Transversion").sum())
    indel_complex = int((df["SNP_Type"] == "Indel / Complex").sum())
    ti_tv_ratio = round(transitions / transversions, 3) if transversions > 0 else (
        float("inf") if transitions > 0 else 0.0
    )

    summary = {
        "total_snps": total,
        "annotated_snps": annotated,
        "unannotated_snps": unannotated,
        "annotated_pct": _pct(annotated, total),
        "unannotated_pct": _pct(unannotated, total),
        "affected_genes": affected_genes,
        "unique_positions": unique_positions,
        "categories": categories,
        "transitions": transitions,
        "transversions": transversions,
        "indel_complex": indel_complex,
        "ti_tv_ratio": ti_tv_ratio,
    }

    logger.info(
        "Summary statistics: %d total SNPs, %d genes affected, Ti/Tv=%s",
        total, affected_genes, ti_tv_ratio,
    )
    return summary


def compute_filter_impact(before_count: int, after_count: int) -> dict:
    """Small helper for the QC-filtering UI (§31): report before/after counts."""
    return {
        "records_before": before_count,
        "records_after": after_count,
        "records_removed": before_count - after_count,
        "pct_retained": _pct(after_count, before_count),
    }
