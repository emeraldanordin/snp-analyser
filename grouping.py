"""
grouping.py
-----------
Two independent grouping systems over the SNP-level DataFrame:

1. Gene-level grouping: aggregates every SNP that falls in the same gene,
   WITHOUT merging individual SNPs -- each stays a distinct row in the
   detail lists.

2. Position-level grouping: aggregates by (CHROM, POS, REF), which is a
   different key from gene grouping. This is what lets us recognise that
   9304 G>A and 9304 G>C are the same genomic position with two different
   alternative alleles, without conflating that with gene-level grouping.

Both functions are pure: input SNP DataFrame -> output summary DataFrame.
Neither mutates the input.
"""

from __future__ import annotations

import pandas as pd

from utils import get_logger

logger = get_logger(__name__)

NOT_ANNOTATED = "Not annotated"


def group_by_gene(df: pd.DataFrame) -> pd.DataFrame:
    """
    Group the SNP table by Gene_Name and compute per-gene summary stats.
    SNPs with Gene_Name == 'Not annotated' are grouped together under
    that label so they remain visible/countable, but are clearly not a
    real gene.

    Returns one row per gene with:
        Gene_Name, Gene_ID, Total_SNPs, Unique_Positions,
        Missense, Synonymous, Nonsense, Frameshift, Other_Effects,
        Transitions, Transversions, Ti_Tv_Ratio,
        Positions (list), DNA_Changes (list), Protein_Changes (list),
        Effects (list)
    """
    if df.empty:
        return pd.DataFrame()

    def _agg_gene(g: pd.DataFrame) -> pd.Series:
        total = len(g)
        transitions = int((g["SNP_Type"] == "Transition").sum())
        transversions = int((g["SNP_Type"] == "Transversion").sum())
        ti_tv = round(transitions / transversions, 3) if transversions > 0 else float("inf") if transitions > 0 else 0.0

        return pd.Series({
            "Gene_ID": g["Gene_ID"].iloc[0] if g["Gene_ID"].nunique() == 1 else "; ".join(sorted(g["Gene_ID"].dropna().unique())),
            "Total_SNPs": total,
            "Unique_Positions": int(g["POS"].nunique()),
            "Missense": int((g["Mutation_Category"] == "Missense").sum()),
            "Synonymous": int((g["Mutation_Category"] == "Silent / Synonymous").sum()),
            "Nonsense": int((g["Mutation_Category"] == "Nonsense / Stop gained").sum()),
            "Frameshift": int((g["Mutation_Category"] == "Frameshift").sum()),
            "Other_Effects": int((~g["Mutation_Category"].isin(
                ["Missense", "Silent / Synonymous", "Nonsense / Stop gained", "Frameshift"]
            )).sum()),
            "Transitions": transitions,
            "Transversions": transversions,
            "Ti_Tv_Ratio": ti_tv,
            "Positions": sorted(g["POS"].unique().tolist()),
            "DNA_Changes": g["DNA_Change"].tolist(),
            "Protein_Changes": g["HGVS_p"].tolist(),
            "Effects": g["Effect_Raw"].tolist(),
        })

    grouped = (
        df.groupby("Gene_Name", dropna=False)
        .apply(_agg_gene, include_groups=False)
        .reset_index()
        .rename(columns={"Gene_Name": "Gene_Name"})
    )
    grouped = grouped.sort_values("Total_SNPs", ascending=False).reset_index(drop=True)

    logger.info("Gene-level grouping: %d genes/groups from %d SNPs", len(grouped), len(df))
    return grouped


def group_by_position(df: pd.DataFrame) -> pd.DataFrame:
    """
    Group the SNP table by (CHROM, POS, REF) -- genomic position, NOT gene.
    This correctly treats e.g. 9304 G>A and 9304 G>C as the same position
    with two alternative alleles, and computes per-allele isolate counts
    and percentages when a 'Sample' column with multiple distinct isolates
    is present.

    Returns one row per (CHROM, POS, REF) with:
        CHROM, POS, REF, Num_Alt_Alleles, Alt_Alleles (list),
        Allele_Isolate_Counts (dict allele -> count),
        Allele_Isolate_Percentages (dict allele -> %),
        Gene_Name, Effect_Raw, Mutation_Category
        (gene/effect are taken from the first row per allele; if a
        position has >1 distinct gene/effect across alleles, all are
        listed rather than silently picking one)
    """
    if df.empty:
        return pd.DataFrame()

    total_isolates = df["Sample"].nunique() if "Sample" in df.columns else 1

    def _agg_position(g: pd.DataFrame) -> pd.Series:
        alt_alleles = sorted(g["ALT"].dropna().unique().tolist())

        allele_counts = {}
        allele_pcts = {}
        for allele in alt_alleles:
            isolates_with_allele = g.loc[g["ALT"] == allele, "Sample"].nunique() if "Sample" in g.columns else 1
            allele_counts[allele] = int(isolates_with_allele)
            allele_pcts[allele] = round(100 * isolates_with_allele / total_isolates, 1) if total_isolates else 0.0

        genes = sorted(g["Gene_Name"].dropna().unique().tolist())
        effects = sorted(g["Effect_Raw"].dropna().unique().tolist())
        categories = sorted(g["Mutation_Category"].dropna().unique().tolist())

        return pd.Series({
            "Num_Alt_Alleles": len(alt_alleles),
            "Alt_Alleles": alt_alleles,
            "Allele_Isolate_Counts": allele_counts,
            "Allele_Isolate_Percentages": allele_pcts,
            "Gene_Name": genes[0] if len(genes) == 1 else "; ".join(genes) if genes else NOT_ANNOTATED,
            "Effect_Raw": effects[0] if len(effects) == 1 else "; ".join(effects) if effects else NOT_ANNOTATED,
            "Mutation_Category": categories[0] if len(categories) == 1 else "; ".join(categories) if categories else "Unknown / Not annotated",
        })

    grouped = (
        df.groupby(["CHROM", "POS", "REF"], dropna=False)
        .apply(_agg_position, include_groups=False)
        .reset_index()
    )
    grouped = grouped.sort_values(["CHROM", "POS"]).reset_index(drop=True)

    logger.info(
        "Position-level grouping: %d unique positions from %d SNPs (%d total isolates)",
        len(grouped), len(df), total_isolates,
    )
    return grouped
