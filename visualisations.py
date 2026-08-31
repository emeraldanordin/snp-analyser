"""
visualisations.py
------------------
Builds every Plotly figure used in the app and in the HTML/PDF report,
so the chart-construction logic exists in exactly one place. Each
function takes the SNP DataFrame (and sometimes the gene summary) and
returns a plotly.graph_objects.Figure.
"""

from __future__ import annotations

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

SCIENTIFIC_TEMPLATE = "simple_white"


def fig_effect_distribution(df: pd.DataFrame) -> go.Figure:
    """Mutation effect (raw SnpEff term) distribution -- bar chart."""
    counts = df["Effect_Raw"].value_counts().reset_index()
    counts.columns = ["Effect", "Count"]
    fig = px.bar(
        counts, x="Effect", y="Count", template=SCIENTIFIC_TEMPLATE,
        title="Mutation Effect Distribution (raw SnpEff terms)",
    )
    fig.update_layout(xaxis_tickangle=-45, xaxis_title="Effect", yaxis_title="Number of SNPs")
    return fig


def fig_category_distribution(df: pd.DataFrame) -> go.Figure:
    """Simplified mutation category distribution -- pie chart."""
    counts = df["Mutation_Category"].value_counts().reset_index()
    counts.columns = ["Category", "Count"]
    fig = px.pie(
        counts, names="Category", values="Count", template=SCIENTIFIC_TEMPLATE,
        title="Mutation Category Distribution",
    )
    return fig


def fig_transition_transversion(df: pd.DataFrame) -> go.Figure:
    """Transition vs. transversion counts -- bar chart."""
    counts = df["SNP_Type"].value_counts().reset_index()
    counts.columns = ["SNP_Type", "Count"]
    fig = px.bar(
        counts, x="SNP_Type", y="Count", template=SCIENTIFIC_TEMPLATE,
        title="Transition vs. Transversion", color="SNP_Type",
    )
    fig.update_layout(xaxis_title="", yaxis_title="Number of SNPs", showlegend=False)
    return fig


def fig_snps_per_gene(gene_df: pd.DataFrame, top_n: int | None = None) -> go.Figure:
    """SNPs per gene -- bar chart, optionally limited to top N genes."""
    data = gene_df.sort_values("Total_SNPs", ascending=False)
    if top_n:
        data = data.head(top_n)
    fig = px.bar(
        data, x="Gene_Name", y="Total_SNPs", template=SCIENTIFIC_TEMPLATE,
        title=f"SNPs per Gene{f' (Top {top_n})' if top_n else ''}",
    )
    fig.update_layout(xaxis_tickangle=-45, xaxis_title="Gene", yaxis_title="Number of SNPs")
    return fig


def fig_top_genes(gene_df: pd.DataFrame, top_n: int = 15) -> go.Figure:
    """Top N genes by SNP count -- horizontal bar chart for readability."""
    data = gene_df.sort_values("Total_SNPs", ascending=False).head(top_n)
    fig = px.bar(
        data.sort_values("Total_SNPs"), x="Total_SNPs", y="Gene_Name", orientation="h",
        template=SCIENTIFIC_TEMPLATE, title=f"Top {top_n} Genes by Number of SNPs",
    )
    fig.update_layout(xaxis_title="Number of SNPs", yaxis_title="Gene")
    return fig


def fig_snp_distribution_along_chromosome(df: pd.DataFrame, bin_size: int = 50_000) -> go.Figure:
    """SNP density along the chromosome, binned by position -- histogram."""
    fig = px.histogram(
        df, x="POS", nbins=max(int(df["POS"].max() / bin_size), 10) if len(df) else 10,
        template=SCIENTIFIC_TEMPLATE,
        title="SNP Distribution Along Chromosome",
    )
    fig.update_layout(xaxis_title="Genomic Position", yaxis_title="Number of SNPs")
    return fig


def fig_missense_vs_synonymous(df: pd.DataFrame) -> go.Figure:
    """Missense vs synonymous counts -- simple comparison bar chart."""
    subset = df[df["Mutation_Category"].isin(["Missense", "Silent / Synonymous"])]
    counts = subset["Mutation_Category"].value_counts().reset_index()
    counts.columns = ["Category", "Count"]
    fig = px.bar(
        counts, x="Category", y="Count", template=SCIENTIFIC_TEMPLATE,
        title="Missense vs. Synonymous SNPs", color="Category",
    )
    fig.update_layout(xaxis_title="", yaxis_title="Number of SNPs", showlegend=False)
    return fig


def fig_isolate_sharing(sharing_matrix: pd.DataFrame) -> go.Figure | None:
    """SNP sharing between isolates -- heatmap. Returns None if fewer than
    2 isolates (nothing meaningful to show)."""
    if sharing_matrix is None or sharing_matrix.empty:
        return None
    fig = px.imshow(
        sharing_matrix.astype(int), text_auto=True, template=SCIENTIFIC_TEMPLATE,
        title="SNP Sharing Between Isolates", color_continuous_scale="Blues",
        labels=dict(color="Shared SNPs"),
    )
    return fig


def build_all_figures(df: pd.DataFrame, gene_df: pd.DataFrame, sharing_matrix: pd.DataFrame | None = None) -> dict[str, go.Figure]:
    """Build every figure at once, keyed by name -- used by report.py to
    embed the full figure set into the HTML/PDF report."""
    figs = {
        "effect_distribution": fig_effect_distribution(df),
        "category_distribution": fig_category_distribution(df),
        "transition_transversion": fig_transition_transversion(df),
        "snps_per_gene": fig_snps_per_gene(gene_df, top_n=25),
        "top_genes": fig_top_genes(gene_df, top_n=15),
        "chromosome_distribution": fig_snp_distribution_along_chromosome(df),
        "missense_vs_synonymous": fig_missense_vs_synonymous(df),
    }
    if sharing_matrix is not None and not sharing_matrix.empty:
        sharing_fig = fig_isolate_sharing(sharing_matrix)
        if sharing_fig is not None:
            figs["isolate_sharing"] = sharing_fig
    return figs


# TWO-STRAIN COMPARISON VISUALISATIONS

def fig_snp_counts_by_strain(df1: pd.DataFrame, df2: pd.DataFrame, name1: str, name2: str) -> go.Figure:
    """SNP counts by strain -- simple side-by-side bar chart."""
    data = pd.DataFrame({"Strain": [name1, name2], "Total_SNPs": [len(df1), len(df2)]})
    fig = px.bar(
        data, x="Strain", y="Total_SNPs", template=SCIENTIFIC_TEMPLATE,
        title="Total SNP Count by Strain", color="Strain",
    )
    fig.update_layout(showlegend=False, yaxis_title="Number of SNPs")
    return fig


def fig_shared_vs_specific(summary: dict, name1: str, name2: str) -> go.Figure:
    """Shared vs. strain-specific SNP positions -- bar chart. This is the
    Plotly-based equivalent of a Venn diagram (spec §18): rather than an
    unreliable 2-circle overlap render, a labeled bar chart communicates
    the same three categories unambiguously."""
    data = pd.DataFrame({
        "Category": [f"{name1}-specific", "Shared", f"{name2}-specific"],
        "Positions": [
            summary.get(f"{name1}_specific", 0),
            summary.get("identical_mutations", 0) + summary.get("same_position_different_mutation", 0)
            + summary.get("reference_discrepancies", 0),
            summary.get(f"{name2}_specific", 0),
        ],
    })
    fig = px.bar(
        data, x="Category", y="Positions", template=SCIENTIFIC_TEMPLATE,
        title="Shared vs. Strain-Specific Genomic Positions", color="Category",
    )
    fig.update_layout(showlegend=False, yaxis_title="Number of Positions", xaxis_title="")
    return fig


def fig_mutation_category_by_strain(df1: pd.DataFrame, df2: pd.DataFrame, name1: str, name2: str) -> go.Figure:
    """Mutation category distribution, grouped side-by-side per strain."""
    c1 = df1["Mutation_Category"].value_counts().reset_index()
    c1.columns = ["Category", "Count"]
    c1["Strain"] = name1
    c2 = df2["Mutation_Category"].value_counts().reset_index()
    c2.columns = ["Category", "Count"]
    c2["Strain"] = name2
    combined = pd.concat([c1, c2], ignore_index=True)

    fig = px.bar(
        combined, x="Category", y="Count", color="Strain", barmode="group",
        template=SCIENTIFIC_TEMPLATE, title="Mutation Category Distribution by Strain",
    )
    fig.update_layout(xaxis_tickangle=-45, xaxis_title="", yaxis_title="Number of SNPs")
    return fig


def fig_missense_synonymous_comparison(df1: pd.DataFrame, df2: pd.DataFrame, name1: str, name2: str) -> go.Figure:
    """Missense vs synonymous counts, compared side-by-side per strain."""
    rows = []
    for name, df in [(name1, df1), (name2, df2)]:
        rows.append({"Strain": name, "Category": "Missense", "Count": int((df["Mutation_Category"] == "Missense").sum())})
        rows.append({"Strain": name, "Category": "Silent / Synonymous", "Count": int((df["Mutation_Category"] == "Silent / Synonymous").sum())})
    data = pd.DataFrame(rows)
    fig = px.bar(
        data, x="Category", y="Count", color="Strain", barmode="group",
        template=SCIENTIFIC_TEMPLATE, title="Missense vs. Synonymous SNPs by Strain",
    )
    fig.update_layout(xaxis_title="", yaxis_title="Number of SNPs")
    return fig


def fig_titv_comparison(df1: pd.DataFrame, df2: pd.DataFrame, name1: str, name2: str) -> go.Figure:
    """Transition/transversion counts, compared side-by-side per strain."""
    rows = []
    for name, df in [(name1, df1), (name2, df2)]:
        rows.append({"Strain": name, "SNP_Type": "Transition", "Count": int((df["SNP_Type"] == "Transition").sum())})
        rows.append({"Strain": name, "SNP_Type": "Transversion", "Count": int((df["SNP_Type"] == "Transversion").sum())})
    data = pd.DataFrame(rows)
    fig = px.bar(
        data, x="SNP_Type", y="Count", color="Strain", barmode="group",
        template=SCIENTIFIC_TEMPLATE, title="Transition / Transversion Counts by Strain",
    )
    fig.update_layout(xaxis_title="", yaxis_title="Number of SNPs")
    return fig


def fig_shared_genes_venn_style(gene_comparison_df: pd.DataFrame, name1: str, name2: str) -> go.Figure:
    """Shared vs strain-specific GENES (distinct from shared positions
    above) -- same Venn-style bar chart approach applied at the gene level."""
    if gene_comparison_df.empty:
        return go.Figure()

    col1 = f"{name1}_SNPs"
    col2 = f"{name2}_SNPs"
    shared = int(gene_comparison_df["Shared_Gene"].sum())
    only1 = int(((gene_comparison_df[col1] > 0) & (gene_comparison_df[col2] == 0)).sum())
    only2 = int(((gene_comparison_df[col2] > 0) & (gene_comparison_df[col1] == 0)).sum())

    data = pd.DataFrame({
        "Category": [f"{name1}-only genes", "Shared genes", f"{name2}-only genes"],
        "Genes": [only1, shared, only2],
    })
    fig = px.bar(
        data, x="Category", y="Genes", template=SCIENTIFIC_TEMPLATE,
        title="Shared vs. Strain-Specific Affected Genes", color="Category",
    )
    fig.update_layout(showlegend=False, yaxis_title="Number of Genes", xaxis_title="")
    return fig


def fig_top_genes_comparison(gene_comparison_df: pd.DataFrame, name1: str, name2: str, top_n: int = 15) -> go.Figure:
    """Top genes by combined SNP count across both strains, shown as
    grouped bars per strain."""
    if gene_comparison_df.empty:
        return go.Figure()

    col1 = f"{name1}_SNPs"
    col2 = f"{name2}_SNPs"
    data = gene_comparison_df.copy()
    data["Combined"] = data[col1] + data[col2]
    top = data.sort_values("Combined", ascending=False).head(top_n)

    melted = top.melt(
        id_vars=["Gene_Name"], value_vars=[col1, col2],
        var_name="Strain", value_name="SNPs",
    )
    melted["Strain"] = melted["Strain"].str.replace("_SNPs", "", regex=False)

    fig = px.bar(
        melted, x="Gene_Name", y="SNPs", color="Strain", barmode="group",
        template=SCIENTIFIC_TEMPLATE, title=f"Top {top_n} Genes by Combined SNP Count",
    )
    fig.update_layout(xaxis_tickangle=-45, xaxis_title="Gene", yaxis_title="Number of SNPs")
    return fig


def fig_chromosome_distribution_by_strain(df1: pd.DataFrame, df2: pd.DataFrame, name1: str, name2: str, bin_size: int = 50_000) -> go.Figure:
    """SNP density along the chromosome, overlaid per strain."""
    d1 = df1[["POS"]].copy()
    d1["Strain"] = name1
    d2 = df2[["POS"]].copy()
    d2["Strain"] = name2
    combined = pd.concat([d1, d2], ignore_index=True)

    max_pos = combined["POS"].max() if len(combined) else 1
    fig = px.histogram(
        combined, x="POS", color="Strain", barmode="overlay", opacity=0.6,
        nbins=max(int(max_pos / bin_size), 10), template=SCIENTIFIC_TEMPLATE,
        title="SNP Distribution Along Chromosome by Strain",
    )
    fig.update_layout(xaxis_title="Genomic Position", yaxis_title="Number of SNPs")
    return fig


def build_comparison_figures(
    df1: pd.DataFrame, df2: pd.DataFrame,
    gene_comparison_df: pd.DataFrame, summary: dict,
    name1: str, name2: str,
) -> dict[str, go.Figure]:
    """Build the full set of comparison-page figures at once."""
    return {
        "snp_counts_by_strain": fig_snp_counts_by_strain(df1, df2, name1, name2),
        "shared_vs_specific": fig_shared_vs_specific(summary, name1, name2),
        "category_by_strain": fig_mutation_category_by_strain(df1, df2, name1, name2),
        "missense_synonymous_comparison": fig_missense_synonymous_comparison(df1, df2, name1, name2),
        "titv_comparison": fig_titv_comparison(df1, df2, name1, name2),
        "shared_genes": fig_shared_genes_venn_style(gene_comparison_df, name1, name2),
        "top_genes_comparison": fig_top_genes_comparison(gene_comparison_df, name1, name2),
        "chromosome_distribution_by_strain": fig_chromosome_distribution_by_strain(df1, df2, name1, name2),
    }

# MULTI-STRAIN COMPARISON VISUALISATIONS


def fig_snps_per_strain_n(strains: dict[str, dict]) -> go.Figure:
    """Total SNP count per strain -- bar chart, any number of strains."""
    strain_names = list(strains.keys())
    counts = [len(strains[name]["df"]) for name in strain_names]
    data = pd.DataFrame({"Strain": strain_names, "Total_SNPs": counts})
    fig = px.bar(
        data, x="Strain", y="Total_SNPs", template=SCIENTIFIC_TEMPLATE,
        title="Total SNP Count by Strain",
    )
    fig.update_layout(xaxis_tickangle=-45, yaxis_title="Number of SNPs", showlegend=False)
    return fig


def fig_pairwise_sharing_heatmap(pairwise_matrix: pd.DataFrame) -> go.Figure:
    """
    NxN heatmap of shared genomic positions between every pair of
    strains. This is the N-strain replacement for the two-strain
    Venn-style bar chart -- it scales cleanly to 10+ strains where a
    Venn diagram or grouped bar chart would become unreadable.
    """
    if pairwise_matrix is None or pairwise_matrix.empty:
        return go.Figure()
    fig = px.imshow(
        pairwise_matrix.astype(int), text_auto=True, template=SCIENTIFIC_TEMPLATE,
        title="Pairwise Shared Genomic Positions Between Strains",
        color_continuous_scale="Blues", labels=dict(color="Shared Positions"),
    )
    fig.update_layout(xaxis_title="", yaxis_title="")
    return fig


def fig_core_genome_summary(summary: dict) -> go.Figure:
    """
    Core genome vs. accessory/partial vs. strain-specific -- the N-strain
    replacement for the two-strain 'shared vs specific' bar. Breaks the
    total position count into: core (identical), core (variable),
    partial (shared by some but not all), and total strain-specific.
    """
    total_unique = sum(summary["strain_unique_counts"].values())
    data = pd.DataFrame({
        "Category": ["Core, identical", "Core, variable", "Partial sharing", "Strain-specific (combined)"],
        "Positions": [
            summary["core_identical"],
            summary["core_variable"],
            summary["partial_positions"],
            total_unique,
        ],
    })
    fig = px.bar(
        data, x="Category", y="Positions", template=SCIENTIFIC_TEMPLATE,
        title=f"Genomic Position Categories Across {summary['num_strains']} Strains",
        color="Category",
    )
    fig.update_layout(showlegend=False, xaxis_title="", yaxis_title="Number of Positions")
    return fig


def fig_strain_unique_counts_n(summary: dict) -> go.Figure:
    """Strain-specific position counts, one bar per strain -- shows which
    strains carry the most/fewest unique variants."""
    strain_names = list(summary["strain_unique_counts"].keys())
    counts = list(summary["strain_unique_counts"].values())
    data = pd.DataFrame({"Strain": strain_names, "Unique_Positions": counts})
    fig = px.bar(
        data, x="Strain", y="Unique_Positions", template=SCIENTIFIC_TEMPLATE,
        title="Strain-Specific (Unique) Genomic Positions",
    )
    fig.update_layout(xaxis_tickangle=-45, yaxis_title="Number of Unique Positions", showlegend=False)
    return fig


def fig_genes_by_strain_count(gene_comparison_df: pd.DataFrame, strain_names: list[str]) -> go.Figure:
    """
    Histogram: how many genes are affected in exactly k of the N
    strains. A gene affected in all N strains is 'core'; a gene affected
    in only 1 is strain-specific. This is the gene-level analog of the
    core-genome-size concept, useful at a glance for 10+ strains.
    """
    if gene_comparison_df.empty:
        return go.Figure()
    counts = gene_comparison_df["Num_Strains_Affected"].value_counts().sort_index().reset_index()
    counts.columns = ["Num_Strains_Affected", "Num_Genes"]
    fig = px.bar(
        counts, x="Num_Strains_Affected", y="Num_Genes", template=SCIENTIFIC_TEMPLATE,
        title="Genes by Number of Strains Affected",
    )
    fig.update_layout(
        xaxis_title=f"Number of strains affected (out of {len(strain_names)})",
        yaxis_title="Number of genes", showlegend=False,
    )
    return fig


def fig_presence_absence_heatmap(classified_matrix: pd.DataFrame, strain_names: list[str], max_positions: int = 300) -> go.Figure:
    """
    Presence/absence heatmap: genomic positions (rows, limited to
    max_positions for renderability) x strains (columns), 1 = variant
    present, 0 = absent. Sorted by Num_Strains_Present so core positions
    cluster together visually. For very large position counts (10+
    strains can easily produce 5,000+ unique positions), this is capped
    and the caller should note the cap in surrounding text.
    """
    if classified_matrix.empty:
        return go.Figure()

    display = classified_matrix.sort_values("Num_Strains_Present", ascending=False).head(max_positions)
    z = []
    labels = []
    for _, row in display.iterrows():
        z.append([1 if s in row["Strains_Present"] else 0 for s in strain_names])
        labels.append(f"{row['CHROM']}:{row['POS']}")

    fig = go.Figure(data=go.Heatmap(
        z=z, x=strain_names, y=labels, colorscale=[[0, "#eaf2f8"], [1, "#1f5f8b"]],
        showscale=False,
    ))
    fig.update_layout(
        template=SCIENTIFIC_TEMPLATE,
        title=f"Presence/Absence Matrix (top {len(display)} positions by strain count)",
        xaxis_title="Strain", yaxis_title="Genomic Position",
        yaxis=dict(showticklabels=len(display) <= 60),  # avoid unreadable label spam
        height=max(400, min(20 * len(display), 1400)),
    )
    return fig


def build_multi_strain_figures(
    strains: dict[str, dict],
    classified_matrix: pd.DataFrame,
    gene_comparison_df: pd.DataFrame,
    summary: dict,
    pairwise_matrix: pd.DataFrame,
) -> dict[str, go.Figure]:
    """Build the full set of N-strain comparison figures at once."""
    strain_names = list(strains.keys())
    return {
        "snps_per_strain": fig_snps_per_strain_n(strains),
        "pairwise_sharing": fig_pairwise_sharing_heatmap(pairwise_matrix),
        "core_genome_summary": fig_core_genome_summary(summary),
        "strain_unique_counts": fig_strain_unique_counts_n(summary),
        "genes_by_strain_count": fig_genes_by_strain_count(gene_comparison_df, strain_names),
        "presence_absence": fig_presence_absence_heatmap(classified_matrix, strain_names),
    }
