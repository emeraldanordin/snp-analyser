"""
comparison.py
-------------
Multi-isolate comparison logic. Works on the same long-format SNP
DataFrame used everywhere else (one row per SNP x Sample), so it doesn't
matter whether the isolate came from a single multi-sample VCF or were
concatenated from several single-sample VCFs by main.py/app.py.

Distinguishes explicitly between:
  1. Same gene, multiple SNPs           -> handled by grouping.group_by_gene
  2. Same genomic position, across isolate -> isolate_position_summary()
  3. Different ALT alleles at the same position -> grouping.group_by_position
     (allele breakdown) combined with isolate counts here.
"""

from __future__ import annotations

import pandas as pd

from utils import get_logger

logger = get_logger(__name__)


def get_isolates(df: pd.DataFrame) -> list[str]:
    if "Sample" not in df.columns:
        return []
    return sorted(df["Sample"].dropna().unique().tolist())


def snp_sharing_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each isolate pair, count SNPs (by snp_id) shared between them.
    Returns a square DataFrame (isolate x isolate) of shared-SNP counts,
    useful for the 'SNP sharing between isolates' visualisation.
    """
    isolates = get_isolates(df)
    if len(isolates) < 2:
        return pd.DataFrame()

    isolate_snp_sets = {
        iso: set(df.loc[df["Sample"] == iso, "snp_id"]) for iso in isolates
    }

    matrix = pd.DataFrame(index=isolates, columns=isolates, dtype=int)
    for iso_a in isolates:
        for iso_b in isolates:
            matrix.loc[iso_a, iso_b] = len(
                isolate_snp_sets[iso_a] & isolate_snp_sets[iso_b]
            )
    return matrix


def classify_snp_sharing(df: pd.DataFrame) -> pd.DataFrame:
    """
    For every distinct snp_id, determine how many isolates carry it and
    whether it's 'Shared' (>1 isolate) or 'Isolate-specific' (exactly 1).

    Returns one row per snp_id with:
        snp_id, CHROM, POS, REF, ALT, Gene_Name,
        num_isolates, isolates (list), sharing_status
    """
    if df.empty or "Sample" not in df.columns:
        return pd.DataFrame()

    total_isolates = df["Sample"].nunique()

    def _agg(g: pd.DataFrame) -> pd.Series:
        isolates = sorted(g["Sample"].dropna().unique().tolist())
        n = len(isolates)
        return pd.Series({
            "CHROM": g["CHROM"].iloc[0],
            "POS": g["POS"].iloc[0],
            "REF": g["REF"].iloc[0],
            "ALT": g["ALT"].iloc[0],
            "Gene_Name": g["Gene_Name"].iloc[0],
            "num_isolates": n,
            "isolates": isolates,
            "sharing_status": "Shared" if n > 1 else "Isolate-specific",
            "pct_isolates": round(100 * n / total_isolates, 1) if total_isolates else 0.0,
        })

    result = (
        df.groupby("snp_id", dropna=False)
        .apply(_agg, include_groups=False)
        .reset_index()
        .sort_values(["CHROM", "POS"])
        .reset_index(drop=True)
    )
    return result


def conserved_vs_variable_positions(df: pd.DataFrame) -> pd.DataFrame:
    """
    For every genomic position (CHROM, POS, REF), determine whether it's
    'Conserved' (every isolate that has a call there shows the same ALT)
    or 'Variable' (different isolates show different ALT alleles at the
    same position -- i.e. more than one ALT allele observed).
    """
    if df.empty:
        return pd.DataFrame()

    def _agg(g: pd.DataFrame) -> pd.Series:
        alt_alleles = sorted(g["ALT"].dropna().unique().tolist())
        return pd.Series({
            "Num_Alt_Alleles": len(alt_alleles),
            "Alt_Alleles": alt_alleles,
            "Status": "Variable" if len(alt_alleles) > 1 else "Conserved",
        })

    result = (
        df.groupby(["CHROM", "POS", "REF"], dropna=False)
        .apply(_agg, include_groups=False)
        .reset_index()
        .sort_values(["CHROM", "POS"])
        .reset_index(drop=True)
    )
    return result


def isolate_summary_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Per-isolate SNP counts, split by mutation category -- a quick
    per-sample overview table for the UI."""
    if df.empty or "Sample" not in df.columns:
        return pd.DataFrame()

    summary = (
        df.groupby("Sample")
        .agg(
            Total_SNPs=("snp_id", "count"),
            Missense=("Mutation_Category", lambda s: (s == "Missense").sum()),
            Synonymous=("Mutation_Category", lambda s: (s == "Silent / Synonymous").sum()),
            Nonsense=("Mutation_Category", lambda s: (s == "Nonsense / Stop gained").sum()),
            Frameshift=("Mutation_Category", lambda s: (s == "Frameshift").sum()),
            Transitions=("SNP_Type", lambda s: (s == "Transition").sum()),
            Transversions=("SNP_Type", lambda s: (s == "Transversion").sum()),
        )
        .reset_index()
    )
    return summary

# TWO-STRAIN COMPARISON 

REFERENCE_DISCREPANCY = "Reference allele discrepancy"
IDENTICAL_MUTATION = "Identical mutation"
SAME_POSITION_DIFFERENT_MUTATION = "Same position, different mutation"
STRAIN1_SPECIFIC = "{strain1}-specific"
STRAIN2_SPECIFIC = "{strain2}-specific"


def _position_level_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse one strain's SNP table to one row per (CHROM, POS): the set
    of REF alleles seen (normally exactly one), the set of ALT alleles
    seen (handles multi-allelic sites correctly per §25), and the full
    original row(s) for that position kept as a list of dicts so no
    original information is discarded.
    """
    if df.empty:
        return pd.DataFrame(columns=["CHROM", "POS", "REF_set", "ALT_set", "Gene_Name", "records"])

    def _agg(g: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "REF_set": sorted(g["REF"].dropna().unique().tolist()),
            "ALT_set": sorted(g["ALT"].dropna().unique().tolist()),
            "Gene_Name": g["Gene_Name"].iloc[0],
            "records": g.to_dict(orient="records"),
        })

    return (
        df.groupby(["CHROM", "POS"], dropna=False)
        .apply(_agg, include_groups=False)
        .reset_index()
    )


def _classify_position_status(row: pd.Series, strain1_name: str, strain2_name: str) -> str:
    """Apply the priority-ordered comparison_status rules to one merged
    position row (from the outer-merged position summaries)."""
    has1 = isinstance(row.get("REF_set_1"), list) and len(row["REF_set_1"]) > 0
    has2 = isinstance(row.get("REF_set_2"), list) and len(row["REF_set_2"]) > 0

    if has1 and not has2:
        return f"{strain1_name}-specific"
    if has2 and not has1:
        return f"{strain2_name}-specific"
    if not has1 and not has2:
        return "Unknown"  # should not occur from an outer merge

    ref1, ref2 = set(row["REF_set_1"]), set(row["REF_set_2"])
    if ref1 != ref2:
        return REFERENCE_DISCREPANCY

    alt1, alt2 = set(row["ALT_set_1"]), set(row["ALT_set_2"])
    if alt1 == alt2:
        return IDENTICAL_MUTATION
    return SAME_POSITION_DIFFERENT_MUTATION


def build_position_comparison(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    strain1_name: str,
    strain2_name: str,
) -> pd.DataFrame:
    """
    Core comparison engine (spec §5-§10, §13, §39): build one row per
    unique (CHROM, POS) across BOTH strains, linked by genomic coordinate.

    Each row retains the COMPLETE original record(s) from both strains
    (not just REF/ALT) via the '{strain}_records' list-of-dict columns,
    so the full annotation, quality, and genotype information for both
    strains stays reachable from the comparison table -- nothing is
    reduced to position/REF/ALT only.

    Multi-allelic sites are handled by comparing ALT allele SETS rather
    than single values, so A>G vs A>T at the same position is correctly
    distinguished rather than merged or falsely matched.
    """
    if df1.empty and df2.empty:
        return pd.DataFrame()

    pos1 = _position_level_summary(df1)
    pos2 = _position_level_summary(df2)

    merged = pd.merge(
        pos1, pos2, on=["CHROM", "POS"], how="outer",
        suffixes=("_1", "_2"), indicator=False,
    )

    # Gene name: prefer whichever strain has it; if both have it and they
    # differ (shouldn't normally happen for the same coordinate, but VCFs
    # can disagree), keep both, separated, rather than silently picking one.
    def _resolve_gene(row):
        g1 = row.get("Gene_Name_1")
        g2 = row.get("Gene_Name_2")
        g1_valid = isinstance(g1, str) and g1 not in (None, "Not annotated") and pd.notna(g1)
        g2_valid = isinstance(g2, str) and g2 not in (None, "Not annotated") and pd.notna(g2)
        if g1_valid and g2_valid:
            return g1 if g1 == g2 else f"{g1} / {g2}"
        if g1_valid:
            return g1
        if g2_valid:
            return g2
        return "Not annotated"

    merged["Gene_Name"] = merged.apply(_resolve_gene, axis=1)
    merged["Comparison_Status"] = merged.apply(
        lambda r: _classify_position_status(r, strain1_name, strain2_name), axis=1
    )

    # Rename the per-strain columns to use the actual strain names, so the
    # output table is self-documenting rather than exposing internal "_1"/"_2".
    merged = merged.rename(columns={
        "REF_set_1": f"{strain1_name}_REF_alleles",
        "ALT_set_1": f"{strain1_name}_ALT_alleles",
        "records_1": f"{strain1_name}_records",
        "REF_set_2": f"{strain2_name}_REF_alleles",
        "ALT_set_2": f"{strain2_name}_ALT_alleles",
        "records_2": f"{strain2_name}_records",
    })

    # Fill missing (strain-specific) side with empty lists rather than NaN,
    # so downstream code can check `len(...) == 0` uniformly.
    for col in [
        f"{strain1_name}_REF_alleles", f"{strain1_name}_ALT_alleles", f"{strain1_name}_records",
        f"{strain2_name}_REF_alleles", f"{strain2_name}_ALT_alleles", f"{strain2_name}_records",
    ]:
        merged[col] = merged[col].apply(lambda v: v if isinstance(v, list) else [])

    merged = merged.drop(columns=["Gene_Name_1", "Gene_Name_2"], errors="ignore")
    merged = merged.sort_values(["CHROM", "POS"]).reset_index(drop=True)

    logger.info(
        "Position comparison built: %d unique positions (%s: %d, %s: %d, shared coords: %d)",
        len(merged), strain1_name, len(pos1), strain2_name, len(pos2),
        int(((merged[f"{strain1_name}_REF_alleles"].apply(len) > 0) &
             (merged[f"{strain2_name}_REF_alleles"].apply(len) > 0)).sum()),
    )
    return merged


def build_gene_comparison(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    strain1_name: str,
    strain2_name: str,
) -> pd.DataFrame:
    """
    Gene-level comparison (spec §11-§12): for every gene appearing in
    EITHER strain, report each strain's SNP count, shared vs.
    strain-specific POSITIONS within that gene (distinct from shared
    SNPs -- a shared gene does not imply a shared position), and
    per-strain missense/synonymous/transition/transversion counts.
    """
    if df1.empty and df2.empty:
        return pd.DataFrame()

    def _gene_positions(df: pd.DataFrame) -> dict[str, set]:
        if df.empty:
            return {}
        return {
            gene: set(g["POS"].unique())
            for gene, g in df[df["Gene_Name"] != "Not annotated"].groupby("Gene_Name")
        }

    def _gene_counts(df: pd.DataFrame) -> dict[str, dict]:
        if df.empty:
            return {}
        result = {}
        for gene, g in df[df["Gene_Name"] != "Not annotated"].groupby("Gene_Name"):
            result[gene] = {
                "total": len(g),
                "missense": int((g["Mutation_Category"] == "Missense").sum()),
                "synonymous": int((g["Mutation_Category"] == "Silent / Synonymous").sum()),
                "transitions": int((g["SNP_Type"] == "Transition").sum()),
                "transversions": int((g["SNP_Type"] == "Transversion").sum()),
            }
        return result

    pos1_by_gene = _gene_positions(df1)
    pos2_by_gene = _gene_positions(df2)
    counts1 = _gene_counts(df1)
    counts2 = _gene_counts(df2)

    all_genes = sorted(set(pos1_by_gene) | set(pos2_by_gene))
    rows = []
    for gene in all_genes:
        p1 = pos1_by_gene.get(gene, set())
        p2 = pos2_by_gene.get(gene, set())
        c1 = counts1.get(gene, {"total": 0, "missense": 0, "synonymous": 0, "transitions": 0, "transversions": 0})
        c2 = counts2.get(gene, {"total": 0, "missense": 0, "synonymous": 0, "transitions": 0, "transversions": 0})

        rows.append({
            "Gene_Name": gene,
            f"{strain1_name}_SNPs": c1["total"],
            f"{strain2_name}_SNPs": c2["total"],
            "Shared_Positions": len(p1 & p2),
            f"{strain1_name}_Specific_Positions": len(p1 - p2),
            f"{strain2_name}_Specific_Positions": len(p2 - p1),
            f"{strain1_name}_Missense": c1["missense"],
            f"{strain2_name}_Missense": c2["missense"],
            f"{strain1_name}_Synonymous": c1["synonymous"],
            f"{strain2_name}_Synonymous": c2["synonymous"],
            f"{strain1_name}_Transitions": c1["transitions"],
            f"{strain2_name}_Transitions": c2["transitions"],
            f"{strain1_name}_Transversions": c1["transversions"],
            f"{strain2_name}_Transversions": c2["transversions"],
            "Shared_Gene": bool(p1) and bool(p2),
        })

    gene_comparison = pd.DataFrame(rows).sort_values(
        [f"{strain1_name}_SNPs", f"{strain2_name}_SNPs"], ascending=False
    ).reset_index(drop=True)

    logger.info(
        "Gene comparison built: %d genes total, %d shared between %s and %s",
        len(gene_comparison), int(gene_comparison["Shared_Gene"].sum()), strain1_name, strain2_name,
    )
    return gene_comparison


def comparison_summary(
    position_comparison_df: pd.DataFrame,
    gene_comparison_df: pd.DataFrame,
    strain1_name: str,
    strain2_name: str,
) -> dict:
    """
    Build the overall comparison summary dict (spec §15, §38): total
    unique positions, shared/specific counts, identical mutations,
    same-position-different-allele counts, reference discrepancies,
    and shared gene count.
    """
    if position_comparison_df.empty:
        return {
            "total_unique_positions": 0,
            "shared_positions": 0,
            f"{strain1_name}_specific": 0,
            f"{strain2_name}_specific": 0,
            "identical_mutations": 0,
            "same_position_different_mutation": 0,
            "reference_discrepancies": 0,
            "shared_genes": 0,
        }

    status_counts = position_comparison_df["Comparison_Status"].value_counts().to_dict()

    return {
        "total_unique_positions": len(position_comparison_df),
        "shared_positions": status_counts.get(IDENTICAL_MUTATION, 0)
            + status_counts.get(SAME_POSITION_DIFFERENT_MUTATION, 0)
            + status_counts.get(REFERENCE_DISCREPANCY, 0),
        f"{strain1_name}_specific": status_counts.get(f"{strain1_name}-specific", 0),
        f"{strain2_name}_specific": status_counts.get(f"{strain2_name}-specific", 0),
        "identical_mutations": status_counts.get(IDENTICAL_MUTATION, 0),
        "same_position_different_mutation": status_counts.get(SAME_POSITION_DIFFERENT_MUTATION, 0),
        "reference_discrepancies": status_counts.get(REFERENCE_DISCREPANCY, 0),
        "shared_genes": int(gene_comparison_df["Shared_Gene"].sum()) if not gene_comparison_df.empty else 0,
    }


def format_dna_change(ref_alleles: list, alt_alleles: list) -> str:
    """
    Format a strain's REF/ALT allele lists at one position into a
    readable nucleotide-change string, e.g. ['G'], ['A'] -> 'G>A'.

    Handles the cases that come up in practice:
    - No variant at this position (empty lists) -> '-'
    - The normal single-allele case -> 'G>A'
    - Multi-allelic sites (rare, but the ALT set can have >1 entry) ->
      'G>A, G>C' -- every observed change is listed rather than picking
      one arbitrarily.
    """
    if not ref_alleles or not alt_alleles:
        return "-"
    ref = ref_alleles[0] if len(ref_alleles) == 1 else "/".join(sorted(ref_alleles))
    changes = [f"{ref}>{alt}" for alt in sorted(alt_alleles)]
    return ", ".join(changes)


def add_dna_change_columns(matrix: pd.DataFrame, strain_names: list[str]) -> pd.DataFrame:
    """
    Add a '{strain}_DNA_Change' display column to a SNP matrix (from
    build_snp_matrix or build_position_comparison) for every strain,
    computed from that strain's existing REF/ALT allele-list columns.
    Does not modify or remove the underlying REF/ALT list columns --
    this is purely an additional, more readable column for display,
    exports, and reports.
    """
    matrix = matrix.copy()
    for name in strain_names:
        ref_col = f"{name}_REF_alleles"
        alt_col = f"{name}_ALT_alleles"
        if ref_col in matrix.columns and alt_col in matrix.columns:
            matrix[f"{name}_DNA_Change"] = matrix.apply(
                lambda row: format_dna_change(row[ref_col], row[alt_col]), axis=1
            )
    return matrix


def get_position_detail(
    position_comparison_df: pd.DataFrame,
    chrom: str,
    pos: int,
) -> dict | None:
    """
    Retrieve the full comparison detail for one genomic position, for
    the Position Detail view (spec §19-§20). Returns None if the
    position isn't in the comparison table at all.
    """
    match = position_comparison_df[
        (position_comparison_df["CHROM"] == chrom) & (position_comparison_df["POS"] == pos)
    ]
    if match.empty:
        return None
    return match.iloc[0].to_dict()


def filter_comparison(
    position_comparison_df: pd.DataFrame,
    status_filter: list[str] | None = None,
    gene_filter: list[str] | None = None,
    search_text: str | None = None,
) -> pd.DataFrame:
    """Apply the comparison-page filters (spec §34) to the position
    comparison table. Never mutates the input."""
    filtered = position_comparison_df.copy()
    if status_filter:
        filtered = filtered[filtered["Comparison_Status"].isin(status_filter)]
    if gene_filter:
        filtered = filtered[filtered["Gene_Name"].isin(gene_filter)]
    if search_text:
        mask = filtered.astype(str).apply(
            lambda col: col.str.contains(search_text, case=False, na=False)
        ).any(axis=1)
        filtered = filtered[mask]
    return filtered



# MULTI-STRAIN COMPARISON 

CORE_IDENTICAL = "Core, identical"
CORE_VARIABLE = "Core, variable"
PARTIAL = "Partial"


def build_snp_matrix(strains: dict[str, dict]) -> pd.DataFrame:
    """
    Build the N-way genomic position matrix: one row per unique
    (CHROM, POS) across ALL strains, with one pair of columns per strain
    ('{strain}_REF_alleles', '{strain}_ALT_alleles', '{strain}_records')
    holding that strain's complete data at that position (empty lists if
    the strain has no variant there).

    This is the N-strain generalisation of build_position_comparison().
    No strain's original DataFrame is read more than once or modified.
    """
    strain_names = list(strains.keys())
    if len(strain_names) < 2:
        raise ValueError("build_snp_matrix requires at least 2 strains.")

    per_strain_summaries = {
        name: _position_level_summary(data["df"]) for name, data in strains.items()
    }

    # Outer-merge all per-strain position summaries together on (CHROM, POS).
    matrix = None
    for name, summary in per_strain_summaries.items():
        renamed = summary.rename(columns={
            "REF_set": f"{name}_REF_alleles",
            "ALT_set": f"{name}_ALT_alleles",
            "records": f"{name}_records",
            "Gene_Name": f"{name}_Gene_Name",
        })
        if matrix is None:
            matrix = renamed
        else:
            matrix = pd.merge(matrix, renamed, on=["CHROM", "POS"], how="outer")

    # Fill missing (strain-absent) cells with empty lists, not NaN.
    for name in strain_names:
        for suffix in ("_REF_alleles", "_ALT_alleles", "_records"):
            col = f"{name}{suffix}"
            matrix[col] = matrix[col].apply(lambda v: v if isinstance(v, list) else [])

    # Resolve a single Gene_Name per position across whichever strains have it.
    gene_cols = [f"{name}_Gene_Name" for name in strain_names]

    def _resolve_gene(row):
        genes = {
            row[c] for c in gene_cols
            if isinstance(row[c], str) and row[c] not in (None, "Not annotated") and pd.notna(row[c])
        }
        if not genes:
            return "Not annotated"
        return genes.pop() if len(genes) == 1 else " / ".join(sorted(genes))

    matrix["Gene_Name"] = matrix.apply(_resolve_gene, axis=1)
    matrix = matrix.drop(columns=gene_cols)

    matrix = matrix.sort_values(["CHROM", "POS"]).reset_index(drop=True)
    logger.info("SNP matrix built: %d unique positions across %d strains", len(matrix), len(strain_names))
    return matrix


def classify_snp_matrix(matrix: pd.DataFrame, strain_names: list[str]) -> pd.DataFrame:
    """
    Add Comparison_Status, Num_Strains_Present, and Strains_Present
    columns to a SNP matrix built by build_snp_matrix(). Kept as a
    separate step from build_snp_matrix() so the matrix itself and its
    classification can be tested/reasoned about independently.
    """
    matrix = matrix.copy()
    n_total = len(strain_names)

    def _classify(row) -> pd.Series:
        present = [
            name for name in strain_names
            if len(row[f"{name}_REF_alleles"]) > 0
        ]
        n_present = len(present)

        if n_present == 1:
            status = f"{present[0]}-unique"
        elif n_present == n_total:
            refs = set()
            alts = set()
            for name in present:
                refs.update(row[f"{name}_REF_alleles"])
                alts_by_strain = frozenset(row[f"{name}_ALT_alleles"])
                alts.add(alts_by_strain)
            if len(refs) > 1:
                status = REFERENCE_DISCREPANCY
            elif len(alts) == 1:
                status = CORE_IDENTICAL
            else:
                status = CORE_VARIABLE
        else:
            status = f"{PARTIAL} ({n_present}/{n_total})"

        return pd.Series({
            "Num_Strains_Present": n_present,
            "Strains_Present": present,
            "Comparison_Status": status,
        })

    classification = matrix.apply(_classify, axis=1)
    matrix = pd.concat([matrix, classification], axis=1)
    return matrix


def build_gene_comparison_n_way(strains: dict[str, dict]) -> pd.DataFrame:
    """
    N-way gene-level comparison: for every gene appearing in ANY strain,
    report each strain's SNP count in that gene (one column per strain),
    how many strains are affected in that gene, and how many genomic
    positions within that gene are shared by ALL affected strains vs.
    unique to a subset.
    """
    strain_names = list(strains.keys())

    def _gene_positions(df: pd.DataFrame) -> dict[str, set]:
        if df.empty:
            return {}
        return {
            gene: set(g["POS"].unique())
            for gene, g in df[df["Gene_Name"] != "Not annotated"].groupby("Gene_Name")
        }

    def _gene_counts(df: pd.DataFrame) -> dict[str, int]:
        if df.empty:
            return {}
        return df[df["Gene_Name"] != "Not annotated"].groupby("Gene_Name").size().to_dict()

    positions_by_strain = {name: _gene_positions(data["df"]) for name, data in strains.items()}
    counts_by_strain = {name: _gene_counts(data["df"]) for name, data in strains.items()}

    all_genes = sorted(set().union(*[set(p.keys()) for p in positions_by_strain.values()]))

    rows = []
    for gene in all_genes:
        row = {"Gene_Name": gene}
        strains_with_gene = []
        gene_position_sets = []

        for name in strain_names:
            count = counts_by_strain[name].get(gene, 0)
            row[f"{name}_SNPs"] = count
            if count > 0:
                strains_with_gene.append(name)
                gene_position_sets.append(positions_by_strain[name][gene])

        row["Num_Strains_Affected"] = len(strains_with_gene)
        row["Strains_Affected"] = strains_with_gene

        if gene_position_sets:
            shared_by_all_affected = set.intersection(*gene_position_sets)
            union_positions = set.union(*gene_position_sets)
            row["Shared_Positions_All_Affected"] = len(shared_by_all_affected)
            row["Total_Unique_Positions"] = len(union_positions)
        else:
            row["Shared_Positions_All_Affected"] = 0
            row["Total_Unique_Positions"] = 0

        row["Core_Gene"] = len(strains_with_gene) == len(strain_names)
        rows.append(row)

    gene_comparison = pd.DataFrame(rows)
    snp_cols = [f"{name}_SNPs" for name in strain_names]
    gene_comparison["Total_SNPs_All_Strains"] = gene_comparison[snp_cols].sum(axis=1)
    gene_comparison = gene_comparison.sort_values(
        "Total_SNPs_All_Strains", ascending=False
    ).reset_index(drop=True)

    logger.info(
        "N-way gene comparison built: %d genes total, %d core (affected in all %d strains)",
        len(gene_comparison), int(gene_comparison["Core_Gene"].sum()), len(strain_names),
    )
    return gene_comparison


def comparison_summary_n_way(
    classified_matrix: pd.DataFrame,
    gene_comparison_df: pd.DataFrame,
    strain_names: list[str],
) -> dict:
    """
    Overall N-way comparison summary: total unique positions, core
    genome size (identical + variable), per-strain-unique counts,
    partial-sharing counts, reference discrepancies, core genes, and
    per-strain total SNP counts.
    """
    if classified_matrix.empty:
        return {
            "num_strains": len(strain_names),
            "total_unique_positions": 0,
            "core_identical": 0,
            "core_variable": 0,
            "core_positions_total": 0,
            "reference_discrepancies": 0,
            "partial_positions": 0,
            "strain_unique_counts": {name: 0 for name in strain_names},
            "core_genes": 0,
            "total_genes": 0,
        }

    status_counts = classified_matrix["Comparison_Status"].value_counts().to_dict()

    strain_unique_counts = {
        name: status_counts.get(f"{name}-unique", 0) for name in strain_names
    }
    partial_positions = sum(
        count for status, count in status_counts.items() if status.startswith(PARTIAL)
    )

    return {
        "num_strains": len(strain_names),
        "total_unique_positions": len(classified_matrix),
        "core_identical": status_counts.get(CORE_IDENTICAL, 0),
        "core_variable": status_counts.get(CORE_VARIABLE, 0),
        "core_positions_total": status_counts.get(CORE_IDENTICAL, 0) + status_counts.get(CORE_VARIABLE, 0),
        "reference_discrepancies": status_counts.get(REFERENCE_DISCREPANCY, 0),
        "partial_positions": partial_positions,
        "strain_unique_counts": strain_unique_counts,
        "core_genes": int(gene_comparison_df["Core_Gene"].sum()) if not gene_comparison_df.empty else 0,
        "total_genes": len(gene_comparison_df),
    }


def pairwise_sharing_matrix_n_way(strains: dict[str, dict]) -> pd.DataFrame:
    """
    NxN matrix of shared genomic-position counts between every pair of
    strains (by CHROM+POS, not just identical mutations) -- generalises
    snp_sharing_matrix() (which works on Sample within one DataFrame) to
    N independently-parsed strain DataFrames. Used for the sharing
    heatmap visualisation, which scales far better than a Venn diagram
    once N is more than 3-4.
    """
    strain_names = list(strains.keys())
    position_sets = {
        name: set(zip(data["df"]["CHROM"], data["df"]["POS"]))
        for name, data in strains.items()
    }

    matrix = pd.DataFrame(index=strain_names, columns=strain_names, dtype=int)
    for a in strain_names:
        for b in strain_names:
            matrix.loc[a, b] = len(position_sets[a] & position_sets[b])
    return matrix


def get_position_detail_n_way(
    classified_matrix: pd.DataFrame,
    chrom: str,
    pos: int,
    strain_names: list[str],
) -> dict | None:
    """Retrieve the full comparison detail for one genomic position
    across all N strains -- the N-strain Position Detail view."""
    match = classified_matrix[
        (classified_matrix["CHROM"] == chrom) & (classified_matrix["POS"] == pos)
    ]
    if match.empty:
        return None
    return match.iloc[0].to_dict()


def filter_snp_matrix(
    classified_matrix: pd.DataFrame,
    status_filter: list[str] | None = None,
    gene_filter: list[str] | None = None,
    strain_filter: list[str] | None = None,
    search_text: str | None = None,
) -> pd.DataFrame:
    """Filter the N-way SNP matrix by comparison status, gene, strains
    present, or free text. Never mutates the input."""
    filtered = classified_matrix.copy()
    if status_filter:
        filtered = filtered[filtered["Comparison_Status"].isin(status_filter)]
    if gene_filter:
        filtered = filtered[filtered["Gene_Name"].isin(gene_filter)]
    if strain_filter:
        filtered = filtered[
            filtered["Strains_Present"].apply(lambda present: any(s in present for s in strain_filter))
        ]
    if search_text:
        display_cols = [c for c in filtered.columns if not c.endswith("_records")]
        mask = filtered[display_cols].astype(str).apply(
            lambda col: col.str.contains(search_text, case=False, na=False)
        ).any(axis=1)
        filtered = filtered[mask]
    return filtered
