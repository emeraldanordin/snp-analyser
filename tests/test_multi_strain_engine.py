"""
tests/test_multi_strain_engine.py
------------------------------------
Tests for the N-strain (multi-strain) comparison engine in comparison.py.
Covers synthetic cases at N=4 (core/partial/unique classification) plus
consistency checks against the real bundled STB7A/STB36A/STB_DEMO3 VCFs.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from comparison import (
    CORE_IDENTICAL,
    CORE_VARIABLE,
    build_gene_comparison_n_way,
    build_snp_matrix,
    classify_snp_matrix,
    comparison_summary_n_way,
    filter_snp_matrix,
    get_position_detail_n_way,
    pairwise_sharing_matrix_n_way,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
STB7A_VCF = PROJECT_ROOT / "data" / "STB7A.vcf"
STB36A_VCF = PROJECT_ROOT / "data" / "STB36A.vcf"
DEMO3_VCF = PROJECT_ROOT / "data" / "STB_DEMO3.vcf"


def _row(pos, ref, alt, gene="Not annotated", effect="Not annotated",
         category="Unknown / Not annotated", snp_type="Transition"):
    return {
        "CHROM": "Chromosome", "POS": pos, "REF": ref, "ALT": alt,
        "Gene_Name": gene, "Effect_Raw": effect, "Mutation_Category": category,
        "SNP_Type": snp_type, "HGVS_p": "Not annotated",
        "snp_id": f"Chromosome:{pos}:{ref}:{alt}", "Sample": "isolate1",
        "QUAL": 200.0, "INFO_DP": 100,
    }


def _make_strains(strain_positions: dict[str, list[tuple]]) -> dict[str, dict]:
    """Build a synthetic `strains` dict (as pipeline.process_multiple_strains
    would return) from {strain_name: [(pos, ref, alt, gene), ...]}."""
    strains = {}
    for name, records in strain_positions.items():
        rows = [_row(pos, ref, alt, gene=gene) for (pos, ref, alt, gene) in records]
        strains[name] = {"df": pd.DataFrame(rows), "parse_report": None}
    return strains


# ---------------------------------------------------------------------------
# 4-strain synthetic test: core, partial, and unique positions
# ---------------------------------------------------------------------------
def test_four_strain_core_partial_unique_classification():
    strains = _make_strains({
        "A": [(100, "G", "A", "geneX"), (200, "C", "T", "geneY"), (500, "A", "G", "geneZ")],
        "B": [(100, "G", "A", "geneX"), (200, "C", "T", "geneY")],
        "C": [(100, "G", "A", "geneX"), (300, "T", "C", "geneW")],
        "D": [(100, "G", "A", "geneX")],
    })

    matrix = build_snp_matrix(strains)
    classified = classify_snp_matrix(matrix, ["A", "B", "C", "D"])

    # Position 100: all 4 strains, same REF/ALT -> Core, identical
    row_100 = classified[classified["POS"] == 100].iloc[0]
    assert row_100["Comparison_Status"] == CORE_IDENTICAL
    assert row_100["Num_Strains_Present"] == 4

    # Position 200: A and B only -> Partial (2/4)
    row_200 = classified[classified["POS"] == 200].iloc[0]
    assert row_200["Comparison_Status"] == "Partial (2/4)"
    assert set(row_200["Strains_Present"]) == {"A", "B"}

    # Position 500: A only -> A-unique
    row_500 = classified[classified["POS"] == 500].iloc[0]
    assert row_500["Comparison_Status"] == "A-unique"

    # Position 300: C only -> C-unique
    row_300 = classified[classified["POS"] == 300].iloc[0]
    assert row_300["Comparison_Status"] == "C-unique"


def test_core_variable_when_alleles_differ_across_all_strains():
    """All strains have a variant at the position, but not the same ALT --
    must be 'Core, variable', not 'Core, identical'."""
    strains = _make_strains({
        "A": [(100, "G", "A", "geneX")],
        "B": [(100, "G", "C", "geneX")],
        "C": [(100, "G", "A", "geneX")],
    })
    matrix = build_snp_matrix(strains)
    classified = classify_snp_matrix(matrix, ["A", "B", "C"])
    row_100 = classified[classified["POS"] == 100].iloc[0]
    assert row_100["Comparison_Status"] == CORE_VARIABLE
    assert row_100["Comparison_Status"] != CORE_IDENTICAL


def test_reference_discrepancy_at_n_way():
    """Strains disagree on REF at the same position -- must be flagged,
    not silently classified as core-variable."""
    strains = _make_strains({
        "A": [(100, "G", "A", "geneX")],
        "B": [(100, "A", "T", "geneX")],  # different REF
        "C": [(100, "G", "C", "geneX")],
    })
    matrix = build_snp_matrix(strains)
    classified = classify_snp_matrix(matrix, ["A", "B", "C"])
    row_100 = classified[classified["POS"] == 100].iloc[0]
    assert row_100["Comparison_Status"] == "Reference allele discrepancy"


# ---------------------------------------------------------------------------
# Gene-level N-way comparison
# ---------------------------------------------------------------------------
def test_gene_comparison_n_way_core_gene_detection():
    strains = _make_strains({
        "A": [(100, "G", "A", "gyrA"), (200, "C", "T", "gyrA")],
        "B": [(150, "T", "C", "gyrA")],  # different position, same gene
        "C": [(999, "A", "G", "otherGene")],  # no gyrA at all
    })
    gene_cmp = build_gene_comparison_n_way(strains)

    gyra_row = gene_cmp[gene_cmp["Gene_Name"] == "gyrA"].iloc[0]
    assert gyra_row["A_SNPs"] == 2
    assert gyra_row["B_SNPs"] == 1
    assert gyra_row["C_SNPs"] == 0
    assert gyra_row["Num_Strains_Affected"] == 2
    assert gyra_row["Core_Gene"] == False  # not affected in strain C  # noqa: E712

    other_row = gene_cmp[gene_cmp["Gene_Name"] == "otherGene"].iloc[0]
    assert other_row["Num_Strains_Affected"] == 1


# ---------------------------------------------------------------------------
# Summary consistency
# ---------------------------------------------------------------------------
def test_summary_totals_are_internally_consistent():
    strains = _make_strains({
        "A": [(100, "G", "A", "geneX"), (200, "C", "T", "geneY"), (500, "A", "G", "geneZ")],
        "B": [(100, "G", "A", "geneX"), (200, "C", "T", "geneY")],
        "C": [(100, "G", "A", "geneX"), (300, "T", "C", "geneW")],
        "D": [(100, "G", "A", "geneX")],
    })
    matrix = build_snp_matrix(strains)
    classified = classify_snp_matrix(matrix, ["A", "B", "C", "D"])
    gene_cmp = build_gene_comparison_n_way(strains)
    summary = comparison_summary_n_way(classified, gene_cmp, ["A", "B", "C", "D"])

    total_from_categories = (
        summary["core_identical"] + summary["core_variable"]
        + summary["partial_positions"] + sum(summary["strain_unique_counts"].values())
        + summary["reference_discrepancies"]
    )
    assert total_from_categories == summary["total_unique_positions"]


def test_position_detail_n_way_retrieval():
    strains = _make_strains({
        "A": [(100, "G", "A", "geneX")],
        "B": [(100, "G", "A", "geneX")],
    })
    matrix = build_snp_matrix(strains)
    classified = classify_snp_matrix(matrix, ["A", "B"])
    detail = get_position_detail_n_way(classified, "Chromosome", 100, ["A", "B"])
    assert detail is not None
    assert len(detail["A_records"]) == 1
    assert len(detail["B_records"]) == 1
    assert detail["A_records"][0]["REF"] == "G"


def test_filter_snp_matrix_by_status_and_strain():
    strains = _make_strains({
        "A": [(100, "G", "A", "geneX"), (500, "A", "G", "geneZ")],
        "B": [(100, "G", "A", "geneX")],
    })
    matrix = build_snp_matrix(strains)
    classified = classify_snp_matrix(matrix, ["A", "B"])

    core_only = filter_snp_matrix(classified, status_filter=[CORE_IDENTICAL])
    assert len(core_only) == 1
    assert core_only.iloc[0]["POS"] == 100

    a_specific = filter_snp_matrix(classified, strain_filter=["A"])
    assert set(a_specific["POS"]) == {100, 500}


# ---------------------------------------------------------------------------
# N=2 must reduce to identical results as the dedicated two-strain engine
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not (STB7A_VCF.exists() and STB36A_VCF.exists()),
    reason="Bundled strain VCFs not present",
)
def test_n_way_engine_matches_two_strain_engine_at_n_equals_2():
    from pipeline import process_two_strains, process_multiple_strains
    from comparison import build_position_comparison, build_gene_comparison, comparison_summary

    df1, df2, r1, r2 = process_two_strains(STB7A_VCF, STB36A_VCF, "STB7A", "STB36A")
    two_way_pos = build_position_comparison(df1, df2, "STB7A", "STB36A")
    two_way_gene = build_gene_comparison(df1, df2, "STB7A", "STB36A")
    two_way_summary = comparison_summary(two_way_pos, two_way_gene, "STB7A", "STB36A")

    strains = process_multiple_strains([STB7A_VCF, STB36A_VCF], ["STB7A", "STB36A"])
    n_way_matrix = build_snp_matrix(strains)
    n_way_classified = classify_snp_matrix(n_way_matrix, ["STB7A", "STB36A"])
    n_way_gene = build_gene_comparison_n_way(strains)
    n_way_summary = comparison_summary_n_way(n_way_classified, n_way_gene, ["STB7A", "STB36A"])

    assert n_way_summary["total_unique_positions"] == two_way_summary["total_unique_positions"]
    assert n_way_summary["core_positions_total"] == two_way_summary["shared_positions"]
    assert n_way_summary["strain_unique_counts"]["STB7A"] == two_way_summary["STB7A_specific"]
    assert n_way_summary["strain_unique_counts"]["STB36A"] == two_way_summary["STB36A_specific"]


# ---------------------------------------------------------------------------
# Real 3-strain data (STB7A, STB36A, STB_DEMO3 -- a synthetic third strain
# derived from real data, used to genuinely exercise N=3 logic)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not (STB7A_VCF.exists() and STB36A_VCF.exists() and DEMO3_VCF.exists()),
    reason="Bundled 3-strain VCFs not present",
)
def test_real_three_strain_comparison_is_consistent():
    from pipeline import process_multiple_strains

    strains = process_multiple_strains(
        [STB7A_VCF, STB36A_VCF, DEMO3_VCF], ["STB7A", "STB36A", "STB_DEMO3"]
    )
    matrix = build_snp_matrix(strains)
    classified = classify_snp_matrix(matrix, list(strains.keys()))
    gene_cmp = build_gene_comparison_n_way(strains)
    summary = comparison_summary_n_way(classified, gene_cmp, list(strains.keys()))

    # Genuine 3-way partial sharing must appear (not collapsed to core/unique only)
    assert summary["partial_positions"] > 0

    # Every unique position accounted for exactly once
    total = (
        summary["core_identical"] + summary["core_variable"]
        + summary["partial_positions"] + sum(summary["strain_unique_counts"].values())
        + summary["reference_discrepancies"]
    )
    assert total == summary["total_unique_positions"]

    # No original strain DataFrame mutated by building the comparison
    for name, data in strains.items():
        assert "Comparison_Status" not in data["df"].columns


def test_pairwise_sharing_matrix_diagonal_equals_self_count():
    strains = _make_strains({
        "A": [(100, "G", "A", "geneX"), (200, "C", "T", "geneY")],
        "B": [(100, "G", "A", "geneX")],
    })
    pairwise = pairwise_sharing_matrix_n_way(strains)
    assert pairwise.loc["A", "A"] == 2
    assert pairwise.loc["B", "B"] == 1
    assert pairwise.loc["A", "B"] == 1
    assert pairwise.loc["B", "A"] == 1  # symmetric
