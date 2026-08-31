"""
tests/test_comparison_engine.py
---------------------------------
Tests for the two-strain comparison engine in comparison.py. Implements
the four synthetic test cases specified for this feature, plus checks
against the real bundled STB7A/STB36A VCFs.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from comparison import (
    build_gene_comparison,
    build_position_comparison,
    comparison_summary,
    get_position_detail,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
STB7A_VCF = PROJECT_ROOT / "data" / "STB7A.vcf"
STB36A_VCF = PROJECT_ROOT / "data" / "STB36A.vcf"


def _row(pos, ref, alt, gene="Not annotated", effect="Not annotated",
         category="Unknown / Not annotated", snp_type="Transition", hgvs_p="Not annotated"):
    return {
        "CHROM": "Chromosome", "POS": pos, "REF": ref, "ALT": alt,
        "Gene_Name": gene, "Effect_Raw": effect, "Mutation_Category": category,
        "SNP_Type": snp_type, "HGVS_p": hgvs_p, "snp_id": f"Chromosome:{pos}:{ref}:{alt}",
        "Sample": "isolate1", "QUAL": 200.0, "INFO_DP": 100,
    }


# ---------------------------------------------------------------------------
# Test case 1 (spec §40): shared identical mutation + strain-specific SNPs
# ---------------------------------------------------------------------------
def test_case_1_shared_identical_and_strain_specific():
    stb7a = pd.DataFrame([
        _row(456, "C", "T"),
        _row(768, "T", "C"),
        _row(987, "G", "A"),
    ])
    stb20a = pd.DataFrame([
        _row(436, "G", "A"),
        _row(565, "C", "T"),
        _row(987, "G", "A"),
    ])

    cmp_df = build_position_comparison(stb7a, stb20a, "STB7A", "STB20A")

    # Shared position 987
    row_987 = cmp_df[cmp_df["POS"] == 987].iloc[0]
    assert row_987["Comparison_Status"] == "Identical mutation"

    # STB7A-specific: 456, 768
    stb7a_specific = set(cmp_df[cmp_df["Comparison_Status"] == "STB7A-specific"]["POS"])
    assert stb7a_specific == {456, 768}

    # STB20A-specific: 436, 565
    stb20a_specific = set(cmp_df[cmp_df["Comparison_Status"] == "STB20A-specific"]["POS"])
    assert stb20a_specific == {436, 565}


# ---------------------------------------------------------------------------
# Test case 2 (spec §41): same position, different ALT -- must NOT be
# classified as identical.
# ---------------------------------------------------------------------------
def test_case_2_same_position_different_mutation():
    stb7a = pd.DataFrame([_row(987, "G", "A")])
    stb20a = pd.DataFrame([_row(987, "G", "C")])

    cmp_df = build_position_comparison(stb7a, stb20a, "STB7A", "STB20A")
    row_987 = cmp_df[cmp_df["POS"] == 987].iloc[0]

    assert row_987["Comparison_Status"] == "Same position, different mutation"
    assert row_987["Comparison_Status"] != "Identical mutation"


# ---------------------------------------------------------------------------
# Test case 3 (spec §42): same gene, different positions -- shared gene but
# NO shared genomic position.
# ---------------------------------------------------------------------------
def test_case_3_shared_gene_different_positions():
    stb7a = pd.DataFrame([_row(987, "G", "A", gene="gyrA")])
    stb20a = pd.DataFrame([_row(1200, "G", "A", gene="gyrA")])

    pos_cmp = build_position_comparison(stb7a, stb20a, "STB7A", "STB20A")
    gene_cmp = build_gene_comparison(stb7a, stb20a, "STB7A", "STB20A")

    # No shared genomic position
    shared_positions = pos_cmp[
        pos_cmp["Comparison_Status"].isin(["Identical mutation", "Same position, different mutation"])
    ]
    assert shared_positions.empty

    # But gyrA IS a shared gene
    gyra_row = gene_cmp[gene_cmp["Gene_Name"] == "gyrA"].iloc[0]
    assert gyra_row["Shared_Gene"] is True or gyra_row["Shared_Gene"] == True  # noqa: E712
    assert gyra_row["Shared_Positions"] == 0
    assert gyra_row["STB7A_Specific_Positions"] == 1
    assert gyra_row["STB20A_Specific_Positions"] == 1

    # STB7A-specific position is 987, STB20A-specific is 1200
    stb7a_specific = set(pos_cmp[pos_cmp["Comparison_Status"] == "STB7A-specific"]["POS"])
    stb20a_specific = set(pos_cmp[pos_cmp["Comparison_Status"] == "STB20A-specific"]["POS"])
    assert stb7a_specific == {987}
    assert stb20a_specific == {1200}


# ---------------------------------------------------------------------------
# Test case 4 (spec §43): one strain has a SNP, the other has none at that
# position -- must not fabricate a record for the missing strain.
# ---------------------------------------------------------------------------
def test_case_4_no_false_record_for_missing_strain():
    stb7a = pd.DataFrame([_row(456, "C", "T")])
    stb20a = pd.DataFrame([_row(999, "A", "G")])  # unrelated position

    cmp_df = build_position_comparison(stb7a, stb20a, "STB7A", "STB20A")

    row_456 = cmp_df[cmp_df["POS"] == 456].iloc[0]
    assert row_456["Comparison_Status"] == "STB7A-specific"
    # STB20A side must be empty, not a fabricated record
    assert row_456["STB20A_records"] == []
    assert row_456["STB20A_REF_alleles"] == []
    assert row_456["STB20A_ALT_alleles"] == []

    row_999 = cmp_df[cmp_df["POS"] == 999].iloc[0]
    assert row_999["Comparison_Status"] == "STB20A-specific"
    assert row_999["STB7A_records"] == []


# ---------------------------------------------------------------------------
# Reference allele discrepancy detection
# ---------------------------------------------------------------------------
def test_reference_allele_discrepancy_detected():
    stb7a = pd.DataFrame([_row(987, "G", "A")])
    stb20a = pd.DataFrame([_row(987, "A", "T")])  # different REF at same position

    cmp_df = build_position_comparison(stb7a, stb20a, "STB7A", "STB20A")
    row_987 = cmp_df[cmp_df["POS"] == 987].iloc[0]

    assert row_987["Comparison_Status"] == "Reference allele discrepancy"
    assert row_987["Comparison_Status"] != "Identical mutation"


# ---------------------------------------------------------------------------
# Original data preservation: complete records must remain retrievable
# ---------------------------------------------------------------------------
def test_original_data_preserved_in_comparison():
    stb7a = pd.DataFrame([_row(987, "G", "A", gene="gyrA", effect="missense_variant",
                                category="Missense", hgvs_p="p.Gly668Asp")])
    stb20a = pd.DataFrame([_row(987, "G", "A", gene="gyrA", effect="missense_variant",
                                 category="Missense", hgvs_p="p.Gly668Asp")])

    cmp_df = build_position_comparison(stb7a, stb20a, "STB7A", "STB20A")
    detail = get_position_detail(cmp_df, "Chromosome", 987)

    assert detail is not None
    stb7a_record = detail["STB7A_records"][0]
    stb20a_record = detail["STB20A_records"][0]

    # Full original annotation must be reachable, not reduced to POS/REF/ALT
    assert stb7a_record["Gene_Name"] == "gyrA"
    assert stb7a_record["HGVS_p"] == "p.Gly668Asp"
    assert stb7a_record["QUAL"] == 200.0
    assert stb20a_record["Gene_Name"] == "gyrA"
    assert stb20a_record["HGVS_p"] == "p.Gly668Asp"


def test_multiallelic_positions_distinguished():
    """A>G in one strain vs A>T in the other at the same position must be
    correctly distinguished, not merged into a false match."""
    stb7a = pd.DataFrame([_row(500, "A", "G")])
    stb20a = pd.DataFrame([_row(500, "A", "T")])

    cmp_df = build_position_comparison(stb7a, stb20a, "STB7A", "STB20A")
    row_500 = cmp_df[cmp_df["POS"] == 500].iloc[0]
    assert row_500["Comparison_Status"] == "Same position, different mutation"


# ---------------------------------------------------------------------------
# Real-data cross-checks against the bundled STB7A/STB36A VCFs
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not (STB7A_VCF.exists() and STB36A_VCF.exists()),
    reason="Bundled strain VCFs not present",
)
def test_real_strains_comparison_summary_is_internally_consistent():
    from pipeline import process_two_strains

    df1, df2, _, _ = process_two_strains(STB7A_VCF, STB36A_VCF, "STB7A", "STB36A")
    pos_cmp = build_position_comparison(df1, df2, "STB7A", "STB36A")
    gene_cmp = build_gene_comparison(df1, df2, "STB7A", "STB36A")
    summary = comparison_summary(pos_cmp, gene_cmp, "STB7A", "STB36A")

    # Every unique position must fall into exactly one status category.
    status_total = (
        summary["identical_mutations"]
        + summary["same_position_different_mutation"]
        + summary["reference_discrepancies"]
        + summary["STB7A_specific"]
        + summary["STB36A_specific"]
    )
    assert status_total == summary["total_unique_positions"]

    # Shared positions = identical + different-allele + ref-discrepancy
    assert summary["shared_positions"] == (
        summary["identical_mutations"]
        + summary["same_position_different_mutation"]
        + summary["reference_discrepancies"]
    )

    # No original SNPs lost: total unique positions >= either strain's own count
    assert summary["total_unique_positions"] >= df1["POS"].nunique()
    assert summary["total_unique_positions"] >= df2["POS"].nunique()


@pytest.mark.skipif(
    not (STB7A_VCF.exists() and STB36A_VCF.exists()),
    reason="Bundled strain VCFs not present",
)
def test_real_strains_original_dataframes_untouched():
    """Building the comparison must never mutate the original strain
    DataFrames (spec: individual data is never lost/modified)."""
    from pipeline import process_two_strains

    df1, df2, _, _ = process_two_strains(STB7A_VCF, STB36A_VCF, "STB7A", "STB36A")
    df1_rows_before, df2_rows_before = len(df1), len(df2)
    df1_cols_before = set(df1.columns)
    df2_cols_before = set(df2.columns)

    build_position_comparison(df1, df2, "STB7A", "STB36A")
    build_gene_comparison(df1, df2, "STB7A", "STB36A")

    assert len(df1) == df1_rows_before
    assert len(df2) == df2_rows_before
    assert set(df1.columns) == df1_cols_before
    assert set(df2.columns) == df2_cols_before
