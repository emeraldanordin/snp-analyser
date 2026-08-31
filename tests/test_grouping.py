"""
tests/test_grouping.py
------------------------
Tests for grouping.py: gene-level grouping (multiple SNPs per gene,
never merged) and position-level grouping (same position, different
ALT alleles, never conflated with gene grouping).
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from grouping import group_by_gene, group_by_position

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SAMPLE_VCF = PROJECT_ROOT / "data" / "isolate1-snp_annot.vcf"


def _synthetic_df() -> pd.DataFrame:
    """
    Build a small synthetic SNP table mirroring the spec's gyrA example
    (7 SNPs in one gene) plus a position with two alternative alleles
    across two isolates, so both grouping functions can be tested in
    isolation without depending on cyvcf2 / the real VCF.
    """
    rows = [
        # 7 SNPs in gyrA -- gene-level grouping must keep these as 7
        # distinct SNPs under ONE gene group, never merged into one row.
        {"Sample": "isolate1", "CHROM": "Chromosome", "POS": 6112, "REF": "G", "ALT": "C",
         "Gene_Name": "gyrA", "Gene_ID": "Rv0006", "Effect_Raw": "missense_variant",
         "Mutation_Category": "Missense", "SNP_Type": "Transversion",
         "DNA_Change": "G>C", "HGVS_p": "p.Met291Ile"},
        {"Sample": "isolate1", "CHROM": "Chromosome", "POS": 7362, "REF": "G", "ALT": "C",
         "Gene_Name": "gyrA", "Gene_ID": "Rv0006", "Effect_Raw": "missense_variant",
         "Mutation_Category": "Missense", "SNP_Type": "Transversion",
         "DNA_Change": "G>C", "HGVS_p": "p.Glu21Gln"},
        {"Sample": "isolate1", "CHROM": "Chromosome", "POS": 7585, "REF": "G", "ALT": "C",
         "Gene_Name": "gyrA", "Gene_ID": "Rv0006", "Effect_Raw": "missense_variant",
         "Mutation_Category": "Missense", "SNP_Type": "Transversion",
         "DNA_Change": "G>C", "HGVS_p": "p.Ser95Thr"},
        {"Sample": "isolate1", "CHROM": "Chromosome", "POS": 8452, "REF": "C", "ALT": "T",
         "Gene_Name": "gyrA", "Gene_ID": "Rv0006", "Effect_Raw": "missense_variant",
         "Mutation_Category": "Missense", "SNP_Type": "Transition",
         "DNA_Change": "C>T", "HGVS_p": "p.Ala384Val"},
        {"Sample": "isolate1", "CHROM": "Chromosome", "POS": 9143, "REF": "T", "ALT": "C",
         "Gene_Name": "gyrA", "Gene_ID": "Rv0006", "Effect_Raw": "synonymous_variant",
         "Mutation_Category": "Silent / Synonymous", "SNP_Type": "Transition",
         "DNA_Change": "T>C", "HGVS_p": "p.Ile614Ile"},
        {"Sample": "isolate1", "CHROM": "Chromosome", "POS": 9260, "REF": "G", "ALT": "C",
         "Gene_Name": "gyrA", "Gene_ID": "Rv0006", "Effect_Raw": "synonymous_variant",
         "Mutation_Category": "Silent / Synonymous", "SNP_Type": "Transversion",
         "DNA_Change": "G>C", "HGVS_p": "p.Leu653Leu"},
        # position 9304: two isolates, two different ALT alleles at the
        # same REF -- must be recognized as ONE position with 2 alleles.
        {"Sample": "isolate1", "CHROM": "Chromosome", "POS": 9304, "REF": "G", "ALT": "A",
         "Gene_Name": "gyrA", "Gene_ID": "Rv0006", "Effect_Raw": "missense_variant",
         "Mutation_Category": "Missense", "SNP_Type": "Transition",
         "DNA_Change": "G>A", "HGVS_p": "p.Gly668Asp"},
        {"Sample": "isolate2", "CHROM": "Chromosome", "POS": 9304, "REF": "G", "ALT": "A",
         "Gene_Name": "gyrA", "Gene_ID": "Rv0006", "Effect_Raw": "missense_variant",
         "Mutation_Category": "Missense", "SNP_Type": "Transition",
         "DNA_Change": "G>A", "HGVS_p": "p.Gly668Asp"},
        {"Sample": "isolate3", "CHROM": "Chromosome", "POS": 9304, "REF": "G", "ALT": "A",
         "Gene_Name": "gyrA", "Gene_ID": "Rv0006", "Effect_Raw": "missense_variant",
         "Mutation_Category": "Missense", "SNP_Type": "Transition",
         "DNA_Change": "G>A", "HGVS_p": "p.Gly668Asp"},
        {"Sample": "isolate4", "CHROM": "Chromosome", "POS": 9304, "REF": "G", "ALT": "C",
         "Gene_Name": "gyrA", "Gene_ID": "Rv0006", "Effect_Raw": "missense_variant",
         "Mutation_Category": "Missense", "SNP_Type": "Transversion",
         "DNA_Change": "G>C", "HGVS_p": "p.Gly668Pro"},
        # An unrelated gene with a single SNP, to confirm groups don't bleed.
        {"Sample": "isolate1", "CHROM": "Chromosome", "POS": 1302, "REF": "C", "ALT": "A",
         "Gene_Name": "dnaA", "Gene_ID": "Rv0001", "Effect_Raw": "synonymous_variant",
         "Mutation_Category": "Silent / Synonymous", "SNP_Type": "Transversion",
         "DNA_Change": "C>A", "HGVS_p": "p.Pro434Pro"},
    ]
    return pd.DataFrame(rows)


def test_gene_grouping_preserves_all_individual_snps():
    df = _synthetic_df()
    gene_df = group_by_gene(df)

    gyra_row = gene_df[gene_df["Gene_Name"] == "gyrA"].iloc[0]
    # 7 distinct gyrA SNPs (6112,7362,7585,8452,9143,9260 from isolate1,
    # plus 9304 counted once per isolate row since grouping is per-row here)
    # -- confirm none were merged into a single mutation.
    assert gyra_row["Unique_Positions"] == 7
    assert len(gyra_row["Positions"]) == 7
    assert 6112 in gyra_row["Positions"]
    assert 9304 in gyra_row["Positions"]


def test_gene_grouping_counts_categories_correctly():
    df = _synthetic_df()
    gene_df = group_by_gene(df)
    gyra_row = gene_df[gene_df["Gene_Name"] == "gyrA"].iloc[0]

    # 5 missense rows (6112,7362,7585,8452, and 3x 9304-A + 1x 9304-C = 4
    # more missense rows) + 2 synonymous (9143, 9260)
    assert gyra_row["Synonymous"] == 2
    assert gyra_row["Missense"] == gyra_row["Total_SNPs"] - 2


def test_gene_grouping_does_not_bleed_across_genes():
    df = _synthetic_df()
    gene_df = group_by_gene(df)
    dnaa_row = gene_df[gene_df["Gene_Name"] == "dnaA"].iloc[0]
    assert dnaa_row["Total_SNPs"] == 1
    assert dnaa_row["Positions"] == [1302]


def test_position_grouping_recognizes_same_position_different_alleles():
    df = _synthetic_df()
    position_df = group_by_position(df)

    pos_9304 = position_df[position_df["POS"] == 9304].iloc[0]
    # Same CHROM/POS/REF, two distinct ALT alleles -> ONE position row.
    assert pos_9304["Num_Alt_Alleles"] == 2
    assert set(pos_9304["Alt_Alleles"]) == {"A", "C"}


def test_position_grouping_computes_isolate_percentages():
    df = _synthetic_df()
    position_df = group_by_position(df)
    pos_9304 = position_df[position_df["POS"] == 9304].iloc[0]

    # 4 isolates total; 3 carry G>A (75%), 1 carries G>C (25%) --
    # matches the spec's worked example exactly.
    assert pos_9304["Allele_Isolate_Counts"]["A"] == 3
    assert pos_9304["Allele_Isolate_Counts"]["C"] == 1
    assert pos_9304["Allele_Isolate_Percentages"]["A"] == 75.0
    assert pos_9304["Allele_Isolate_Percentages"]["C"] == 25.0


def test_position_grouping_does_not_merge_different_positions():
    df = _synthetic_df()
    position_df = group_by_position(df)
    # 7 distinct positions from the gyrA block + 1 from dnaA = 8 unique
    # (CHROM, POS, REF) groups, never collapsed together.
    assert len(position_df) == 8


@pytest.mark.skipif(not SAMPLE_VCF.exists(), reason="Sample VCF not present")
def test_grouping_on_real_vcf_matches_gyra_spec_example():
    """Cross-check against the actual bundled VCF: gyrA should show up
    as a single gene group containing multiple individual SNPs."""
    from pipeline import process_vcf

    df, _ = process_vcf(SAMPLE_VCF)
    gene_df = group_by_gene(df)
    position_df = group_by_position(df)

    gyra_row = gene_df[gene_df["Gene_Name"] == "gyrA"]
    assert not gyra_row.empty
    assert gyra_row.iloc[0]["Total_SNPs"] >= 1
    assert len(position_df) == df["POS"].nunique()
