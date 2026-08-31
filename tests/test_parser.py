"""
tests/test_parser.py
---------------------
Tests for vcf_parser.py. Uses the bundled data/isolate1-snp_annot.vcf
as a real-world fixture, plus a synthetic minimal VCF for edge cases.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from vcf_parser import VCFValidationError, parse_vcf, validate_vcf_path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SAMPLE_VCF = PROJECT_ROOT / "data" / "isolate1-snp_annot.vcf"


def test_validate_missing_file(tmp_path):
    with pytest.raises(VCFValidationError):
        validate_vcf_path(tmp_path / "does_not_exist.vcf")


def test_validate_empty_file(tmp_path):
    empty_vcf = tmp_path / "empty.vcf"
    empty_vcf.write_text("")
    with pytest.raises(VCFValidationError):
        validate_vcf_path(empty_vcf)


def test_validate_wrong_extension(tmp_path):
    bad_file = tmp_path / "not_a_vcf.txt"
    bad_file.write_text("some content")
    with pytest.raises(VCFValidationError):
        validate_vcf_path(bad_file)


@pytest.mark.skipif(not SAMPLE_VCF.exists(), reason="Sample VCF not present")
def test_parse_real_vcf_no_crash():
    rows, report = parse_vcf(SAMPLE_VCF)
    assert report.total_records > 0
    assert report.failed_records == 0
    assert len(rows) == report.total_rows
    assert "isolate1" in report.samples


@pytest.mark.skipif(not SAMPLE_VCF.exists(), reason="Sample VCF not present")
def test_parse_real_vcf_has_annotated_and_unannotated():
    rows, report = parse_vcf(SAMPLE_VCF)
    assert report.annotated_rows > 0
    assert report.unannotated_rows > 0
    assert report.annotated_rows + report.unannotated_rows == report.total_rows


@pytest.mark.skipif(not SAMPLE_VCF.exists(), reason="Sample VCF not present")
def test_missing_ann_does_not_crash():
    """Position 1977 in the sample VCF has no ANN field per the spec example."""
    rows, report = parse_vcf(SAMPLE_VCF)
    row_1977 = [r for r in rows if r["POS"] == 1977]
    assert len(row_1977) == 1
    assert row_1977[0]["has_ann"] is False
    assert row_1977[0]["ann_raw"] is None


def test_parse_minimal_synthetic_vcf(tmp_path):
    """A hand-built minimal VCF with one annotated and one unannotated record."""
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=Chromosome,length=1000000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=ANN,Number=.,Type=String,Description="ann">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttest_sample
Chromosome\t100\t.\tC\tA\t200\tPASS\tDP=50;ANN=A|synonymous_variant|LOW|geneX|GX001|transcript|T1|protein_coding|1/1|c.1C>A|p.Pro1Pro|1/100|1/100|1/33||\tGT:AD\t1/1:0,50
Chromosome\t200\t.\tG\tT\t150\tPASS\tDP=40\tGT:AD\t1/1:0,40
"""
    vcf_path = tmp_path / "synthetic.vcf"
    vcf_path.write_text(vcf_content)

    rows, report = parse_vcf(vcf_path)
    assert report.total_records == 2
    assert report.annotated_rows == 1
    assert report.unannotated_rows == 1
    assert report.failed_records == 0
    assert report.samples == ["test_sample"]

    annotated_row = [r for r in rows if r["POS"] == 100][0]
    assert annotated_row["has_ann"] is True

    unannotated_row = [r for r in rows if r["POS"] == 200][0]
    assert unannotated_row["has_ann"] is False
