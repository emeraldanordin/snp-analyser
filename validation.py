"""
validation.py
--------------
Cross-VCF compatibility checks used before comparing two strains.

This module never blocks a comparison from running -- per the spec,
the application must warn rather than silently merge or refuse
incompatible datasets. It returns structured warnings that the caller
(main.py / app.py) is responsible for displaying prominently.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from utils import get_logger
from vcf_parser import ParseReport

logger = get_logger(__name__)


@dataclass
class CompatibilityReport:
    """Result of comparing two strains' VCF metadata for reference-genome
    compatibility before running coordinate-based comparison."""

    compatible: bool = True
    warnings: list[str] = field(default_factory=list)
    strain1_contigs: list[str] = field(default_factory=list)
    strain2_contigs: list[str] = field(default_factory=list)
    shared_contigs: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "compatible": self.compatible,
            "warnings": self.warnings,
            "strain1_contigs": self.strain1_contigs,
            "strain2_contigs": self.strain2_contigs,
            "shared_contigs": self.shared_contigs,
        }


def check_reference_compatibility(
    parse_report1: ParseReport,
    parse_report2: ParseReport,
    strain1_name: str,
    strain2_name: str,
) -> CompatibilityReport:
    """
    Compare contig/chromosome names between two parsed VCFs and warn if
    they look like they may come from different reference assemblies.

    This is a heuristic, name-based check (VCF headers don't reliably
    encode assembly identity) -- it is intentionally conservative: it
    warns rather than blocks, since a real mismatch requires the user's
    judgment, not an automatic decision.
    """
    report = CompatibilityReport()
    report.strain1_contigs = sorted(parse_report1.contigs)
    report.strain2_contigs = sorted(parse_report2.contigs)

    set1 = set(report.strain1_contigs)
    set2 = set(report.strain2_contigs)
    report.shared_contigs = sorted(set1 & set2)

    if not set1 or not set2:
        report.warnings.append(
            "Could not determine contig/chromosome names for one or both "
            "strains; reference compatibility could not be verified."
        )
        report.compatible = False
        return report

    if set1 != set2:
        report.compatible = False
        only1 = sorted(set1 - set2)
        only2 = sorted(set2 - set1)
        detail = []
        if only1:
            detail.append(f"{strain1_name} only: {', '.join(only1)}")
        if only2:
            detail.append(f"{strain2_name} only: {', '.join(only2)}")
        report.warnings.append(
            "These VCF files may have been generated against different "
            "reference assemblies. Position-based comparison may be "
            "unreliable. " + "; ".join(detail)
        )

    logger.info(
        "Reference compatibility check (%s vs %s): compatible=%s, shared_contigs=%d",
        strain1_name, strain2_name, report.compatible, len(report.shared_contigs),
    )
    return report


def check_reference_allele_consistency(
    comparison_df: pd.DataFrame,
) -> dict:
    """
    Summarise how many shared positions show a REF allele mismatch
    between strains (Comparison_Status == 'Reference allele discrepancy').
    This is a lightweight summary helper for the comparison report/UI;
    the actual per-position detection happens in comparison.py.
    """
    if comparison_df.empty or "Comparison_Status" not in comparison_df.columns:
        return {"reference_discrepancies": 0, "positions": []}

    mismatches = comparison_df[comparison_df["Comparison_Status"] == "Reference allele discrepancy"]
    return {
        "reference_discrepancies": len(mismatches),
        "positions": mismatches["POS"].tolist() if not mismatches.empty else [],
    }
