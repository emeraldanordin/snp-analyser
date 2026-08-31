"""
annotation.py
-------------
Parses the SnpEff ANN INFO field into structured columns.

ANN format (per the VCF header definition SnpEff writes):
    Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID |
    Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c |
    HGVS.p | cDNA.pos/cDNA.length | CDS.pos/CDS.length | AA.pos/AA.length |
    Distance | ERRORS/WARNINGS/INFO

A record's ANN field can contain multiple comma-separated annotations
(multiple transcripts/genes affected, or one annotation per ALT allele).
This module:
  1. Splits on ',' to get all individual annotations.
  2. Splits each on '|' into the 16 named sub-fields.
  3. Filters to annotations matching the row's specific ALT allele
     (SnpEff's "Allele" sub-field), since a multi-allelic site's ANN
     block contains annotations for every ALT allele.
  4. Picks a single "primary" annotation to promote into flat columns
     for the main table, while retaining the full list for later display.
"""

from __future__ import annotations

from typing import Any

from utils import get_logger

logger = get_logger(__name__)

ANN_FIELDS = [
    "Allele",
    "Effect_Raw",
    "Impact",
    "Gene_Name",
    "Gene_ID",
    "Feature_Type",
    "Feature_ID",
    "Transcript_BioType",
    "Rank",
    "HGVS_c",
    "HGVS_p",
    "cDNA_pos_len",
    "CDS_pos_len",
    "AA_pos_len",
    "Distance",
    "Warnings",
]

NOT_ANNOTATED = "Not annotated"

# SnpEff impact ranks, used to choose the "most severe" annotation as
# primary when a record maps to multiple genes/transcripts.
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}


def _empty_annotation() -> dict[str, Any]:
    return {
        "Allele": None,
        "Effect_Raw": NOT_ANNOTATED,
        "Impact": "Unknown",
        "Gene_Name": NOT_ANNOTATED,
        "Gene_ID": NOT_ANNOTATED,
        "Feature_Type": NOT_ANNOTATED,
        "Feature_ID": NOT_ANNOTATED,
        "Transcript_BioType": NOT_ANNOTATED,
        "Rank": NOT_ANNOTATED,
        "HGVS_c": NOT_ANNOTATED,
        "HGVS_p": NOT_ANNOTATED,
        "cDNA_pos_len": NOT_ANNOTATED,
        "CDS_pos_len": NOT_ANNOTATED,
        "AA_pos_len": NOT_ANNOTATED,
        "Distance": NOT_ANNOTATED,
        "Warnings": "",
    }


def parse_single_ann_entry(entry: str) -> dict[str, Any]:
    """Split one '|'-delimited ANN sub-record into a field dict."""
    parts = entry.split("|")
    while len(parts) < len(ANN_FIELDS):
        parts.append("")
    parsed = dict(zip(ANN_FIELDS, parts))
    # Normalize blanks to a consistent "not present" representation
    # without inventing biological meaning for genuinely empty fields.
    for key in ("HGVS_p", "Distance", "Rank", "Gene_Name", "Gene_ID"):
        if parsed.get(key, "") == "":
            parsed[key] = NOT_ANNOTATED if key in ("Gene_Name", "Gene_ID") else ""
    return parsed


def parse_ann_field(ann_raw: str | None, alt_allele: str | None = None) -> list[dict[str, Any]]:
    """
    Parse a raw ANN string into a list of annotation dicts.

    If alt_allele is provided, annotations are filtered to those whose
    'Allele' sub-field matches it (handles multi-allelic sites where the
    ANN block covers more than one ALT). If none match (or alt_allele is
    None), all parsed annotations are returned unfiltered.
    """
    if not ann_raw or not str(ann_raw).strip():
        return []

    entries = [e for e in str(ann_raw).split(",") if e.strip()]
    parsed = []
    for entry in entries:
        try:
            parsed.append(parse_single_ann_entry(entry))
        except Exception as exc:
            logger.warning("Could not parse ANN sub-entry '%s': %s", entry, exc)
            continue

    if alt_allele:
        matching = [p for p in parsed if p.get("Allele") == alt_allele]
        if matching:
            return matching
    return parsed


def choose_primary_annotation(annotations: list[dict[str, Any]]) -> dict[str, Any]:
    """
    Choose which annotation to promote to the main flat columns when a
    record has multiple. Defensible rule: prefer the highest-impact
    annotation (HIGH > MODERATE > LOW > MODIFIER); ties broken by order
    of appearance (SnpEff already ranks by relevance internally).
    """
    if not annotations:
        return _empty_annotation()
    return sorted(
        annotations, key=lambda a: IMPACT_RANK.get(a.get("Impact", ""), 99)
    )[0]


def annotate_row(row: dict[str, Any]) -> dict[str, Any]:
    """
    Given one parsed VCF row (from vcf_parser.parse_vcf), return a dict of
    annotation-derived columns to merge into it: the flat primary-annotation
    columns, plus `all_annotations` (full list) and `has_multiple_annotations`.
    """
    ann_raw = row.get("ann_raw")
    alt = row.get("ALT")

    annotations = parse_ann_field(ann_raw, alt_allele=alt)
    primary = choose_primary_annotation(annotations)

    return {
        "Effect_Raw": primary["Effect_Raw"],
        "Impact": primary["Impact"],
        "Gene_Name": primary["Gene_Name"],
        "Gene_ID": primary["Gene_ID"],
        "Feature_Type": primary["Feature_Type"],
        "Feature_ID": primary["Feature_ID"],
        "Transcript_BioType": primary["Transcript_BioType"],
        "Rank": primary["Rank"],
        "HGVS_c": primary["HGVS_c"],
        "HGVS_p": primary["HGVS_p"],
        "cDNA_pos_len": primary["cDNA_pos_len"],
        "CDS_pos_len": primary["CDS_pos_len"],
        "AA_pos_len": primary["AA_pos_len"],
        "Distance": primary["Distance"],
        "Warnings": primary["Warnings"],
        "all_annotations": annotations,
        "has_multiple_annotations": len(annotations) > 1,
        "num_annotations": len(annotations),
    }


def annotate_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Apply annotate_row to every row and merge the results in-place."""
    for row in rows:
        ann_cols = annotate_row(row)
        row.update(ann_cols)
    return rows
