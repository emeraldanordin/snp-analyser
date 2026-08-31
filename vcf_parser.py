"""
vcf_parser.py
-------------
Reads annotated VCF (.vcf / .vcf.gz) files using cyvcf2 and extracts every
biologically relevant field into a list of plain Python dicts (one per
sample x record combination), ready to be handed to annotation.py and
classifier.py.

This module NEVER classifies or interprets biology -- it only reads and
normalises what is literally present in the VCF. Missing fields become
None, never crash the parser, and are always counted/reported.

Design notes:
- Sample names are discovered from the VCF header (vcf.samples). Never
  hard-coded.
- Multi-allelic records (multiple ALT alleles) are expanded into one row
  per ALT allele, since each ALT can have a distinct ANN annotation and a
  distinct genotype interpretation.
- Every row keeps a `snp_id` (CHROM:POS:REF:ALT) so it stays traceable
  back to its original VCF record even after grouping/merging elsewhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import cyvcf2

from utils import get_logger, is_valid_vcf_filename

logger = get_logger(__name__)


@dataclass
class ParseReport:
    """Tracks what happened while parsing, per the data-integrity requirement
    that no VCF record is ever silently discarded."""

    filename: str = ""
    total_records: int = 0          # raw VCF lines/records read
    total_rows: int = 0             # rows after multi-allelic expansion
    annotated_rows: int = 0
    unannotated_rows: int = 0
    failed_records: int = 0
    warnings: list[str] = field(default_factory=list)
    samples: list[str] = field(default_factory=list)
    contigs: list[str] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return {
            "filename": self.filename,
            "total_records": self.total_records,
            "total_rows": self.total_rows,
            "annotated_rows": self.annotated_rows,
            "unannotated_rows": self.unannotated_rows,
            "failed_records": self.failed_records,
            "num_samples": len(self.samples),
            "samples": self.samples,
            "contigs": self.contigs,
            "warning_count": len(self.warnings),
        }


class VCFValidationError(Exception):
    """Raised when a file cannot be opened/parsed as a VCF at all."""


def validate_vcf_path(path: str | Path) -> Path:
    """
    Cheap, fast pre-check before attempting a full parse.
    Raises VCFValidationError with a human-readable message on failure.
    """
    p = Path(path)
    if not p.exists():
        raise VCFValidationError(f"File not found: {p}")
    if p.stat().st_size == 0:
        raise VCFValidationError(f"File is empty: {p}")
    if not is_valid_vcf_filename(p.name):
        raise VCFValidationError(
            f"File does not have a .vcf or .vcf.gz extension: {p.name}"
        )
    return p


def _safe_info_get(info, key: str, default=None):
    """cyvcf2's INFO object raises KeyError (not returns None) for absent
    keys in some versions -- normalise that away."""
    try:
        value = info.get(key)
        return value if value is not None else default
    except KeyError:
        return default


def _gt_string(variant, sample_index: int):
    """
    Return a human-readable genotype string (e.g. '1/1', '0/1', './.')
    for a given sample. cyvcf2's format('GT') returns raw encoded bytes,
    not usable text, so we use variant.gt_bases instead, which cyvcf2
    already renders as REF/ALT-aware strings like 'A/A' or 'A/G'.
    We keep the phased/unphased separator ('/' or '|') as reported.
    """
    try:
        bases = variant.gt_bases
        if bases is None:
            return None
        val = bases[sample_index]
        return val if val else None
    except (IndexError, TypeError, AttributeError):
        return None


def _format_value(variant, key: str, sample_index: int):
    """
    Safely pull a per-sample FORMAT field. Returns None if the field is
    absent from this VCF entirely, or not called for this sample.
    """
    try:
        arr = variant.format(key)
    except Exception:
        return None
    if arr is None:
        return None
    try:
        row = arr[sample_index]
    except (IndexError, TypeError):
        return None
    # cyvcf2 returns numpy arrays; convert to plain python types/lists
    try:
        if hasattr(row, "tolist"):
            row = row.tolist()
        if isinstance(row, (list, tuple)) and len(row) == 1:
            return row[0]
        return row
    except Exception:
        return row


def parse_vcf(path: str | Path) -> tuple[list[dict[str, Any]], ParseReport]:
    """
    Parse a VCF file into a list of row-dicts (one per sample x ALT allele)
    plus a ParseReport summarising what happened.

    Each row dict contains the raw VCF-level fields only:
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, Sample,
        INFO_* (flattened common INFO fields + a raw INFO dict),
        GT, PL, AD,
        ann_raw (the *unparsed* ANN string for this ALT, or None)

    annotation.py is responsible for turning ann_raw into structured
    columns; this function's only job is faithful, complete extraction.
    """
    path = validate_vcf_path(path)
    report = ParseReport(filename=path.name)

    try:
        vcf = cyvcf2.VCF(str(path))
    except Exception as exc:
        raise VCFValidationError(f"Could not open VCF file '{path.name}': {exc}") from exc

    report.samples = list(vcf.samples)
    if not report.samples:
        report.warnings.append(
            "No samples found in VCF header -- proceeding, but genotype "
            "fields (GT/AD/PL) will be unavailable."
        )

    try:
        report.contigs = list(vcf.seqnames)
    except Exception:
        report.warnings.append("Could not read contig/chromosome list from header.")

    rows: list[dict[str, Any]] = []

    for variant in vcf:
        report.total_records += 1
        try:
            chrom = variant.CHROM
            pos = variant.POS
            vid = variant.ID if variant.ID else "."
            ref = variant.REF
            alts = variant.ALT if variant.ALT else []
            qual = variant.QUAL
            filt = variant.FILTER  # None means PASS in cyvcf2

            if not alts:
                # No ALT allele at all -- still record it so nothing is
                # silently dropped, but flag it.
                alts = [None]
                report.warnings.append(
                    f"{chrom}:{pos} has no ALT allele; recorded with ALT=None."
                )

            info_dict = {}
            try:
                info_dict = dict(variant.INFO)
            except Exception:
                report.warnings.append(f"{chrom}:{pos} INFO field could not be read.")

            ann_raw_full = _safe_info_get(variant.INFO, "ANN", default=None)
            # ANN may list one annotation per allele (comma-separated);
            # we keep the *entire* raw string per row here and let
            # annotation.py decide which sub-annotations apply to which ALT.
            has_ann = ann_raw_full is not None and str(ann_raw_full).strip() != ""

            for alt_index, alt in enumerate(alts):
                if report.samples:
                    for sample_index, sample_name in enumerate(report.samples):
                        row = {
                            "snp_id": f"{chrom}:{pos}:{ref}:{alt}",
                            "CHROM": chrom,
                            "POS": pos,
                            "ID": vid,
                            "REF": ref,
                            "ALT": alt,
                            "ALT_INDEX": alt_index,
                            "QUAL": qual,
                            "FILTER": "PASS" if filt is None else filt,
                            "Sample": sample_name,
                            "INFO_DP": _safe_info_get(variant.INFO, "DP"),
                            "INFO_VDB": _safe_info_get(variant.INFO, "VDB"),
                            "INFO_SGB": _safe_info_get(variant.INFO, "SGB"),
                            "INFO_MQSBZ": _safe_info_get(variant.INFO, "MQSBZ"),
                            "INFO_MQ0F": _safe_info_get(variant.INFO, "MQ0F"),
                            "INFO_AC": _safe_info_get(variant.INFO, "AC"),
                            "INFO_AN": _safe_info_get(variant.INFO, "AN"),
                            "INFO_DP4": _safe_info_get(variant.INFO, "DP4"),
                            "INFO_MQ": _safe_info_get(variant.INFO, "MQ"),
                            "INFO_RAW": info_dict,
                            "GT": _gt_string(variant, sample_index),
                            "PL": _format_value(variant, "PL", sample_index),
                            "AD": _format_value(variant, "AD", sample_index),
                            "ann_raw": str(ann_raw_full) if has_ann else None,
                            "has_ann": has_ann,
                        }
                        rows.append(row)
                        report.total_rows += 1
                        if has_ann:
                            report.annotated_rows += 1
                        else:
                            report.unannotated_rows += 1
                else:
                    # No samples in file at all -- still record the site.
                    row = {
                        "snp_id": f"{chrom}:{pos}:{ref}:{alt}",
                        "CHROM": chrom,
                        "POS": pos,
                        "ID": vid,
                        "REF": ref,
                        "ALT": alt,
                        "ALT_INDEX": alt_index,
                        "QUAL": qual,
                        "FILTER": "PASS" if filt is None else filt,
                        "Sample": None,
                        "INFO_DP": _safe_info_get(variant.INFO, "DP"),
                        "INFO_VDB": _safe_info_get(variant.INFO, "VDB"),
                        "INFO_SGB": _safe_info_get(variant.INFO, "SGB"),
                        "INFO_MQSBZ": _safe_info_get(variant.INFO, "MQSBZ"),
                        "INFO_MQ0F": _safe_info_get(variant.INFO, "MQ0F"),
                        "INFO_AC": _safe_info_get(variant.INFO, "AC"),
                        "INFO_AN": _safe_info_get(variant.INFO, "AN"),
                        "INFO_DP4": _safe_info_get(variant.INFO, "DP4"),
                        "INFO_MQ": _safe_info_get(variant.INFO, "MQ"),
                        "INFO_RAW": info_dict,
                        "GT": None,
                        "PL": None,
                        "AD": None,
                        "ann_raw": str(ann_raw_full) if has_ann else None,
                        "has_ann": has_ann,
                    }
                    rows.append(row)
                    report.total_rows += 1
                    if has_ann:
                        report.annotated_rows += 1
                    else:
                        report.unannotated_rows += 1

        except Exception as exc:
            report.failed_records += 1
            msg = f"Failed to parse record #{report.total_records} ({exc})"
            report.warnings.append(msg)
            logger.warning(msg)
            continue

    logger.info(
        "Parsed %s: %d records -> %d rows (%d annotated, %d unannotated, %d failed)",
        report.filename,
        report.total_records,
        report.total_rows,
        report.annotated_rows,
        report.unannotated_rows,
        report.failed_records,
    )

    return rows, report
