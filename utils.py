"""
utils.py
--------
Shared, dependency-light helpers used across the SNP Annotation Analyser:
logging setup and common filesystem path utilities.

Keeping these in one place avoids duplicated logging configuration across
vcf_parser.py, annotation.py, classifier.py, grouping.py, statistics.py,
comparison.py, report.py, and export.py.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path


# Project-level paths (never hard-code paths elsewhere in the app)

PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "output"
TEMPLATES_DIR = PROJECT_ROOT / "templates"
STATIC_DIR = PROJECT_ROOT / "static"

# Accepted VCF file extensions (used by the uploader / CLI validators later)
VALID_VCF_SUFFIXES = (".vcf", ".vcf.gz")


def get_logger(name: str, log_to_file: bool = True) -> logging.Logger:
    """
    Return a configured logger.

    Every module in the pipeline calls `get_logger(__name__)` rather than
    configuring logging itself. This guarantees one consistent log format
    and prevents duplicate handlers being attached if a module is imported
    more than once (e.g. by both app.py and main.py in the same session).

    Parameters
    ----------
    name : str
        Usually __name__ of the calling module.
    log_to_file : bool
        If True, also write logs to output/snp_analyzer.log so a run's
        full processing history is available after the fact (important for
        the "never silently discard a record" requirement -- warnings and
        errors must be traceable after a Streamlit session closes).

    Returns
    -------
    logging.Logger
    """
    logger = logging.getLogger(name)

    if logger.handlers:
        # Already configured (e.g. re-imported); avoid duplicate handlers.
        return logger

    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler: INFO and above, so normal runs aren't noisy.
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if log_to_file:
        try:
            OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(
                OUTPUT_DIR / "snp_analyzer.log", encoding="utf-8"
            )
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        except OSError:
            # If we can't write a log file (e.g. read-only filesystem),
            # fall back to console-only logging rather than crashing.
            logger.warning(
                "Could not create log file in %s; continuing with console logging only.",
                OUTPUT_DIR,
            )

    return logger


def ensure_output_dir(path: Path | str | None = None) -> Path:
    """
    Ensure the given output directory exists and return it as a Path.
    Defaults to the project's `output/` directory if none is given.
    Used by main.py (CLI) and export.py so output locations are never
    hard-coded in more than one place.
    """
    out_path = Path(path) if path is not None else OUTPUT_DIR
    out_path.mkdir(parents=True, exist_ok=True)
    return out_path


def ensure_strain_output_dirs(base_output: Path | str, strain1_name: str, strain2_name: str) -> dict[str, Path]:
    """
    Create (and return) the three-way output layout used by two-strain
    comparison runs:
        <base_output>/<strain1_name>/
        <base_output>/<strain2_name>/
        <base_output>/comparison/
    Individual strain reports/exports go in their own subfolder so
    STB7A's and STB20A's files never collide or overwrite each other,
    and the comparison outputs live in a clearly separate location.
    """
    base = ensure_output_dir(base_output)
    dirs = {
        "strain1": ensure_output_dir(base / strain1_name),
        "strain2": ensure_output_dir(base / strain2_name),
        "comparison": ensure_output_dir(base / "comparison"),
    }
    return dirs


def is_valid_vcf_filename(filename: str) -> bool:
    """
    Basic filename-level check that a file looks like a VCF.
    This is NOT full VCF validation (that happens in vcf_parser.py by
    actually attempting to open the file with cyvcf2) -- it's a fast,
    cheap first check used by the Streamlit uploader and the CLI to give
    an immediate, friendly error before any parsing is attempted.
    """
    lowered = filename.lower()
    return lowered.endswith(".vcf") or lowered.endswith(".vcf.gz")
