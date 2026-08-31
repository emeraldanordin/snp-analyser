"""
report.py
---------
Generates a standalone, self-contained HTML report (no internet
dependency -- Plotly JS is embedded inline) and a publication-friendly
PDF version via WeasyPrint.

WeasyPrint cannot execute JavaScript, so for the PDF path the same
Plotly figures are rendered as static PNGs (via kaleido) and swapped in
before conversion, while the HTML report keeps the fully interactive
Plotly.js version.
"""

from __future__ import annotations

import base64
from datetime import datetime
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader

from utils import STATIC_DIR, TEMPLATES_DIR, ensure_output_dir, get_logger
from vcf_parser import ParseReport

logger = get_logger(__name__)

# Columns to display in the report's SNP table -- kept to the most useful
# subset so the report table isn't unreadably wide; the full data is still
# available via CSV/Excel export.
SNP_TABLE_DISPLAY_COLUMNS = [
    "Sample", "CHROM", "POS", "REF", "ALT", "DNA_Change", "QUAL", "FILTER",
    "INFO_DP", "GT", "AD", "Effect_Raw", "Impact", "Mutation_Category",
    "Gene_Name", "Gene_ID", "Feature_ID", "Transcript_BioType",
    "HGVS_c", "HGVS_p", "AA_pos_len", "SNP_Type",
]

GENE_TABLE_DISPLAY_COLUMNS = [
    "Gene_Name", "Gene_ID", "Total_SNPs", "Unique_Positions", "Missense",
    "Synonymous", "Nonsense", "Frameshift", "Other_Effects",
    "Transitions", "Transversions", "Ti_Tv_Ratio",
]

POSITION_TABLE_DISPLAY_COLUMNS = [
    "CHROM", "POS", "REF", "Num_Alt_Alleles", "Alt_Alleles",
    "Gene_Name", "Effect_Raw", "Mutation_Category",
]


def _df_to_html_table(df: pd.DataFrame, columns: list[str] | None = None, max_rows: int | None = None) -> str:
    """Render a DataFrame as an HTML table, restricted to display columns
    that actually exist, with nested list/dict cells stringified."""
    if df.empty:
        return "<p><em>No data available.</em></p>"

    display_df = df.copy()
    if columns:
        available = [c for c in columns if c in display_df.columns]
        display_df = display_df[available]

    for col in display_df.columns:
        if display_df[col].apply(lambda v: isinstance(v, (list, dict))).any():
            display_df[col] = display_df[col].apply(
                lambda v: "; ".join(map(str, v)) if isinstance(v, list)
                else "; ".join(f"{k}={val}" for k, val in v.items()) if isinstance(v, dict)
                else v
            )

    truncated_note = ""
    if max_rows and len(display_df) > max_rows:
        truncated_note = (
            f"<p><em>Showing first {max_rows} of {len(display_df)} rows. "
            f"Export CSV/Excel for the complete table.</em></p>"
        )
        display_df = display_df.head(max_rows)

    table_html = display_df.to_html(
        index=False, classes="wide-table", border=0, na_rep="-", escape=True
    )
    return truncated_note + table_html


def build_report_context(
    vcf_filename: str,
    snp_df: pd.DataFrame,
    gene_df: pd.DataFrame,
    position_df: pd.DataFrame,
    stats: dict,
    parse_report: ParseReport,
    figures_html: dict[str, str],
    max_snp_table_rows: int | None = 2000,
) -> dict:
    """Assemble the full Jinja2 template context in one place, so both
    the HTML and PDF paths build the report from identical data."""
    css_content = (STATIC_DIR / "style.css").read_text(encoding="utf-8")

    samples = sorted(snp_df["Sample"].dropna().unique().tolist()) if "Sample" in snp_df else []
    contigs = sorted(snp_df["CHROM"].dropna().unique().tolist()) if "CHROM" in snp_df else []

    return {
        "vcf_filename": vcf_filename,
        "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "samples": samples,
        "contigs": contigs,
        "stats": stats,
        "parse_report": parse_report.as_dict(),
        "css_content": css_content,
        "figures": figures_html,
        "gene_table_html": _df_to_html_table(gene_df, GENE_TABLE_DISPLAY_COLUMNS),
        "gene_table_note": f"{len(gene_df)} genes/groups with at least one SNP.",
        "position_table_html": _df_to_html_table(position_df, POSITION_TABLE_DISPLAY_COLUMNS, max_rows=2000),
        "position_table_note": f"{len(position_df)} unique genomic positions.",
        "snp_table_html": _df_to_html_table(snp_df, SNP_TABLE_DISPLAY_COLUMNS, max_rows=max_snp_table_rows),
        "snp_table_note": f"{len(snp_df)} total SNPs.",
    }


def render_html_report(context: dict, output_path: str | Path) -> Path:
    """Render templates/report.html with the given context and write it
    to disk as a standalone HTML file."""
    output_path = Path(output_path)
    ensure_output_dir(output_path.parent)

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)), autoescape=False)
    template = env.get_template("report.html")
    html_content = template.render(**context)

    output_path.write_text(html_content, encoding="utf-8")
    logger.info("Wrote HTML report: %s", output_path)
    return output_path


def render_pdf_report(context: dict, output_path: str | Path, figures_static: dict[str, "any"] | None = None) -> Path:
    """
    Render a PDF version of the report using WeasyPrint. Since WeasyPrint
    cannot execute the interactive Plotly JS, `context['figures']` should
    already contain static <img> HTML (see report.generate_static_figure_html)
    rather than the interactive Plotly.to_html() output used for the HTML
    report.
    """
    try:
        from weasyprint import HTML
    except (ImportError, OSError) as exc:
        # ImportError: weasyprint itself isn't installed.
        # OSError: weasyprint IS installed but can't load its native
        # dependencies (Pango/cairo/gdk-pixbuf/libgobject) at import time --
        # this is the common case on a fresh macOS setup before `brew
        # install pango` has been run. Either way, treat it as "PDF isn't
        # available right now" rather than letting it crash the whole run.
        raise RuntimeError(
            "WeasyPrint could not load its native dependencies "
            "(Pango/cairo/gdk-pixbuf/libgobject) -- PDF generation skipped. "
            "See README.md for the macOS installation steps "
            "(brew install pango cairo gdk-pixbuf libffi)."
        ) from exc

    output_path = Path(output_path)
    ensure_output_dir(output_path.parent)

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)), autoescape=False)
    template = env.get_template("report.html")
    html_content = template.render(**context)

    try:
        HTML(string=html_content, base_url=str(TEMPLATES_DIR)).write_pdf(str(output_path))
    except Exception as exc:
        # Any failure during actual PDF rendering (bad font config, layout
        # error, etc.) should also degrade gracefully rather than take
        # down the whole batch/comparison run.
        raise RuntimeError(f"PDF rendering failed: {exc}") from exc

    logger.info("Wrote PDF report: %s", output_path)
    return output_path



def figures_to_interactive_html(figures: dict[str, "any"]) -> dict[str, str]:
    """Convert a dict of Plotly figures to embedded, JS-inline HTML
    fragments for the interactive HTML report."""
    html_fragments = {}
    include_js = "inline"
    for i, (name, fig) in enumerate(figures.items()):
        html_fragments[name] = fig.to_html(
            full_html=False,
            include_plotlyjs=(include_js if i == 0 else False),
            div_id=f"fig-{name}",
        )
    return html_fragments


def figures_to_static_html(figures: dict[str, "any"]) -> dict[str, str]:
    """
    Convert a dict of Plotly figures to static base64-embedded <img> tags
    for the PDF report (WeasyPrint can't execute Plotly's JS). Requires
    kaleido to be installed for fig.to_image().
    """
    html_fragments = {}
    for name, fig in figures.items():
        try:
            png_bytes = fig.to_image(format="png", scale=2)
            b64 = base64.b64encode(png_bytes).decode("ascii")
            html_fragments[name] = f'<img src="data:image/png;base64,{b64}" style="max-width:100%;">'
        except Exception as exc:
            logger.warning("Could not render static image for figure '%s': %s", name, exc)
            html_fragments[name] = f"<p><em>Figure '{name}' could not be rendered for PDF.</em></p>"
    return html_fragments


def generate_reports(
    vcf_filename: str,
    snp_df: pd.DataFrame,
    gene_df: pd.DataFrame,
    position_df: pd.DataFrame,
    stats: dict,
    parse_report: ParseReport,
    figures: dict[str, "any"],
    output_dir: str | Path,
    make_pdf: bool = True,
) -> dict[str, Path]:
    """
    High-level convenience function used by both main.py and app.py:
    builds the shared context, writes report.html (interactive), and
    optionally report.pdf (static images). Returns paths written.
    """
    out = ensure_output_dir(output_dir)
    results = {}

    interactive_figs = figures_to_interactive_html(figures)
    html_context = build_report_context(
        vcf_filename, snp_df, gene_df, position_df, stats, parse_report, interactive_figs
    )
    results["html"] = render_html_report(html_context, out / "report.html")

    if make_pdf:
        try:
            static_figs = figures_to_static_html(figures)
            pdf_context = build_report_context(
                vcf_filename, snp_df, gene_df, position_df, stats, parse_report, static_figs,
                max_snp_table_rows=500,  # keep PDF a sane size
            )
            results["pdf"] = render_pdf_report(pdf_context, out / "report.pdf")
        except RuntimeError as exc:
            logger.warning("Skipping PDF generation: %s", exc)

    return results


# TWO-STRAIN COMPARISON REPORT

COMPARISON_GENE_DISPLAY_LIMIT = 500
COMPARISON_POSITION_DISPLAY_LIMIT = 2000
COMPARISON_SPECIFIC_DISPLAY_LIMIT = 1000


def _comparison_df_for_display(df: pd.DataFrame, strain1_name: str, strain2_name: str) -> pd.DataFrame:
    """Prepare the position comparison table for HTML rendering: adds a
    readable '{strain}_DNA_Change' column (e.g. 'G>A') per strain, and
    drops the bulky '{strain}_records' list-of-dict columns plus the raw
    REF/ALT allele-list columns (superseded by DNA_Change for display).
    The full records and raw allele lists remain in the CSV export and
    in the in-memory DataFrame for the detail view."""
    from comparison import add_dna_change_columns

    display = add_dna_change_columns(df, [strain1_name, strain2_name])
    drop_cols = [
        f"{strain1_name}_records", f"{strain2_name}_records",
        f"{strain1_name}_REF_alleles", f"{strain1_name}_ALT_alleles",
        f"{strain2_name}_REF_alleles", f"{strain2_name}_ALT_alleles",
    ]
    return display.drop(columns=[c for c in drop_cols if c in display.columns])


def build_comparison_report_context(
    strain1_name: str,
    strain2_name: str,
    vcf1_filename: str,
    vcf2_filename: str,
    stats1: dict,
    stats2: dict,
    position_comparison_df: pd.DataFrame,
    gene_comparison_df: pd.DataFrame,
    summary: dict,
    figures_html: dict[str, str],
    compatibility_warnings: list[str] | None = None,
) -> dict:
    """Assemble the Jinja2 context for the comparison report template."""
    css_content = (STATIC_DIR / "style.css").read_text(encoding="utf-8")

    display_df = _comparison_df_for_display(position_comparison_df, strain1_name, strain2_name)
    strain1_specific = display_df[display_df["Comparison_Status"] == f"{strain1_name}-specific"]
    strain2_specific = display_df[display_df["Comparison_Status"] == f"{strain2_name}-specific"]

    return {
        "strain1_name": strain1_name,
        "strain2_name": strain2_name,
        "vcf1_filename": vcf1_filename,
        "vcf2_filename": vcf2_filename,
        "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "compatibility_warnings": compatibility_warnings or [],
        "stats1": stats1,
        "stats2": stats2,
        "summary": summary,
        "css_content": css_content,
        "figures": figures_html,
        "gene_table_html": _df_to_html_table(gene_comparison_df, max_rows=COMPARISON_GENE_DISPLAY_LIMIT),
        "gene_table_note": f"{len(gene_comparison_df)} genes affected in either strain.",
        "position_comparison_html": _df_to_html_table(display_df, max_rows=COMPARISON_POSITION_DISPLAY_LIMIT),
        "position_table_note": f"{len(position_comparison_df)} unique genomic positions across both strains.",
        "strain1_specific_html": _df_to_html_table(strain1_specific, max_rows=COMPARISON_SPECIFIC_DISPLAY_LIMIT),
        "strain2_specific_html": _df_to_html_table(strain2_specific, max_rows=COMPARISON_SPECIFIC_DISPLAY_LIMIT),
    }


def render_comparison_html_report(context: dict, output_path: str | Path) -> Path:
    """Render templates/comparison_report.html with the given context."""
    output_path = Path(output_path)
    ensure_output_dir(output_path.parent)

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)), autoescape=False)
    template = env.get_template("comparison_report.html")
    html_content = template.render(**context)

    output_path.write_text(html_content, encoding="utf-8")
    logger.info("Wrote comparison HTML report: %s", output_path)
    return output_path


def render_comparison_pdf_report(context: dict, output_path: str | Path) -> Path:
    """Render a PDF version of the comparison report via WeasyPrint."""
    try:
        from weasyprint import HTML
    except (ImportError, OSError) as exc:
        raise RuntimeError(
            "WeasyPrint could not load its native dependencies "
            "(Pango/cairo/gdk-pixbuf/libgobject) -- PDF generation skipped. "
            "See README.md for the macOS installation steps."
        ) from exc

    output_path = Path(output_path)
    ensure_output_dir(output_path.parent)

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)), autoescape=False)
    template = env.get_template("comparison_report.html")
    html_content = template.render(**context)

    try:
        HTML(string=html_content, base_url=str(TEMPLATES_DIR)).write_pdf(str(output_path))
    except Exception as exc:
        raise RuntimeError(f"PDF rendering failed: {exc}") from exc

    logger.info("Wrote comparison PDF report: %s", output_path)
    return output_path



def generate_comparison_reports(
    strain1_name: str,
    strain2_name: str,
    vcf1_filename: str,
    vcf2_filename: str,
    stats1: dict,
    stats2: dict,
    position_comparison_df: pd.DataFrame,
    gene_comparison_df: pd.DataFrame,
    summary: dict,
    figures: dict[str, "any"],
    output_dir: str | Path,
    make_pdf: bool = True,
    compatibility_warnings: list[str] | None = None,
) -> dict[str, Path]:
    """
    High-level convenience function for the comparison report, mirroring
    generate_reports() above but for the two-strain comparison layer.
    Writes <strain1>_<strain2>_comparison.html and .pdf.
    """
    out = ensure_output_dir(output_dir)
    results = {}
    prefix = f"{strain1_name}_{strain2_name}_comparison"

    interactive_figs = figures_to_interactive_html(figures)
    html_context = build_comparison_report_context(
        strain1_name, strain2_name, vcf1_filename, vcf2_filename, stats1, stats2,
        position_comparison_df, gene_comparison_df, summary, interactive_figs,
        compatibility_warnings,
    )
    results["html"] = render_comparison_html_report(html_context, out / f"{prefix}.html")

    if make_pdf:
        try:
            static_figs = figures_to_static_html(figures)
            pdf_context = build_comparison_report_context(
                strain1_name, strain2_name, vcf1_filename, vcf2_filename, stats1, stats2,
                position_comparison_df, gene_comparison_df, summary, static_figs,
                compatibility_warnings,
            )
            results["pdf"] = render_comparison_pdf_report(pdf_context, out / f"{prefix}.pdf")
        except RuntimeError as exc:
            logger.warning("Skipping comparison PDF generation: %s", exc)

    return results


# MULTI-STRAIN COMPARISON REPORT


MULTI_STRAIN_MATRIX_DISPLAY_LIMIT = 2000
MULTI_STRAIN_GENE_DISPLAY_LIMIT = 500


def _snp_matrix_for_display(matrix: pd.DataFrame, strain_names: list[str]) -> pd.DataFrame:
    """Prepare the SNP matrix for HTML rendering: adds a readable
    '{strain}_DNA_Change' column per strain, drops the bulky
    '{strain}_records' list-of-dict columns and raw REF/ALT allele-list
    columns (superseded by DNA_Change for display), and stringifies the
    Strains_Present list column. Full records/allele lists remain in the
    CSV export."""
    from comparison import add_dna_change_columns

    df = add_dna_change_columns(matrix, strain_names)
    drop_cols = [f"{name}_records" for name in strain_names]
    drop_cols += [f"{name}_REF_alleles" for name in strain_names]
    drop_cols += [f"{name}_ALT_alleles" for name in strain_names]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])
    # Also stringify the Strains_Present list column for cleaner display
    if "Strains_Present" in df.columns:
        df = df.copy()
        df["Strains_Present"] = df["Strains_Present"].apply(lambda v: ", ".join(v) if isinstance(v, list) else v)
    return df


def build_multi_strain_report_context(
    strain_names: list[str],
    vcf_filenames: dict[str, str],
    per_strain_stats: dict[str, dict],
    classified_matrix: pd.DataFrame,
    gene_comparison_df: pd.DataFrame,
    summary: dict,
    figures_html: dict[str, str],
    compatibility_warnings: list[str] | None = None,
) -> dict:
    """Assemble the Jinja2 context for the multi-strain report template."""
    css_content = (STATIC_DIR / "style.css").read_text(encoding="utf-8")
    display_matrix = _snp_matrix_for_display(classified_matrix, strain_names)

    return {
        "strain_names": strain_names,
        "vcf_filenames": vcf_filenames,
        "generated_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "compatibility_warnings": compatibility_warnings or [],
        "per_strain_stats": per_strain_stats,
        "summary": summary,
        "css_content": css_content,
        "figures": figures_html,
        "gene_table_html": _df_to_html_table(gene_comparison_df, max_rows=MULTI_STRAIN_GENE_DISPLAY_LIMIT),
        "gene_table_note": f"{len(gene_comparison_df)} genes affected in at least one of {len(strain_names)} strains.",
        "snp_matrix_html": _df_to_html_table(display_matrix, max_rows=MULTI_STRAIN_MATRIX_DISPLAY_LIMIT),
        "matrix_table_note": f"{len(classified_matrix)} unique genomic positions across {len(strain_names)} strains.",
    }


def render_multi_strain_html_report(context: dict, output_path: str | Path) -> Path:
    """Render templates/multi_strain_report.html with the given context."""
    output_path = Path(output_path)
    ensure_output_dir(output_path.parent)

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)), autoescape=False)
    template = env.get_template("multi_strain_report.html")
    html_content = template.render(**context)

    output_path.write_text(html_content, encoding="utf-8")
    logger.info("Wrote multi-strain HTML report: %s", output_path)
    return output_path


def render_multi_strain_pdf_report(context: dict, output_path: str | Path) -> Path:
    """Render a PDF version of the multi-strain report via WeasyPrint."""
    try:
        from weasyprint import HTML
    except (ImportError, OSError) as exc:
        raise RuntimeError(
            "WeasyPrint could not load its native dependencies "
            "(Pango/cairo/gdk-pixbuf/libgobject) -- PDF generation skipped. "
            "See README.md for the macOS installation steps."
        ) from exc

    output_path = Path(output_path)
    ensure_output_dir(output_path.parent)

    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)), autoescape=False)
    template = env.get_template("multi_strain_report.html")
    html_content = template.render(**context)

    try:
        HTML(string=html_content, base_url=str(TEMPLATES_DIR)).write_pdf(str(output_path))
    except Exception as exc:
        raise RuntimeError(f"PDF rendering failed: {exc}") from exc

    logger.info("Wrote multi-strain PDF report: %s", output_path)
    return output_path


def generate_multi_strain_reports(
    strain_names: list[str],
    vcf_filenames: dict[str, str],
    per_strain_stats: dict[str, dict],
    classified_matrix: pd.DataFrame,
    gene_comparison_df: pd.DataFrame,
    summary: dict,
    figures: dict[str, "any"],
    output_dir: str | Path,
    make_pdf: bool = True,
    compatibility_warnings: list[str] | None = None,
) -> dict[str, Path]:
    """
    High-level convenience function for the N-strain comparison report.
    Writes <N>strains_comparison.html and .pdf (or a joined-names prefix
    for small N, matching export.export_multi_strain_comparison's naming).
    """
    out = ensure_output_dir(output_dir)
    results = {}
    prefix = "_".join(strain_names) if len(strain_names) <= 4 else f"{len(strain_names)}strains"
    filename = f"{prefix}_comparison"

    interactive_figs = figures_to_interactive_html(figures)
    html_context = build_multi_strain_report_context(
        strain_names, vcf_filenames, per_strain_stats, classified_matrix,
        gene_comparison_df, summary, interactive_figs, compatibility_warnings,
    )
    results["html"] = render_multi_strain_html_report(html_context, out / f"{filename}.html")

    if make_pdf:
        try:
            static_figs = figures_to_static_html(figures)
            pdf_context = build_multi_strain_report_context(
                strain_names, vcf_filenames, per_strain_stats, classified_matrix,
                gene_comparison_df, summary, static_figs, compatibility_warnings,
            )
            results["pdf"] = render_multi_strain_pdf_report(pdf_context, out / f"{filename}.pdf")
        except RuntimeError as exc:
            logger.warning("Skipping multi-strain PDF generation: %s", exc)

    return results
