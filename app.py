"""
app.py
------
Streamlit web interface for the SNP Annotation Analyzer.

This file contains ONLY UI/display logic -- every analysis step (parsing,
annotation, classification, grouping, statistics, comparison, export,
reporting) lives in its own module and is imported here. This is what
lets main.py (CLI) reuse identical logic without importing Streamlit.

Run with:
    streamlit run app.py
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from comparison import (
    add_dna_change_columns,
    build_gene_comparison,
    build_gene_comparison_n_way,
    build_position_comparison,
    build_snp_matrix,
    classify_snp_matrix,
    classify_snp_sharing,
    comparison_summary,
    comparison_summary_n_way,
    conserved_vs_variable_positions,
    filter_comparison,
    filter_snp_matrix,
    get_isolates,
    get_position_detail,
    get_position_detail_n_way,
    isolate_summary_counts,
    pairwise_sharing_matrix_n_way,
    snp_sharing_matrix,
)
from export import (
    export_comparison_all,
    export_csv,
    export_excel,
    export_multi_strain_comparison,
    export_strain_labeled,
)
from grouping import group_by_gene, group_by_position
from pipeline import (
    apply_quality_filters,
    process_multiple_strains,
    process_multiple_vcfs,
    process_two_strains,
    process_vcf,
)
from report import generate_comparison_reports, generate_multi_strain_reports, generate_reports
from statistics import compute_filter_impact, compute_summary_statistics
from utils import get_logger
from validation import check_reference_compatibility
from visualisations import (
    build_comparison_figures,
    build_multi_strain_figures,
    fig_category_distribution,
    fig_effect_distribution,
    fig_isolate_sharing,
    fig_missense_vs_synonymous,
    fig_snp_distribution_along_chromosome,
    fig_snps_per_gene,
    fig_top_genes,
    fig_transition_transversion,
)
from vcf_parser import VCFValidationError

logger = get_logger(__name__)

st.set_page_config(page_title="SNP Annotation Analyser", layout="wide")

def _init_state():
    defaults = {
        "raw_df": None,          # unfiltered, fully processed SNP table
        "parse_reports": [],
        "vcf_names": [],
        # Strain Comparison page state (pairwise, 2 strains)
        "strain1_df": None,
        "strain2_df": None,
        "strain1_name": None,
        "strain2_name": None,
        "strain1_parse_report": None,
        "strain2_parse_report": None,
        "strain1_vcf_name": None,
        "strain2_vcf_name": None,
        "compatibility_warnings": [],
        # Multi-Strain Comparison page state (N strains)
        "multi_strains": None,       # dict: name -> {"df":..., "parse_report":...}
        "multi_vcf_filenames": None,  # dict: name -> original filename
        "multi_compat_warnings": [],
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


def _save_uploaded_files(uploaded_files) -> list[Path]:
    """Write Streamlit's in-memory UploadedFile objects to a temp dir so
    cyvcf2 (which needs a real file path) can read them. Files never
    leave this machine."""
    tmp_dir = Path(tempfile.mkdtemp(prefix="snp_analyser_"))
    paths = []
    for uf in uploaded_files:
        dest = tmp_dir / uf.name
        dest.write_bytes(uf.getbuffer())
        paths.append(dest)
    return paths


# ---------------------------------------------------------------------------
# Section renderers
# ---------------------------------------------------------------------------
def render_upload_section():
    st.title("🧬 SNP Annotation Analyser")
    st.caption(
        "Local, annotated-VCF SNP analysis. All processing happens on your "
        "machine -- no data is uploaded to an external server or API."
    )

    st.header("1. Upload VCF")
    uploaded_files = st.file_uploader(
        "Upload annotated VCF (.vcf or .vcf.gz). Upload multiple files to "
        "compare isolates.",
        type=["vcf", "gz"],
        accept_multiple_files=True,
    )

    if not uploaded_files:
        st.info("Upload one or more VCF files to begin.")
        return

    if st.button("Parse & Analyse", type="primary"):
        with st.spinner("Parsing and analysing VCF file(s)..."):
            try:
                paths = _save_uploaded_files(uploaded_files)
                if len(paths) == 1:
                    df, parse_report = process_vcf(paths[0])
                    parse_reports = [parse_report]
                else:
                    df, parse_reports = process_multiple_vcfs(paths)

                if df.empty:
                    st.error("No SNP records could be extracted from the uploaded file(s).")
                    return

                st.session_state["raw_df"] = df
                st.session_state["parse_reports"] = parse_reports
                st.session_state["vcf_names"] = [p.name for p in paths]
                st.success(f"Parsed {len(df)} SNP rows from {len(paths)} file(s).")
            except VCFValidationError as exc:
                st.error(f"VCF validation failed: {exc}")
            except Exception as exc:
                logger.exception("Unexpected error during parsing")
                st.error(f"Unexpected error while parsing: {exc}")

    if st.session_state["raw_df"] is not None:
        df = st.session_state["raw_df"]
        reports = st.session_state["parse_reports"]
        with st.expander("File details", expanded=True):
            cols = st.columns(4)
            cols[0].metric("Files", len(st.session_state["vcf_names"]))
            cols[1].metric("Samples", df["Sample"].nunique() if "Sample" in df else 0)
            cols[2].metric("Variants (rows)", len(df))
            cols[3].metric("Chromosome(s)", df["CHROM"].nunique() if "CHROM" in df else 0)

            st.write("**Filenames:**", ", ".join(st.session_state["vcf_names"]))
            st.write("**Samples:**", ", ".join(sorted(df["Sample"].dropna().unique().tolist())))
            st.write("**Chromosome/contig:**", ", ".join(sorted(df["CHROM"].dropna().unique().tolist())))

            st.write("**Parsing integrity:**")
            report_df = pd.DataFrame([r.as_dict() for r in reports])
            st.dataframe(report_df, use_container_width=True)


def render_filters_sidebar(df: pd.DataFrame) -> pd.DataFrame:
    st.sidebar.header("Filters")

    apply_qc = st.sidebar.checkbox("Apply quality filters", value=False)
    min_qual, min_dp = None, None
    if apply_qc:
        min_qual = st.sidebar.number_input("Minimum QUAL", min_value=0.0, value=0.0, step=1.0)
        min_dp = st.sidebar.number_input("Minimum DP", min_value=0, value=0, step=1)

    genes = sorted(df["Gene_Name"].dropna().unique().tolist())
    gene_filter = st.sidebar.multiselect("Gene", genes)

    categories = sorted(df["Mutation_Category"].dropna().unique().tolist())
    category_filter = st.sidebar.multiselect("Mutation category", categories)

    effects = sorted(df["Effect_Raw"].dropna().unique().tolist())
    effect_filter = st.sidebar.multiselect("Effect", effects)

    snp_types = sorted(df["SNP_Type"].dropna().unique().tolist())
    snp_type_filter = st.sidebar.multiselect("SNP type", snp_types)

    samples = sorted(df["Sample"].dropna().unique().tolist()) if "Sample" in df else []
    sample_filter = st.sidebar.multiselect("Sample", samples)

    search_text = st.sidebar.text_input("Text search (gene, position, effect...)")

    filtered = df.copy()
    records_before = len(filtered)

    if apply_qc:
        filtered = apply_quality_filters(filtered, min_qual=min_qual, min_dp=min_dp)
    if gene_filter:
        filtered = filtered[filtered["Gene_Name"].isin(gene_filter)]
    if category_filter:
        filtered = filtered[filtered["Mutation_Category"].isin(category_filter)]
    if effect_filter:
        filtered = filtered[filtered["Effect_Raw"].isin(effect_filter)]
    if snp_type_filter:
        filtered = filtered[filtered["SNP_Type"].isin(snp_type_filter)]
    if sample_filter:
        filtered = filtered[filtered["Sample"].isin(sample_filter)]
    if search_text:
        mask = filtered.astype(str).apply(
            lambda col: col.str.contains(search_text, case=False, na=False)
        ).any(axis=1)
        filtered = filtered[mask]

    impact = compute_filter_impact(records_before, len(filtered))
    st.sidebar.markdown(
        f"**{impact['records_after']} / {impact['records_before']}** records shown "
        f"({impact['pct_retained']}%)"
    )

    return filtered


def render_summary_cards(df: pd.DataFrame):
    st.header("2. Summary Statistics")
    stats = compute_summary_statistics(df)

    cols = st.columns(7)
    cols[0].metric("Total SNPs", stats["total_snps"])
    cols[1].metric("Missense", stats["categories"].get("Missense", {}).get("count", 0))
    cols[2].metric("Synonymous", stats["categories"].get("Silent / Synonymous", {}).get("count", 0))
    cols[3].metric("Nonsense", stats["categories"].get("Nonsense / Stop gained", {}).get("count", 0))
    cols[4].metric("Transitions", stats["transitions"])
    cols[5].metric("Transversions", stats["transversions"])
    cols[6].metric("Affected genes", stats["affected_genes"])

    st.caption(f"Ti/Tv ratio: **{stats['ti_tv_ratio']}** &nbsp;|&nbsp; "
               f"Annotated: **{stats['annotated_snps']}** ({stats['annotated_pct']}%) &nbsp;|&nbsp; "
               f"Unannotated: **{stats['unannotated_snps']}** ({stats['unannotated_pct']}%)")

    return stats


def render_snp_table(df: pd.DataFrame):
    st.header("3. Interactive SNP Table")
    display_cols = [
        "Sample", "CHROM", "POS", "REF", "ALT", "DNA_Change", "QUAL", "FILTER",
        "INFO_DP", "GT", "AD", "Effect_Raw", "Impact", "Mutation_Category",
        "Gene_Name", "Gene_ID", "HGVS_c", "HGVS_p", "SNP_Type",
        "has_multiple_annotations",
    ]
    available_cols = [c for c in display_cols if c in df.columns]
    st.dataframe(df[available_cols], use_container_width=True, height=400)

    col1, col2 = st.columns(2)
    with col1:
        csv_bytes = df.to_csv(index=False).encode("utf-8")
        st.download_button("⬇️ Download SNP table (CSV)", csv_bytes, "snp_table.csv", "text/csv")
    with col2:
        import io
        buf = io.BytesIO()
        export_df = df.copy()
        for c in export_df.columns:
            if export_df[c].apply(lambda v: isinstance(v, (list, dict))).any():
                export_df[c] = export_df[c].astype(str)
        export_df.to_excel(buf, index=False, engine="openpyxl")
        st.download_button(
            "⬇️ Download SNP table (Excel)", buf.getvalue(), "snp_table.xlsx",
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )


def render_gene_explorer(df: pd.DataFrame):
    st.header("4. Gene Explorer")
    gene_df = group_by_gene(df)

    display_cols = [
        "Gene_Name", "Gene_ID", "Total_SNPs", "Unique_Positions", "Missense",
        "Synonymous", "Nonsense", "Frameshift", "Other_Effects",
        "Transitions", "Transversions", "Ti_Tv_Ratio",
    ]
    st.dataframe(gene_df[display_cols], use_container_width=True, height=350)

    st.download_button(
        "⬇️ Download gene summary (CSV)",
        gene_df.drop(columns=["Positions", "DNA_Changes", "Protein_Changes", "Effects"]).to_csv(index=False).encode("utf-8"),
        "gene_summary.csv", "text/csv",
    )

    st.subheader("Inspect a gene")
    gene_names = gene_df["Gene_Name"].tolist()
    selected_gene = st.selectbox("Select a gene", gene_names) if gene_names else None
    if selected_gene:
        gene_snps = df[df["Gene_Name"] == selected_gene].sort_values("POS")
        st.write(f"**{selected_gene}** -- {len(gene_snps)} SNP(s)")
        st.dataframe(
            gene_snps[["POS", "DNA_Change", "Effect_Raw", "Mutation_Category", "HGVS_p", "SNP_Type"]],
            use_container_width=True,
        )
    return gene_df


def render_position_explorer(df: pd.DataFrame):
    st.header("5. Position Explorer")
    position_df = group_by_position(df)

    st.write("Search a genomic position:")
    search_pos = st.number_input("Position", min_value=0, value=0, step=1)
    if search_pos:
        match = position_df[position_df["POS"] == search_pos]
        if match.empty:
            st.warning(f"No SNP found at position {search_pos}.")
        else:
            row = match.iloc[0]
            st.write(f"**Chromosome:** {row['CHROM']}  |  **Position:** {row['POS']}  |  **Reference:** {row['REF']}")
            st.write(f"**Gene:** {row['Gene_Name']}  |  **Effect:** {row['Effect_Raw']}  |  **Category:** {row['Mutation_Category']}")
            st.write("**Alleles observed:**")
            for allele in row["Alt_Alleles"]:
                count = row["Allele_Isolate_Counts"].get(allele, 0)
                pct = row["Allele_Isolate_Percentages"].get(allele, 0.0)
                st.write(f"- {row['REF']}>{allele}: {count} isolate(s) ({pct}%)")

    with st.expander("Full position-level summary table"):
        display_cols = ["CHROM", "POS", "REF", "Num_Alt_Alleles", "Alt_Alleles", "Gene_Name", "Effect_Raw", "Mutation_Category"]
        st.dataframe(position_df[display_cols], use_container_width=True, height=350)
        st.download_button(
            "⬇️ Download position summary (CSV)",
            position_df.drop(columns=["Alt_Alleles", "Allele_Isolate_Counts", "Allele_Isolate_Percentages"]).to_csv(index=False).encode("utf-8"),
            "position_summary.csv", "text/csv",
        )

    return position_df


def render_multi_isolate_section(df: pd.DataFrame):
    isolates = get_isolates(df)
    if len(isolates) < 2:
        return None

    st.header("6. Multi-Isolate Comparison")
    st.write(f"Comparing **{len(isolates)}** isolates: {', '.join(isolates)}")

    iso_summary = isolate_summary_counts(df)
    st.dataframe(iso_summary, use_container_width=True)

    sharing = classify_snp_sharing(df)
    shared_count = int((sharing["sharing_status"] == "Shared").sum())
    specific_count = int((sharing["sharing_status"] == "Isolate-specific").sum())
    col1, col2 = st.columns(2)
    col1.metric("Shared SNPs", shared_count)
    col2.metric("Isolate-specific SNPs", specific_count)

    conserved = conserved_vs_variable_positions(df)
    conserved_count = int((conserved["Status"] == "Conserved").sum())
    variable_count = int((conserved["Status"] == "Variable").sum())
    col3, col4 = st.columns(2)
    col3.metric("Conserved positions", conserved_count)
    col4.metric("Variable positions", variable_count)

    matrix = snp_sharing_matrix(df)
    fig = fig_isolate_sharing(matrix)
    if fig:
        st.plotly_chart(fig, use_container_width=True)

    with st.expander("Full SNP sharing table"):
        st.dataframe(sharing, use_container_width=True)

    return matrix


def render_visualizations(df: pd.DataFrame, gene_df: pd.DataFrame):
    st.header("7. Visualisations")

    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(fig_effect_distribution(df), use_container_width=True)
    with col2:
        st.plotly_chart(fig_category_distribution(df), use_container_width=True)

    col3, col4 = st.columns(2)
    with col3:
        st.plotly_chart(fig_transition_transversion(df), use_container_width=True)
    with col4:
        st.plotly_chart(fig_missense_vs_synonymous(df), use_container_width=True)

    st.plotly_chart(fig_top_genes(gene_df, top_n=15), use_container_width=True)
    st.plotly_chart(fig_snps_per_gene(gene_df, top_n=30), use_container_width=True)
    st.plotly_chart(fig_snp_distribution_along_chromosome(df), use_container_width=True)


def render_report_section(df, gene_df, position_df, stats, parse_report, vcf_filename, sharing_matrix):
    st.header("8. Full Report")
    st.write("Generate a standalone HTML report and a publication-ready PDF report.")

    make_pdf = st.checkbox("Also generate PDF (requires WeasyPrint + Chrome for chart images)", value=True)

    if st.button("Generate Reports"):
        with st.spinner("Building report... this may take a moment for large VCFs."):
            from visualisations import build_all_figures
            figures = build_all_figures(df, gene_df, sharing_matrix)
            try:
                results = generate_reports(
                    vcf_filename, df, gene_df, position_df, stats, parse_report,
                    figures, "output", make_pdf=make_pdf,
                )
                st.success("Report(s) generated.")
                if "html" in results:
                    st.download_button(
                        "⬇️ Download HTML report",
                        Path(results["html"]).read_bytes(),
                        "report.html", "text/html",
                    )
                if "pdf" in results:
                    st.download_button(
                        "⬇️ Download PDF report",
                        Path(results["pdf"]).read_bytes(),
                        "report.pdf", "application/pdf",
                    )
            except Exception as exc:
                logger.exception("Report generation failed")
                st.error(f"Report generation failed: {exc}")


def render_individual_analysis_page():
    """The original single-VCF analysis workflow, unchanged in behavior --
    now presented as one page within the multi-page navigation."""
    render_upload_section()

    if st.session_state["raw_df"] is None:
        return

    full_df = st.session_state["raw_df"]
    filtered_df = render_filters_sidebar(full_df)

    if filtered_df.empty:
        st.warning("No records match the current filters.")
        return

    stats = render_summary_cards(filtered_df)
    render_snp_table(filtered_df)
    gene_df = render_gene_explorer(filtered_df)
    position_df = render_position_explorer(filtered_df)
    sharing_matrix = render_multi_isolate_section(filtered_df)
    render_visualizations(filtered_df, gene_df)

    vcf_filename = ", ".join(st.session_state["vcf_names"])
    parse_report = st.session_state["parse_reports"][0]
    render_report_section(filtered_df, gene_df, position_df, stats, parse_report, vcf_filename, sharing_matrix)



# STRAIN COMPARISON PAGE

def render_comparison_upload_section():
    st.title("🧬🆚🧬 Strain Comparison")
    st.caption(
        "Compare two independently-analysed strains, linked by genomic coordinate "
        "(CHROM + POS). Neither strain's individual data is modified or reduced."
    )

    st.header("1. Upload Two VCFs")
    col1, col2 = st.columns(2)
    with col1:
        vcf1_file = st.file_uploader("Strain 1 VCF", type=["vcf", "gz"], key="strain1_uploader")
        name1_input = st.text_input(
            "Strain 1 name (optional)", value="",
            placeholder="e.g. your_strain_name1 -- leave blank to use the filename",
            key="strain1_name_input",
        )
    with col2:
        vcf2_file = st.file_uploader("Strain 2 VCF", type=["vcf", "gz"], key="strain2_uploader")
        name2_input = st.text_input(
            "Strain 2 name (optional)", value="",
            placeholder="e.g. your_strain_name2 -- leave blank to use the filename",
            key="strain2_name_input",
        )

    st.caption(
        "⚠️ Both VCFs' sample columns may contain the same internal sample name "
        "(e.g. both say 'strain1') -- provide distinct strain names above so the "
        "comparison can tell them apart. If left blank, the filename is used instead."
    )

    if not (vcf1_file and vcf2_file):
        st.info("Upload both VCF files to begin.")
        return

    if st.button("Analyse & Compare", type="primary"):
        with st.spinner("Parsing and analysing both strains independently..."):
            try:
                paths = _save_uploaded_files([vcf1_file, vcf2_file])
                strain1_name = name1_input.strip() or Path(vcf1_file.name).stem
                strain2_name = name2_input.strip() or Path(vcf2_file.name).stem

                if strain1_name == strain2_name:
                    st.error(
                        f"Strain 1 and Strain 2 have the same name ('{strain1_name}'). "
                        "Please provide distinct names so the comparison can distinguish them."
                    )
                    return

                df1, df2, report1, report2 = process_two_strains(paths[0], paths[1], strain1_name, strain2_name)

                if df1.empty or df2.empty:
                    empty_strain = strain1_name if df1.empty else strain2_name
                    st.error(f"No SNP records could be extracted from {empty_strain}'s VCF.")
                    return

                compat = check_reference_compatibility(report1, report2, strain1_name, strain2_name)

                st.session_state["strain1_df"] = df1
                st.session_state["strain2_df"] = df2
                st.session_state["strain1_name"] = strain1_name
                st.session_state["strain2_name"] = strain2_name
                st.session_state["strain1_parse_report"] = report1
                st.session_state["strain2_parse_report"] = report2
                st.session_state["strain1_vcf_name"] = vcf1_file.name
                st.session_state["strain2_vcf_name"] = vcf2_file.name
                st.session_state["compatibility_warnings"] = compat.warnings

                st.success(f"Analysed {strain1_name} ({len(df1)} SNPs) and {strain2_name} ({len(df2)} SNPs).")
            except VCFValidationError as exc:
                st.error(f"VCF validation failed: {exc}")
            except Exception as exc:
                logger.exception("Unexpected error during two-strain analysis")
                st.error(f"Unexpected error: {exc}")

    if st.session_state["compatibility_warnings"]:
        for w in st.session_state["compatibility_warnings"]:
            st.warning(f"⚠️ {w}")

    if st.session_state["strain1_df"] is not None:
        df1 = st.session_state["strain1_df"]
        df2 = st.session_state["strain2_df"]
        n1, n2 = st.session_state["strain1_name"], st.session_state["strain2_name"]
        with st.expander("Strain file details", expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                st.write(f"**{n1}**")
                st.write(f"Filename: {st.session_state['strain1_vcf_name']}")
                st.write(f"Sample name in VCF: {', '.join(st.session_state['strain1_parse_report'].samples)}")
                st.write(f"Total SNPs: {len(df1)}")
            with col2:
                st.write(f"**{n2}**")
                st.write(f"Filename: {st.session_state['strain2_vcf_name']}")
                st.write(f"Sample name in VCF: {', '.join(st.session_state['strain2_parse_report'].samples)}")
                st.write(f"Total SNPs: {len(df2)}")


def render_comparison_summary_cards(summary: dict, n1: str, n2: str):
    st.header("2. Comparison Summary")
    cols = st.columns(4)
    cols[0].metric("Total unique positions", summary["total_unique_positions"])
    cols[1].metric("Shared positions", summary["shared_positions"])
    cols[2].metric(f"{n1}-specific", summary[f"{n1}_specific"])
    cols[3].metric(f"{n2}-specific", summary[f"{n2}_specific"])

    cols2 = st.columns(4)
    cols2[0].metric("Identical mutations", summary["identical_mutations"])
    cols2[1].metric("Same position, diff. mutation", summary["same_position_different_mutation"])
    cols2[2].metric("Reference discrepancies", summary["reference_discrepancies"])
    cols2[3].metric("Shared genes", summary["shared_genes"])


def render_comparison_filters_sidebar(position_comparison_df: pd.DataFrame, n1: str, n2: str) -> pd.DataFrame:
    st.sidebar.header("Comparison Filters")

    all_statuses = sorted(position_comparison_df["Comparison_Status"].unique().tolist())
    status_filter = st.sidebar.multiselect("Comparison status", all_statuses)

    genes = sorted(position_comparison_df["Gene_Name"].dropna().unique().tolist())
    gene_filter = st.sidebar.multiselect("Gene", genes)

    search_text = st.sidebar.text_input("Search (gene, position, protein/DNA change...)")

    filtered = filter_comparison(position_comparison_df, status_filter or None, gene_filter or None, search_text or None)
    st.sidebar.markdown(f"**{len(filtered)} / {len(position_comparison_df)}** positions shown")
    return filtered


def render_comparison_table(position_comparison_df: pd.DataFrame, n1: str, n2: str):
    st.header("3. Position-Level Comparison Table")
    display = add_dna_change_columns(position_comparison_df, [n1, n2])
    display_cols = [
        "CHROM", "POS", "Gene_Name", "Comparison_Status",
        f"{n1}_DNA_Change", f"{n2}_DNA_Change",
    ]
    available = [c for c in display_cols if c in display.columns]
    st.dataframe(display[available], use_container_width=True, height=400)


def render_position_detail_view(position_comparison_df: pd.DataFrame, n1: str, n2: str):
    st.header("4. Position Detail View")
    search_pos = st.number_input("Enter a genomic position to inspect", min_value=0, value=0, step=1)
    if not search_pos:
        return

    chroms = position_comparison_df["CHROM"].unique().tolist()
    chrom = chroms[0] if len(chroms) == 1 else st.selectbox("Chromosome/contig", chroms)

    detail = get_position_detail(position_comparison_df, chrom, search_pos)
    if detail is None:
        st.warning(f"No SNP found at position {search_pos} in either strain.")
        return

    st.subheader(f"Position {search_pos} -- {detail['Comparison_Status']}")
    col1, col2 = st.columns(2)

    for col, strain_name in [(col1, n1), (col2, n2)]:
        with col:
            st.write(f"**{strain_name}**")
            records = detail.get(f"{strain_name}_records", [])
            if not records:
                st.write("_No variant at this position._")
                continue
            for rec in records:
                st.write(f"REF: {rec.get('REF')}  |  ALT: {rec.get('ALT')}  |  DNA change: {rec.get('DNA_Change')}")
                st.write(f"QUAL: {rec.get('QUAL')}  |  DP: {rec.get('INFO_DP')}  |  GT: {rec.get('GT')}")
                st.write(f"Gene: {rec.get('Gene_Name')}  |  Effect: {rec.get('Effect_Raw')}")
                st.write(f"Category: {rec.get('Mutation_Category')}  |  Protein: {rec.get('HGVS_p')}")
                st.markdown("---")


def render_gene_comparison_table(gene_comparison_df: pd.DataFrame, n1: str, n2: str):
    st.header("5. Gene-Level Comparison")
    display_cols = [
        "Gene_Name", f"{n1}_SNPs", f"{n2}_SNPs", "Shared_Positions",
        f"{n1}_Specific_Positions", f"{n2}_Specific_Positions",
        f"{n1}_Missense", f"{n2}_Missense", "Shared_Gene",
    ]
    available = [c for c in display_cols if c in gene_comparison_df.columns]
    st.dataframe(gene_comparison_df[available], use_container_width=True, height=350)


def render_comparison_visualizations(df1, df2, gene_comparison_df, summary, n1, n2):
    st.header("6. Comparison Visualisations")
    figs = build_comparison_figures(df1, df2, gene_comparison_df, summary, n1, n2)

    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(figs["snp_counts_by_strain"], use_container_width=True)
    with col2:
        st.plotly_chart(figs["shared_vs_specific"], use_container_width=True)

    st.plotly_chart(figs["category_by_strain"], use_container_width=True)

    col3, col4 = st.columns(2)
    with col3:
        st.plotly_chart(figs["missense_synonymous_comparison"], use_container_width=True)
    with col4:
        st.plotly_chart(figs["titv_comparison"], use_container_width=True)

    col5, col6 = st.columns(2)
    with col5:
        st.plotly_chart(figs["shared_genes"], use_container_width=True)
    with col6:
        st.plotly_chart(figs["top_genes_comparison"], use_container_width=True)

    st.plotly_chart(figs["chromosome_distribution_by_strain"], use_container_width=True)
    return figs


def render_comparison_report_section(df1, df2, gene_comparison_df, position_comparison_df, summary, stats1, stats2, n1, n2, figures):
    st.header("7. Comparison Report")
    make_pdf = st.checkbox("Also generate PDF (requires WeasyPrint + Chrome)", value=True, key="comparison_pdf_checkbox")

    if st.button("Generate Comparison Reports"):
        with st.spinner("Building comparison report..."):
            try:
                results = generate_comparison_reports(
                    n1, n2, st.session_state["strain1_vcf_name"], st.session_state["strain2_vcf_name"],
                    stats1, stats2, position_comparison_df, gene_comparison_df, summary, figures,
                    "output/comparison", make_pdf=make_pdf,
                    compatibility_warnings=st.session_state["compatibility_warnings"],
                )
                st.success("Comparison report(s) generated.")
                if "html" in results:
                    st.download_button(
                        "⬇️ Download comparison HTML report",
                        Path(results["html"]).read_bytes(),
                        f"{n1}_{n2}_comparison.html", "text/html",
                    )
                if "pdf" in results:
                    st.download_button(
                        "⬇️ Download comparison PDF report",
                        Path(results["pdf"]).read_bytes(),
                        f"{n1}_{n2}_comparison.pdf", "application/pdf",
                    )
            except Exception as exc:
                logger.exception("Comparison report generation failed")
                st.error(f"Report generation failed: {exc}")

    st.subheader("Export comparison tables")
    col1, col2, col3 = st.columns(3)
    with col1:
        shared = position_comparison_df[position_comparison_df["Comparison_Status"].isin(
            ["Identical mutation", "Same position, different mutation", "Reference allele discrepancy"]
        )]
        st.download_button(
            "⬇️ Shared SNPs (CSV)",
            shared.drop(columns=[c for c in shared.columns if c.endswith("_records")]).to_csv(index=False).encode("utf-8"),
            f"{n1}_{n2}_shared_SNPs.csv", "text/csv",
        )
    with col2:
        s1_specific = position_comparison_df[position_comparison_df["Comparison_Status"] == f"{n1}-specific"]
        st.download_button(
            f"⬇️ {n1}-specific (CSV)",
            s1_specific.drop(columns=[c for c in s1_specific.columns if c.endswith("_records")]).to_csv(index=False).encode("utf-8"),
            f"{n1}_{n2}_{n1}_specific.csv", "text/csv",
        )
    with col3:
        s2_specific = position_comparison_df[position_comparison_df["Comparison_Status"] == f"{n2}-specific"]
        st.download_button(
            f"⬇️ {n2}-specific (CSV)",
            s2_specific.drop(columns=[c for c in s2_specific.columns if c.endswith("_records")]).to_csv(index=False).encode("utf-8"),
            f"{n1}_{n2}_{n2}_specific.csv", "text/csv",
        )

    st.download_button(
        "⬇️ Gene comparison (CSV)",
        gene_comparison_df.to_csv(index=False).encode("utf-8"),
        f"{n1}_{n2}_gene_comparison.csv", "text/csv",
    )


def render_strain_comparison_page():
    render_comparison_upload_section()

    if st.session_state["strain1_df"] is None or st.session_state["strain2_df"] is None:
        return

    df1 = st.session_state["strain1_df"]
    df2 = st.session_state["strain2_df"]
    n1 = st.session_state["strain1_name"]
    n2 = st.session_state["strain2_name"]

    with st.spinner("Building comparison..."):
        position_comparison_df = build_position_comparison(df1, df2, n1, n2)
        gene_comparison_df = build_gene_comparison(df1, df2, n1, n2)
        summary = comparison_summary(position_comparison_df, gene_comparison_df, n1, n2)
        stats1 = compute_summary_statistics(df1)
        stats2 = compute_summary_statistics(df2)

    render_comparison_summary_cards(summary, n1, n2)
    filtered_comparison_df = render_comparison_filters_sidebar(position_comparison_df, n1, n2)
    render_comparison_table(filtered_comparison_df, n1, n2)
    render_position_detail_view(position_comparison_df, n1, n2)
    render_gene_comparison_table(gene_comparison_df, n1, n2)
    figures = render_comparison_visualizations(df1, df2, gene_comparison_df, summary, n1, n2)
    render_comparison_report_section(df1, df2, gene_comparison_df, position_comparison_df, summary, stats1, stats2, n1, n2, figures)


# ===========================================================================
# MULTI-STRAIN COMPARISON PAGE (N strains, 2 to 10+)
# ===========================================================================

def render_multi_strain_upload_section():
    st.title("🧬🔗 Multi-Strain Comparison")
    st.caption(
        "Compare any number of strains (2, 10, 50...) at once, linked by genomic "
        "coordinate. Each strain is analysed completely independently first; "
        "no strain's data is modified or reduced by the comparison."
    )

    st.header("1. Upload VCFs")
    uploaded_files = st.file_uploader(
        "Upload 2 or more annotated VCFs (.vcf or .vcf.gz) -- one per strain",
        type=["vcf", "gz"], accept_multiple_files=True, key="multi_strain_uploader",
    )

    if not uploaded_files:
        st.info("Upload 2 or more VCF files to begin.")
        return

    if len(uploaded_files) < 2:
        st.warning("Upload at least 2 files to run a comparison.")
        return

    st.caption(
        "⚠️ Different strains' VCFs often share the same internal sample name "
        "(e.g. all say 'isolate1') -- give each file below a distinct strain "
        "name so the comparison can tell them apart. Defaults are derived from "
        "the filenames."
    )

    st.subheader(f"Name each of the {len(uploaded_files)} strains")
    strain_name_inputs = []
    cols = st.columns(2)
    for i, uf in enumerate(uploaded_files):
        default_name = Path(uf.name).stem
        with cols[i % 2]:
            name = st.text_input(
                f"Strain name for `{uf.name}`", value=default_name, key=f"multi_strain_name_{i}",
            )
            strain_name_inputs.append(name.strip())

    if st.button("Analyse & Compare All Strains", type="primary"):
        if len(set(strain_name_inputs)) != len(strain_name_inputs):
            counts = {n: strain_name_inputs.count(n) for n in set(strain_name_inputs)}
            dupes = [n for n, c in counts.items() if c > 1]
            st.error(f"Strain names must be unique. Duplicate name(s): {', '.join(dupes)}")
            return
        if any(not n for n in strain_name_inputs):
            st.error("Every strain needs a non-empty name.")
            return

        with st.spinner(f"Parsing and analysing {len(uploaded_files)} strains independently..."):
            try:
                paths = _save_uploaded_files(uploaded_files)
                strains = process_multiple_strains(paths, strain_name_inputs)

                empty = [name for name, data in strains.items() if data["df"].empty]
                if empty:
                    st.error(f"No SNP records could be extracted for: {', '.join(empty)}")
                    return

                # Pairwise reference compatibility checks across all strains
                warnings = []
                names = list(strains.keys())
                for i, name_a in enumerate(names):
                    for name_b in names[i + 1:]:
                        compat = check_reference_compatibility(
                            strains[name_a]["parse_report"], strains[name_b]["parse_report"],
                            name_a, name_b,
                        )
                        warnings.extend(compat.warnings)

                st.session_state["multi_strains"] = strains
                st.session_state["multi_vcf_filenames"] = {
                    name: uf.name for name, uf in zip(strain_name_inputs, uploaded_files)
                }
                st.session_state["multi_compat_warnings"] = warnings

                total_snps = sum(len(data["df"]) for data in strains.values())
                st.success(f"Analysed {len(strains)} strains, {total_snps} total SNPs.")
            except VCFValidationError as exc:
                st.error(f"VCF validation failed: {exc}")
            except ValueError as exc:
                st.error(str(exc))
            except Exception as exc:
                logger.exception("Unexpected error during multi-strain analysis")
                st.error(f"Unexpected error: {exc}")

    if st.session_state["multi_compat_warnings"]:
        for w in st.session_state["multi_compat_warnings"]:
            st.warning(f"⚠️ {w}")

    if st.session_state["multi_strains"] is not None:
        strains = st.session_state["multi_strains"]
        with st.expander("Strain file details", expanded=True):
            rows = []
            for name, data in strains.items():
                rows.append({
                    "Strain": name,
                    "Filename": st.session_state["multi_vcf_filenames"].get(name, ""),
                    "VCF sample name": ", ".join(data["parse_report"].samples),
                    "Total SNPs": len(data["df"]),
                })
            st.dataframe(pd.DataFrame(rows), use_container_width=True)


def render_multi_strain_summary_cards(summary: dict):
    st.header("2. Comparison Summary")
    cols = st.columns(4)
    cols[0].metric("Total unique positions", summary["total_unique_positions"])
    cols[1].metric(f"Core positions (all {summary['num_strains']})", summary["core_positions_total"])
    cols[2].metric("Partial sharing", summary["partial_positions"])
    cols[3].metric("Reference discrepancies", summary["reference_discrepancies"])

    cols2 = st.columns(4)
    cols2[0].metric("Core, identical", summary["core_identical"])
    cols2[1].metric("Core, variable", summary["core_variable"])
    cols2[2].metric("Core genes", summary["core_genes"])
    cols2[3].metric("Total genes (any strain)", summary["total_genes"])

    st.subheader("Strain-specific position counts")
    unique_df = pd.DataFrame(
        list(summary["strain_unique_counts"].items()), columns=["Strain", "Unique_Positions"]
    )
    st.dataframe(unique_df, use_container_width=True)


def render_multi_strain_filters_sidebar(classified_matrix: pd.DataFrame, strain_names: list[str]) -> pd.DataFrame:
    st.sidebar.header("Comparison Filters")

    all_statuses = sorted(classified_matrix["Comparison_Status"].unique().tolist())
    status_filter = st.sidebar.multiselect("Comparison status", all_statuses)

    genes = sorted(classified_matrix["Gene_Name"].dropna().unique().tolist())
    gene_filter = st.sidebar.multiselect("Gene", genes)

    strain_filter = st.sidebar.multiselect("Present in strain(s)", strain_names)

    search_text = st.sidebar.text_input("Search (gene, position...)")

    filtered = filter_snp_matrix(
        classified_matrix, status_filter or None, gene_filter or None,
        strain_filter or None, search_text or None,
    )
    st.sidebar.markdown(f"**{len(filtered)} / {len(classified_matrix)}** positions shown")
    return filtered


def render_multi_strain_matrix_table(classified_matrix: pd.DataFrame, strain_names: list[str]):
    st.header("3. SNP Matrix (Position-Level Comparison)")
    display = add_dna_change_columns(classified_matrix, strain_names)
    display = display.copy()
    display["Strains_Present"] = display["Strains_Present"].apply(
        lambda v: ", ".join(v) if isinstance(v, list) else v
    )
    display_cols = ["CHROM", "POS", "Gene_Name", "Comparison_Status", "Num_Strains_Present", "Strains_Present"]
    display_cols += [f"{name}_DNA_Change" for name in strain_names]
    available = [c for c in display_cols if c in display.columns]
    st.dataframe(display[available], use_container_width=True, height=400)


def render_multi_strain_position_detail(classified_matrix: pd.DataFrame, strain_names: list[str]):
    st.header("4. Position Detail View")
    search_pos = st.number_input("Enter a genomic position to inspect", min_value=0, value=0, step=1, key="multi_pos_search")
    if not search_pos:
        return

    chroms = classified_matrix["CHROM"].unique().tolist()
    chrom = chroms[0] if len(chroms) == 1 else st.selectbox("Chromosome/contig", chroms, key="multi_chrom_select")

    detail = get_position_detail_n_way(classified_matrix, chrom, search_pos, strain_names)
    if detail is None:
        st.warning(f"No SNP found at position {search_pos} in any strain.")
        return

    st.subheader(f"Position {search_pos} -- {detail['Comparison_Status']}")
    n_cols = min(len(strain_names), 4)
    cols = st.columns(n_cols)
    for i, strain_name in enumerate(strain_names):
        with cols[i % n_cols]:
            st.write(f"**{strain_name}**")
            records = detail.get(f"{strain_name}_records", [])
            if not records:
                st.write("_No variant._")
                continue
            for rec in records:
                st.write(f"{rec.get('REF')}>{rec.get('ALT')}")
                st.write(f"Gene: {rec.get('Gene_Name')}")
                st.write(f"Effect: {rec.get('Effect_Raw')}")
                st.write(f"DP: {rec.get('INFO_DP')} | QUAL: {rec.get('QUAL')}")
                st.markdown("---")


def render_multi_strain_gene_table(gene_comparison_df: pd.DataFrame, strain_names: list[str]):
    st.header("5. Gene-Level Comparison")
    snp_cols = [f"{name}_SNPs" for name in strain_names]
    display_cols = ["Gene_Name"] + snp_cols + ["Num_Strains_Affected", "Core_Gene"]
    available = [c for c in display_cols if c in gene_comparison_df.columns]
    st.dataframe(gene_comparison_df[available], use_container_width=True, height=350)


def render_multi_strain_visualizations(strains: dict, classified_matrix, gene_comparison_df, summary, pairwise_matrix):
    st.header("6. Comparison Visualisations")
    figs = build_multi_strain_figures(strains, classified_matrix, gene_comparison_df, summary, pairwise_matrix)

    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(figs["snps_per_strain"], use_container_width=True)
    with col2:
        st.plotly_chart(figs["core_genome_summary"], use_container_width=True)

    st.plotly_chart(figs["pairwise_sharing"], use_container_width=True)

    col3, col4 = st.columns(2)
    with col3:
        st.plotly_chart(figs["strain_unique_counts"], use_container_width=True)
    with col4:
        st.plotly_chart(figs["genes_by_strain_count"], use_container_width=True)

    st.caption("Presence/absence matrix is capped to the top 300 positions by strain count for renderability.")
    st.plotly_chart(figs["presence_absence"], use_container_width=True)
    return figs


def render_multi_strain_report_section(strains, classified_matrix, gene_comparison_df, summary, pairwise_matrix, strain_names, figures):
    st.header("7. Comparison Report")
    make_pdf = st.checkbox("Also generate PDF (requires WeasyPrint + Chrome)", value=True, key="multi_strain_pdf_checkbox")

    if st.button("Generate Multi-Strain Reports"):
        with st.spinner("Building multi-strain report..."):
            try:
                per_strain_stats = {name: compute_summary_statistics(strains[name]["df"]) for name in strain_names}
                results = generate_multi_strain_reports(
                    strain_names, st.session_state["multi_vcf_filenames"], per_strain_stats,
                    classified_matrix, gene_comparison_df, summary, figures,
                    "output/multi_strain", make_pdf=make_pdf,
                    compatibility_warnings=st.session_state["multi_compat_warnings"],
                )
                st.success("Multi-strain report(s) generated.")
                prefix = "_".join(strain_names) if len(strain_names) <= 4 else f"{len(strain_names)}strains"
                if "html" in results:
                    st.download_button(
                        "⬇️ Download comparison HTML report",
                        Path(results["html"]).read_bytes(),
                        f"{prefix}_comparison.html", "text/html",
                    )
                if "pdf" in results:
                    st.download_button(
                        "⬇️ Download comparison PDF report",
                        Path(results["pdf"]).read_bytes(),
                        f"{prefix}_comparison.pdf", "application/pdf",
                    )
            except Exception as exc:
                logger.exception("Multi-strain report generation failed")
                st.error(f"Report generation failed: {exc}")

    st.subheader("Export comparison tables")
    col1, col2, col3 = st.columns(3)
    strain_records_cols = [c for c in classified_matrix.columns if c.endswith("_records")]
    export_matrix = classified_matrix.drop(columns=strain_records_cols)
    with col1:
        st.download_button(
            "⬇️ Full SNP matrix (CSV)",
            export_matrix.to_csv(index=False).encode("utf-8"),
            "snp_matrix.csv", "text/csv",
        )
    with col2:
        core_only = export_matrix[export_matrix["Comparison_Status"].isin(["Core, identical", "Core, variable"])]
        st.download_button(
            "⬇️ Core genome positions (CSV)",
            core_only.to_csv(index=False).encode("utf-8"),
            "core_genome_positions.csv", "text/csv",
        )
    with col3:
        st.download_button(
            "⬇️ Gene comparison (CSV)",
            gene_comparison_df.to_csv(index=False).encode("utf-8"),
            "gene_comparison.csv", "text/csv",
        )


def render_multi_strain_comparison_page():
    render_multi_strain_upload_section()

    if st.session_state["multi_strains"] is None:
        return

    strains = st.session_state["multi_strains"]
    strain_names = list(strains.keys())

    with st.spinner("Building N-way comparison..."):
        matrix = build_snp_matrix(strains)
        classified_matrix = classify_snp_matrix(matrix, strain_names)
        gene_comparison_df = build_gene_comparison_n_way(strains)
        summary = comparison_summary_n_way(classified_matrix, gene_comparison_df, strain_names)
        pairwise_matrix = pairwise_sharing_matrix_n_way(strains)

    render_multi_strain_summary_cards(summary)
    filtered_matrix = render_multi_strain_filters_sidebar(classified_matrix, strain_names)
    render_multi_strain_matrix_table(filtered_matrix, strain_names)
    render_multi_strain_position_detail(classified_matrix, strain_names)
    render_multi_strain_gene_table(gene_comparison_df, strain_names)
    figures = render_multi_strain_visualizations(strains, classified_matrix, gene_comparison_df, summary, pairwise_matrix)
    render_multi_strain_report_section(strains, classified_matrix, gene_comparison_df, summary, pairwise_matrix, strain_names, figures)


def render_home_page():
    st.title("🧬 SNP Annotation Analyser")
    st.caption("Local, annotated-VCF SNP analysis for bacterial genome resequencing.")
    st.markdown(
        """
        Use the sidebar to navigate:

        - **Individual Analysis** -- upload one annotated VCF (one strain, or
          multiple isolates of the same strain) for complete SNP classification,
          gene/position grouping, statistics, and reports.
        - **Strain Comparison** -- upload two VCFs from two independent strains
          for a focused, deep-dive comparison by genomic coordinate.
        - **Multi-Strain Comparison** -- upload any number of strains (2, 10,
          50+) at once and compare them all together: core genome positions
          shared by every strain, partial sharing, strain-specific variants,
          and gene-level comparison across the whole set.

        All processing happens locally on this machine -- no VCF data is
        uploaded to an external server or API.
        """
    )


def main():
    _init_state()

    st.sidebar.title("Navigation")
    page = st.sidebar.radio(
        "Go to",
        ["Home", "Individual Analysis", "Strain Comparison", "Multi-Strain Comparison"],
        label_visibility="collapsed",
    )
    st.sidebar.markdown("---")

    if page == "Home":
        render_home_page()
    elif page == "Individual Analysis":
        render_individual_analysis_page()
    elif page == "Strain Comparison":
        render_strain_comparison_page()
    elif page == "Multi-Strain Comparison":
        render_multi_strain_comparison_page()


if __name__ == "__main__":
    main()
