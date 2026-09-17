# SNP Annotation Analyser

A local, modular SNP annotation analysis tool for annotated VCF files
(SnpEff `ANN` fields), built for bacterial genome resequencing workflows.
Runs entirely on your machine — no VCF data is ever uploaded to an
external server or API.

## What it does

- Parses annotated `.vcf` / `.vcf.gz` files (SnpEff `ANN` field)
- Classifies mutations biologically (silent, missense, nonsense, frameshift, ...)
- Classifies transitions vs. transversions
- Groups SNPs by gene (preserving every individual SNP)
- Groups SNPs by genomic position (across alternative alleles)
- Compares SNPs across multiple isolates (shared / isolate-specific)
- Produces an interactive Streamlit dashboard with Plotly visualizations
- Exports CSV, Excel, a standalone HTML report, and a publication-ready PDF
- Includes a command-line interface for batch/headless use

## Project structure

```
SNP_Analyser/
├── app.py              # Streamlit web interface
├── main.py             # Command-line interface
├── pipeline.py         # Shared parse -> annotate -> classify orchestration
├── vcf_parser.py       # cyvcf2-based VCF reading & validation
├── annotation.py       # SnpEff ANN field parsing
├── classifier.py       # Mutation category + Ti/Tv classification
├── grouping.py         # Gene-level and position-level grouping
├── statistics.py       # Summary statistics
├── comparison.py       # Multi-isolate comparison
├── visualizations.py   # Plotly chart builders (shared by app + report)
├── report.py           # HTML/PDF report generation
├── export.py           # CSV/Excel export
├── utils.py            # Logging + shared path helpers
├── templates/report.html
├── static/style.css
├── data/                # Place VCF files here for CLI use (optional)
├── output/              # Generated reports/exports land here
└── requirements.txt
```

## Installation (macOS, VS Code terminal / zsh)

1. Check Python version:
   ```zsh
   python3 --version
   ```
   You need Python 3.9 or newer (3.10/3.11 recommended).

2. From the `SNP_Analyser/` folder, create and activate a virtual environment:
   ```zsh
   python3 -m venv .venv
   source .venv/bin/activate
   ```

3. Install Xcode's Command Line Tools (needed to compile `cyvcf2`, which
   wraps htslib):
   ```zsh
   xcode-select --install
   ```
   Wait for the GUI installer to finish before continuing.

4. Install Python dependencies:
   ```zsh
   pip install --upgrade pip
   pip install -r requirements.txt
   ```
   `cyvcf2` compiles from source on first install — this can take a
   few minutes and is normal.

5. **WeasyPrint (needed for PDF report generation):** it depends on
   native libraries (Pango, cairo, gdk-pixbuf, libffi) that aren't
   pure-Python. Install Homebrew first if you don't have it
   (https://brew.sh), then:
   ```zsh
   brew install pango cairo gdk-pixbuf libffi
   ```
   If `pip install weasyprint` still fails afterward, close and reopen
   your terminal so Homebrew's paths are picked up, then retry.

6. **Chrome (needed for embedding chart images in the PDF report):**
   the PDF report renders each Plotly chart as a static image via
   `kaleido`, which needs a local Chrome install:
   ```zsh
   plotly_get_chrome
   ```
   You only need this once. If you skip it, the PDF report will still
   generate — chart sections will show a placeholder note instead of
   the image, but every table and statistic is unaffected. The
   interactive HTML report is unaffected either way (it embeds
   interactive Plotly.js charts, not static images).

## Running the Streamlit app

```zsh
streamlit run app.py
```

Your browser opens automatically. Use the sidebar to navigate between three pages:

**Home** — overview of what the app does.

**Individual Analysis** — the original single-VCF workflow:
1. Under "1. Upload VCF", click **Browse files** and select your
   `.vcf` or `.vcf.gz` file. You can select multiple files at once to
   compare isolates of the same strain.
2. Click **Parse & Analyse**.
3. File details, summary statistics, the interactive SNP table, Gene
   Explorer, Position Explorer, visualisations, and (if ≥2 samples are
   present) the multi-isolate comparison section all populate
   automatically.
4. Use the sidebar filters (gene, category, effect, SNP type, sample,
   QUAL/DP, text search) to narrow the table — every section below
   updates to match.
5. Scroll to **8. Full Report** and click **Generate Reports** to build
   the HTML and PDF reports, then use the download buttons.

**Strain Comparison** — compare two independent strains:
1. Upload Strain 1's VCF and Strain 2's VCF.
2. **Give each strain a name** (e.g. `strain1`, `strain2`). This matters:
   two different strains' VCFs can (and often do) use the exact same
   internal sample name — the app cannot tell them apart from the VCF
   alone, so the strain names you type here are what distinguishes them
   everywhere downstream. If left blank, the filename is used instead.
3. Click **Analyse & Compare**. Both strains are parsed and classified
   completely independently first; only then are they linked by genomic
   coordinate (CHROM + POS).
4. The comparison summary, filterable position-level comparison table,
   a **Position Detail View** (type in any position to see both strains'
   complete records side by side), gene-level comparison, and
   visualisations all populate automatically.
5. Use the sidebar filters (comparison status, gene, text search) to
   narrow the comparison table.
6. Scroll to **7. Comparison Report** to generate and download the
   comparison HTML/PDF report, or export individual comparison tables
   (shared SNPs, strain-specific SNPs, gene comparison) as CSV.

**Multi-Strain Comparison** — compare any number of strains (2, 10, 50+) at once:
1. Upload 2 or more VCF files at once (the file uploader accepts
   multiple files in one selection).
2. **Name each strain individually** — a text box appears for every
   uploaded file, pre-filled with a name derived from the filename. Edit
   any of these so every strain has a distinct, meaningful name (same
   reasoning as above: VCFs frequently share the same internal sample
   name across different strains).
3. Click **Analyse & Compare All Strains**. Every strain is parsed and
   classified completely independently first (same pipeline as
   Individual Analysis), then linked into one N-way comparison.
4. Instead of a single "identical/different" status (which only makes
   sense for exactly two strains), each genomic position is classified
   by how many of your N strains carry a variant there:
   - **Core, identical** — every strain has it, same allele
   - **Core, variable** — every strain has a variant there, but not
     the same allele
   - **Partial (k/N)** — some but not all strains have it (exact count shown)
   - **`<strain>`-unique** — only one strain has it
   - **Reference allele discrepancy** — strains disagree on REF itself
5. The comparison summary, filterable SNP matrix, N-way Position Detail
   View (shows every strain's record at a position side by side), and
   gene-level comparison (which strains are affected in each gene, and
   whether a gene is "core" — affected in all of them) all populate
   automatically.
6. Visualisations scale to any N: an NxN pairwise sharing heatmap, a
   core-genome-size summary, per-strain unique-position counts, a
   genes-by-strain-count histogram, and a presence/absence matrix
   (capped to the top 300 positions by strain count, for renderability).
7. Scroll to **7. Comparison Report** to generate the combined HTML/PDF
   report, or export the full SNP matrix, core-genome-only subset, or
   gene comparison as CSV.

## Running the CLI

**Single-VCF mode** (unchanged):
```zsh
python main.py --vcf data/your_file_name.vcf --output output/
```

**Two-strain comparison mode:**
```zsh
python main.py --vcf1 data/your_file_name1.vcf --vcf2 data/your_file_name2.vcf --name1 strain1 --name2 strain2 --output output/
```

**Multiple isolates of the same strain:**
```zsh
python main.py --vcf your_file_name1.vcf your_file_name2.vcf --output output/
```

`--name1`/`--name2` are optional — if omitted, strain names are derived
from the filenames (e.g. `data/your_file_name` → `strain1`). Provide them
explicitly whenever your filenames aren't descriptive, or when two VCFs
share the same internal sample name (common — the VCF's own `#CHROM`
header line often just says something generic like `isolate1` regardless
of which strain it actually is).

`--strain-names` is optional (defaults to filenames), but strongly
recommended for the same reason as two-strain mode: different strains'
VCFs frequently share the same internal sample name. Scales to any
number of files — list as many `--vcfs` paths and matching
`--strain-names` as you have strains.

Produces:
```
output/
├── strain1/          -- individual analysis, same structure as above
├── strain2/
├── strain3/
│   └── (one folder per strain, same files as two-strain mode)
├── multi_strain/
│   ├── <prefix>_snp_matrix.csv               -- full N-way position matrix
│   ├── <prefix>_core_genome_positions.csv    -- positions shared by ALL strains
│   ├── <prefix>_<strain>_unique.csv          -- one file per strain
│   ├── <prefix>_gene_comparison.csv
│   ├── <prefix>_pairwise_sharing_matrix.csv
│   └── <prefix>_comparison.html / .pdf
└── snp_analyser.log
```
(`<prefix>` is all strain names joined for 4 or fewer strains, or
`N strains` for larger sets, to keep filenames manageable at high N.)

## How classification works

- **Mutation category** is derived from the SnpEff `Effect_Raw` term
  (e.g. `missense_variant`) via an explicit lookup table in
  `classifier.py` (`EFFECT_TO_CATEGORY`) — never by substring matching.
  The original SnpEff term is always preserved alongside the simplified
  category. Records with no `ANN` field are labeled
  `Unknown / Not annotated` rather than dropped.
- **Transition / Transversion** classification only applies to clean
  single-nucleotide REF>ALT substitutions (A↔G, C↔T = transitions;
  all other single-base swaps = transversions). Anything that isn't a
  simple 1-base-to-1-base substitution (insertions, deletions, missing
  alleles) is labeled `Indel / Complex` and is never miscounted as a
  transition or transversion.

## How gene-level grouping works

Every SNP keeps its own row — grouping by gene (`grouping.group_by_gene`)
never merges individual mutations. Instead, it aggregates: total SNP
count, unique position count, per-category counts (missense/synonymous/
nonsense/frameshift/other), transition/transversion counts and ratio,
and list columns of every position, DNA change, protein change, and
effect within that gene. 

## How position-level grouping works

Position grouping (`grouping.group_by_position`) uses a different key —
`(CHROM, POS, REF)` — than gene grouping. This correctly recognises that
two different ALT alleles at the same genomic position 
belong to the *same* position, and reports each allele's
isolate count and percentage, without conflating this with gene-level
grouping or merging separate genomic positions together.

## How two-strain comparison works

Two-strain comparison is a separate, additional layer on top of the
existing single-VCF pipeline — it is not a replacement for it. Each
strain's VCF is parsed, annotated, and classified **completely
independently** first, using the exact same pipeline as single-VCF
analysis (`pipeline.process_vcf`), just tagged with a `Strain` column
(`pipeline.process_vcf_with_strain_name`) so two strains that happen to
share the same internal VCF sample name stay distinguishable.

Only after both strains' complete datasets exist does
`comparison.build_position_comparison()` link them, by genomic
coordinate (`CHROM` + `POS`). Neither original dataset is ever modified
by this step — the comparison table stores full original records from
both strains, not just position/REF/ALT. Every shared position is
classified into exactly one status:

- **Identical mutation** — same REF, same ALT, in both strains
- **Same position, different mutation** — same REF, different ALT
  (correctly distinguishes e.g. `G>A` in one strain from `G>C` in the
  other, rather than treating them as a match)
- **Reference allele discrepancy** — the REF itself differs between the
  two VCFs at that position, which usually signals the VCFs may have
  been called against different reference assemblies; the app warns
  about this without attempting to auto-correct it
- **`<strain>`-specific** — the position exists in only one strain's VCF

Multi-allelic sites are compared as *sets* of ALT alleles, not single
values, so a position with two ALT alleles in one strain is compared
correctly against the other strain's allele(s) rather than being
collapsed into one comparison.

Gene-level comparison is intentionally kept separate from position-level
comparison: two strains can share a gene (both have SNPs somewhere in
*gyrA*, for example) without sharing any specific genomic position within
that gene. `comparison.build_gene_comparison()` reports both shared
positions *and* shared genes as distinct numbers, and the app never
conflates the two.

## How multi-strain (N-strain) comparison works

Multi-strain comparison generalises the two-strain engine from exactly
two strains to any number (2, 10, 50+). Every strain is still parsed,
annotated, and classified completely independently first, using the
identical pipeline as single-VCF and two-strain analysis
(`pipeline.process_multiple_strains()`). Only after every strain's
complete dataset exists does `comparison.build_snp_matrix()` link them
by genomic coordinate.

The two-strain engine's five status categories don't generalise cleanly
past N=2 (e.g. "identical" vs. "different" doesn't describe a position
present in 7 of 10 strains), so the N-way engine classifies each
position by **how many of the N strains carry a variant there**:

- **Core, identical** — every strain has it, with the same allele
- **Core, variable** — every strain has a variant there, but not the
  same allele across all of them
- **Partial (k/N)** — some but not all strains have it, with the exact
  count and list of which strains shown
- **`<strain>`-unique** — only one strain has it
- **Reference allele discrepancy** — strains disagree on the REF allele
  itself at that position

Running the N-way engine with exactly 2 strains produces results
mathematically identical to the dedicated two-strain engine (this is
verified directly in `tests/test_multi_strain_engine.py`) — so the two
engines are consistent with each other, not two independent
implementations that could drift apart.

**Visualisations that scale with N:** the two-strain page's Venn-style
bar chart and 2-color grouped bars stop being readable past about 3-4
strains, so multi-strain comparison uses different visualisations
instead: an NxN pairwise sharing heatmap (`fig_pairwise_sharing_heatmap`),
a core-genome-size summary bar, per-strain unique-position counts, a
histogram of genes by how many strains they're affected in, and a
presence/absence matrix (capped to the top 300 positions by strain
count, since a full matrix for 10+ strains can easily exceed several
thousand positions and become unrenderable).

**Gene-level N-way comparison** (`comparison.build_gene_comparison_n_way()`)
reports, per gene, one SNP-count column per strain, how many strains are
affected in that gene at all, and whether the gene is "core" (affected
in every strain) — again, this is the direct generalisation of the
two-strain gene comparison's `Shared_Gene` flag.

**A note on performance:** for realistic bacterial genome VCFs (a few
thousand SNPs per strain), 10+ strains produces a SNP matrix with
roughly 5,000-10,000 unique positions and one column-pair per strain —
this is comfortably within pandas' normal operating range and completes
in well under a minute on a laptop. This hasn't been tested at hundreds
of strains, where the number of pairwise comparisons (`N × (N-1) / 2`)
for the sharing heatmap starts to grow meaningfully, though the core
SNP-matrix construction itself remains a single set of outer merges
regardless of N.

## How to generate HTML / PDF / CSV / Excel

- **HTML report:** self-contained (Plotly.js embedded inline, no
  internet needed to view it), generated via `report.generate_reports()`
  — available from the Streamlit app's "Full Report" section or
  automatically by the CLI.
- **PDF report:** same content, rendered via WeasyPrint with static
  chart images (needs the Chrome setup step above); suitable for
  supplementary material or lab records. Wide tables render in
  landscape sections.
- **CSV / Excel:** available as individual download buttons in the
  Streamlit app (SNP table, gene summary, position summary) or written
  automatically by the CLI into the output directory.

## Troubleshooting

| Problem | Fix |
|---|---|
| `xcrun: error: invalid active developer path` | Run `xcode-select --install`, wait for it to finish, retry `pip install`. |
| `pip install cyvcf2` hangs or is very slow | Normal — it's compiling from source. Give it a few minutes. |
| `weasyprint` import fails / PDF generation errors | Run `brew install pango cairo gdk-pixbuf libffi`, then reopen your terminal and retry. As of this version, a broken WeasyPrint install no longer crashes the whole run — it logs a warning and every other output (CSV, Excel, HTML reports) still completes normally. |
| PDF report shows "could not be rendered" for chart images | Run `plotly_get_chrome` once. The rest of the PDF (tables, stats) is unaffected either way. |
| `ModuleNotFoundError` when running `streamlit run app.py` or `python main.py` | Make sure your virtual environment is activated (`source .venv/bin/activate`) and you're running from inside the `SNP_Analyser/` folder. |
| Streamlit opens but upload fails / errors immediately | Check the error message shown in the app — the uploader validates file extension and content before parsing and will report exactly what's wrong (e.g. empty file, invalid VCF) rather than crashing silently. |
| A VCF record seems to be "missing" from results | It shouldn't be — every parsed record is tracked in the "Parsing integrity" panel (total / annotated / unannotated / failed). Check `output/snp_analyzer.log` for per-record warnings. |
