"""
classifier.py
-------------
Two independent classification systems:

1. Mutation category: maps raw SnpEff `Effect_Raw` terms to simplified,
   human-readable biological categories via an explicit dictionary
   (never substring matching).

2. SNP type: classifies REF>ALT nucleotide substitutions as Transition,
   Transversion, or Indel / Complex.

Both are pure functions operating on already-parsed data -- no VCF or
ANN parsing happens here.
"""

from __future__ import annotations

from utils import get_logger

logger = get_logger(__name__)

NOT_ANNOTATED = "Not annotated"


# 1. Mutation category mapping

""" 
    Keys are SnpEff/Sequence Ontology effect terms (as they literally appear
    in the ANN field). This list covers the standard SnpEff effect vocabulary;
    anything not listed here explicitly falls into "Other" via .get(), never
    via string matching.

""" 

EFFECT_TO_CATEGORY: dict[str, str] = {
    # Silent / Synonymous
    "synonymous_variant": "Silent / Synonymous",
    "start_retained_variant": "Silent / Synonymous",
    "stop_retained_variant": "Silent / Synonymous",

    # Missense
    "missense_variant": "Missense",
    "coding_sequence_variant": "Missense",

    # Nonsense / Stop gained
    "stop_gained": "Nonsense / Stop gained",

    # Stop lost
    "stop_lost": "Stop lost",

    # Start lost
    "start_lost": "Start lost",
    "initiator_codon_variant": "Start lost",

    # Frameshift
    "frameshift_variant": "Frameshift",

    # In-frame insertion / deletion
    "disruptive_inframe_insertion": "In-frame insertion",
    "conservative_inframe_insertion": "In-frame insertion",
    "inframe_insertion": "In-frame insertion",
    "disruptive_inframe_deletion": "In-frame deletion",
    "conservative_inframe_deletion": "In-frame deletion",
    "inframe_deletion": "In-frame deletion",

    # Splice-site
    "splice_acceptor_variant": "Splice-site",
    "splice_donor_variant": "Splice-site",
    "splice_region_variant": "Splice-site",
    "splice_branch_variant": "Splice-site",

    # Regulatory
    "regulatory_region_variant": "Regulatory",
    "TF_binding_site_variant": "Regulatory",
    "5_prime_UTR_premature_start_codon_gain_variant": "Regulatory",
    "5_prime_UTR_variant": "Regulatory",
    "3_prime_UTR_variant": "Regulatory",

    # Intergenic
    "intergenic_region": "Intergenic",

    # Upstream / Downstream
    "upstream_gene_variant": "Upstream",
    "downstream_gene_variant": "Downstream",

    # Intron
    "intron_variant": "Other",

    # Explicit "not annotated" passthrough from annotation.py
    NOT_ANNOTATED: "Unknown / Not annotated",
}


def classify_mutation(effect_raw: str | None) -> str:
    """
    Map a raw SnpEff effect string to a simplified mutation category.

    SnpEff sometimes reports combined effects joined with '&'
    (e.g. 'missense_variant&splice_region_variant'). In that case, we
    classify by the first (primary, most specific) listed effect, which
    is how SnpEff itself orders combined annotations.
    """
    if not effect_raw or effect_raw == NOT_ANNOTATED:
        return "Unknown / Not annotated"

    primary_effect = effect_raw.split("&")[0].strip()
    category = EFFECT_TO_CATEGORY.get(primary_effect)
    if category is None:
        logger.debug("Unmapped SnpEff effect '%s' -> classified as 'Other'", primary_effect)
        return "Other"
    return category


# 2. Transition / Transversion classification

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
TRANSVERSIONS = {
    ("A", "C"), ("C", "A"),
    ("A", "T"), ("T", "A"),
    ("G", "C"), ("C", "G"),
    ("G", "T"), ("T", "G"),
}
VALID_BASES = {"A", "C", "G", "T"}


def classify_snp_type(ref: str | None, alt: str | None) -> str:
    """
    Classify a REF>ALT substitution as 'Transition', 'Transversion', or
    'Indel / Complex' (for anything that isn't a clean single-nucleotide
    substitution -- insertions, deletions, multi-base REF/ALT, or
    missing/ambiguous alleles).
    """
    if not ref or not alt:
        return "Indel / Complex"

    ref = str(ref).strip().upper()
    alt = str(alt).strip().upper()

    if len(ref) != 1 or len(alt) != 1:
        return "Indel / Complex"
    if ref not in VALID_BASES or alt not in VALID_BASES:
        return "Indel / Complex"
    if ref == alt:
        # Not actually a substitution (shouldn't normally occur in a VCF,
        # but guard against it rather than mis-classifying).
        return "Indel / Complex"

    pair = (ref, alt)
    if pair in TRANSITIONS:
        return "Transition"
    if pair in TRANSVERSIONS:
        return "Transversion"
    return "Indel / Complex"  # unreachable given the base set, kept for safety


def dna_change(ref: str | None, alt: str | None) -> str:
    """Build the human-readable DNA_Change column, e.g. 'C>A'. Never
    modifies REF/ALT themselves."""
    r = ref if ref else "?"
    a = alt if alt else "?"
    return f"{r}>{a}"


def classify_row(row: dict) -> dict:
    """Return the classification columns to merge into a parsed+annotated row."""
    ref = row.get("REF")
    alt = row.get("ALT")
    return {
        "DNA_Change": dna_change(ref, alt),
        "SNP_Type": classify_snp_type(ref, alt),
        "Mutation_Category": classify_mutation(row.get("Effect_Raw")),
    }


def classify_rows(rows: list[dict]) -> list[dict]:
    """Apply classify_row to every row and merge results in-place."""
    for row in rows:
        row.update(classify_row(row))
    return rows
