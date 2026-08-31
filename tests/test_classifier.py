"""
tests/test_classifier.py
--------------------------
Tests for classifier.py: mutation category mapping, Ti/Tv classification,
DNA_Change formatting.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from classifier import classify_mutation, classify_snp_type, classify_row, dna_change


def test_synonymous_maps_to_silent():
    assert classify_mutation("synonymous_variant") == "Silent / Synonymous"


def test_missense_maps_correctly():
    assert classify_mutation("missense_variant") == "Missense"


def test_stop_gained_maps_to_nonsense():
    assert classify_mutation("stop_gained") == "Nonsense / Stop gained"


def test_frameshift_maps_correctly():
    assert classify_mutation("frameshift_variant") == "Frameshift"


def test_not_annotated_maps_to_unknown():
    assert classify_mutation("Not annotated") == "Unknown / Not annotated"
    assert classify_mutation(None) == "Unknown / Not annotated"


def test_unrecognized_effect_maps_to_other():
    assert classify_mutation("some_future_snpeff_term_not_in_dict") == "Other"


def test_combined_effect_uses_primary():
    # SnpEff sometimes joins multiple effects with '&'; primary (first) wins.
    assert classify_mutation("missense_variant&splice_region_variant") == "Missense"


def test_transitions():
    assert classify_snp_type("A", "G") == "Transition"
    assert classify_snp_type("G", "A") == "Transition"
    assert classify_snp_type("C", "T") == "Transition"
    assert classify_snp_type("T", "C") == "Transition"


def test_transversions():
    assert classify_snp_type("A", "C") == "Transversion"
    assert classify_snp_type("A", "T") == "Transversion"
    assert classify_snp_type("G", "C") == "Transversion"
    assert classify_snp_type("G", "T") == "Transversion"
    assert classify_snp_type("C", "A") == "Transversion"
    assert classify_snp_type("C", "G") == "Transversion"
    assert classify_snp_type("T", "A") == "Transversion"
    assert classify_snp_type("T", "G") == "Transversion"


def test_indel_not_misclassified_as_transition_or_transversion():
    assert classify_snp_type("C", "CA") == "Indel / Complex"
    assert classify_snp_type("CAG", "C") == "Indel / Complex"
    assert classify_snp_type(None, "A") == "Indel / Complex"
    assert classify_snp_type("A", None) == "Indel / Complex"


def test_dna_change_format():
    assert dna_change("C", "A") == "C>A"
    assert dna_change("T", "C") == "T>C"


def test_classify_row_combines_all_three():
    row = {"REF": "G", "ALT": "C", "Effect_Raw": "missense_variant"}
    result = classify_row(row)
    assert result["DNA_Change"] == "G>C"
    assert result["SNP_Type"] == "Transversion"
    assert result["Mutation_Category"] == "Missense"
