import sys
from pathlib import Path
import pytest
TEST_DIR = Path(__file__).resolve().parent
sys.path.append(str(TEST_DIR.parent))

from reader import FileReader
from sequence import Sequence
from filters import SequenceFilter



def test_sequence_filter_by_length():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    filtered = list(SequenceFilter.by_length(sequences, min_length=3))
    assert len(filtered) == 2
    assert filtered[0].header == ">seq1"
    assert filtered[0].sequence == "ATCG"
    assert filtered[1].header == ">seq2"
    assert filtered[1].sequence == "ATCGATCG"

def test_sequence_filter_by_length_with_none_values():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    filtered = list(SequenceFilter.by_length(sequences, min_length=None, max_length=None))
    assert len(filtered) == 3
    assert filtered[0].header == ">seq1"
    assert filtered[0].sequence == "ATCG"
    assert filtered[1].header == ">seq2"
    assert filtered[1].sequence == "ATCGATCG"
    assert filtered[2].header == ">seq3"
    assert filtered[2].sequence == "AT"

def test_sequence_filter_by_length_with_only_max_length():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    filtered = list(SequenceFilter.by_length(sequences, min_length=None, max_length=4))
    assert len(filtered) == 2
    assert filtered[0].header == ">seq1"
    assert filtered[0].sequence == "ATCG"
    assert filtered[1].header == ">seq3"
    assert filtered[1].sequence == "AT"

def test_sequence_filter_by_length_with_both_min_and_max_length():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    filtered = list(SequenceFilter.by_length(sequences, min_length=2, max_length=4))
    assert len(filtered) == 2
    assert filtered[0].header == ">seq1"
    assert filtered[0].sequence == "ATCG"
    assert filtered[1].header == ">seq3"
    assert filtered[1].sequence == "AT"

def test_sequence_filter_by_gc_content_with_none_values():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "GCGC"),
        Sequence(">seq3", "ATAT"),
    ]
    filtered = list(SequenceFilter.by_gc_content(sequences, min_gc=None, max_gc=None))
    assert len(filtered) == 3
    assert filtered[0].header == ">seq1"
    assert filtered[0].sequence == "ATCG"
    assert filtered[1].header == ">seq2"
    assert filtered[1].sequence == "GCGC"
    assert filtered[2].header == ">seq3"
    assert filtered[2].sequence == "ATAT"

def test_sequence_filter_by_gc_content():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "GCGC"),
        Sequence(">seq3", "ATAT"),
    ]
    filtered = list(SequenceFilter.by_gc_content(sequences, min_gc=0.5))
    assert len(filtered) == 2
    assert filtered[0].header == ">seq1"
    assert filtered[0].sequence == "ATCG"
    assert filtered[1].header == ">seq2"
    assert filtered[1].sequence == "GCGC"

def test_sequence_filter_by_quality():
    sequences = [
        Sequence(">seq1", "ATCG", quality="IIII"),
        Sequence(">seq2", "ATCGATCG", quality="!!!!"),
        Sequence(">seq3", "AT", quality="HHHH"),
    ]
    filtered = list(SequenceFilter.by_quality(sequences, min_avg_quality=30))
    assert len(filtered) == 2
    assert filtered[0].header == ">seq1"
    assert filtered[0].sequence == "ATCG"
    assert filtered[1].header == ">seq3"
    assert filtered[1].sequence == "AT"


def test_sequence_filter_by_quality_with_none_values():
    sequences = [
        Sequence(">seq1", "ATCG", quality="IIII"),
        Sequence(">seq2", "ATCGATCG", quality="!!!!"),
        Sequence(">seq3", "AT", quality="HHHH"),
    ]
    filtered = list(SequenceFilter.by_quality(sequences, min_avg_quality=None))
    assert len(filtered) == 3
    assert filtered[0].header == ">seq1"
    assert filtered[0].sequence == "ATCG"
    assert filtered[1].header == ">seq2"
    assert filtered[1].sequence == "ATCGATCG"
    assert filtered[2].header == ">seq3"
    assert filtered[2].sequence == "AT"

def test_sequence_filter_by_quality_with_fasta_sequences():
    sequences = [
        Sequence(">seq1", "ATCG", quality="IIII"),
        Sequence(">seq2", "ATCGATCG", quality=None),  # This is a FASTA sequence
        Sequence(">seq3", "AT", quality="HHHH"),
    ]
    filtered = list(SequenceFilter.by_quality(sequences, min_avg_quality=30))
    assert len(filtered) == 2  # Only two sequences have quality scores
    assert filtered[0].header == ">seq1"
    assert filtered[0].sequence == "ATCG"
    assert filtered[1].header == ">seq3"
    assert filtered[1].sequence == "AT"

def test_sequence_filter_by_quality_with_no_quality_scores():
    sequences = [
        Sequence(">seq1", "ATCG", quality=None),  # This is a FASTA sequence
        Sequence(">seq2", "ATCGATCG", quality=None),  # This is a FASTA sequence
    ]
    filtered = list(SequenceFilter.by_quality(sequences, min_avg_quality=30))
    assert len(filtered) == 0  # No sequences have quality scores

def test_sequence_filter_by_quality_with_none_sequences():
    with pytest.raises(ValueError):
        list(SequenceFilter.by_quality(None, min_avg_quality=30))

def test_sequence_filter_by_length_with_none_sequences():
    with pytest.raises(ValueError):
        list(SequenceFilter.by_length(None, min_length=3))

def test_sequence_filter_by_gc_content_with_none_sequences():
    with pytest.raises(ValueError):
        list(SequenceFilter.by_gc_content(None, min_gc=0.5))

def test_sequence_filter_custom_with_none_sequences():
    with pytest.raises(ValueError):
        list(SequenceFilter.custom(None, predicate=lambda seq: True))

def test_sequence_filter_custom():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "GCGC"),
        Sequence(">seq3", "ATAT"),
    ]
    filtered = list(SequenceFilter.custom(sequences, predicate=lambda seq: seq.gc_content > 0.5))
    assert len(filtered) == 1
    assert filtered[0].header == ">seq2"
    assert filtered[0].sequence == "GCGC"

