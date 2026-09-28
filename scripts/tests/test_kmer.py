import sys
from pathlib import Path
import pytest

TEST_DIR = Path(__file__).resolve().parent
sys.path.append(str(TEST_DIR.parent))
DATA_DIR = TEST_DIR / "test_data"

from sequence import Sequence
from kmer import KmerAnalyzer

def test_count_kmers():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    kmers = KmerAnalyzer.count_kmers(sequences, 2)
    assert len(kmers) == 4
    valid = ["AT","CG","GA","TC"]
    for kmer in kmers:
        assert kmer in valid

def test_N_most_common():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    kmers = KmerAnalyzer.N_most_common(sequences,2,3)
    assert len(kmers) == 3
    assert ("AT", 4) in kmers
    assert ("CG", 3) in kmers
    assert ("TC", 3) in kmers

def test_find_rare_kmers():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    kmers = KmerAnalyzer.find_rare_kmers(sequences,2,1)
    assert ("GA",1) == kmers[0]

def test_kmer_diversity():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    kmers = KmerAnalyzer.kmer_diversity(sequences,2)
    assert kmers == 4

def test_count_kmers_with_non_positive_int():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    with pytest.raises(ValueError):
        KmerAnalyzer.count_kmers(sequences,-1)

def test_count_kmers_with_no_sequences():
    sequences = []
    with pytest.raises(ValueError):
        KmerAnalyzer.count_kmers(sequences,2)

def test_count_kmers_with_non_positive_top_n():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    with pytest.raises(ValueError):
        KmerAnalyzer.N_most_common(sequences,2,-1)

def test_rare_kmers():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    with pytest.raises(ValueError):
        KmerAnalyzer.find_rare_kmers(sequences,2,-1)
