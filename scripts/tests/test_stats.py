import sys
from pathlib import Path
import pytest

TEST_DIR = Path(__file__).resolve().parent
sys.path.append(str(TEST_DIR.parent))
DATA_DIR = TEST_DIR / "test_data"

from reader import FileReader
from sequence import Sequence
from stats import SequenceStats



def test_length_stats():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    stats = SequenceStats.length_stats(sequences)
    assert stats['count'] == 3
    assert stats['min'] == 2
    assert stats['max'] == 8
    assert stats['mean'] == pytest.approx(4.6667, rel=1e-4)
    assert stats['median'] == 4.0
    assert stats['stdev'] == pytest.approx(3.05505, rel=1e-4)

def test_gc_content_stats():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "GCGC"),
        Sequence(">seq3", "ATAT"),
    ]
    stats = SequenceStats.gc_content_stats(sequences)
    assert stats['count'] == 3
    assert stats['min'] == 0.0
    assert stats['max'] == 1.0
    assert stats['mean'] == pytest.approx(0.5, rel=1e-4)
    assert stats['median'] == 0.5
    assert stats['stdev'] == pytest.approx(0.5, rel=1e-4)

def test_quality_stats():
    sequences = [
        Sequence(">seq1", "ATCG", "#15A"),
        Sequence(">seq2", "GCGC", "#<<A"),
        Sequence(">seq3", "ATAT", None),  # This sequence has no quality scores
    ]
    stats = SequenceStats.quality_stats(sequences)
    assert stats['count'] == 2  # Only two sequences have quality scores
    assert stats['min_quality'] == 2
    assert stats['max_quality'] == 32
    assert stats['avg_quality'] == 19.75

def test_composition_stats():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "GCGC"),
        Sequence(">seq3", "ATAT"),
        Sequence(">seq4", "NGTC")
    ]
    stats = SequenceStats.composition_stats(sequences)
    assert stats['A'] == pytest.approx(.1875, rel=1e-3)
    assert stats['T'] == pytest.approx(.25, rel=1e-3)
    assert stats['G'] == pytest.approx(.25, rel=1e-3)
    assert stats['C'] == pytest.approx(.25, rel=1e-3)
    assert stats['N'] == pytest.approx(.0625, rel=1e-3)

def test_composition_stats_with_no_sequences():
    sequences = []
    with pytest.raises(ValueError):
        SequenceStats.composition_stats(sequences)

def test_quality_stats_with_no_quality_scores():
    sequences = [
        Sequence(">seq1", "ATCG", None),
        Sequence(">seq2", "GCGC", None),
    ]
    with pytest.raises(ValueError, match="No quality scores found in the sequences") as err:
        SequenceStats.quality_stats(sequences)

def test_length_stats_with_empty_list():
    sequences = []
    with pytest.raises(ValueError):
        SequenceStats.length_stats(sequences)

def test_gc_content_stats_with_empty_list():
    sequences = []
    with pytest.raises(ValueError):
        SequenceStats.gc_content_stats(sequences)

