import pytest


def test_length_stats():
    sequences = [
        Sequence(">seq1", "ATCG"),
        Sequence(">seq2", "ATCGATCG"),
        Sequence(">seq3", "AT"),
    ]
    stats = stats.length_stats(sequences)
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
    stats = Stats.gc_content_stats(sequences)
    assert stats['count'] == 3
    assert stats['min'] == 0.0
    assert stats['max'] == 1.0
    assert stats['mean'] == pytest.approx(0.5, rel=1e-4)
    assert stats['median'] == 0.5
    assert stats['stdev'] == pytest.approx(0.57735, rel=1e-4)

def test_quality_stats():
    sequences = [
        Sequence(">seq1", "ATCG", "#15ADFGH"),
        Sequence(">seq2", "GCGC", "#<<A4FFF"),
        Sequence(">seq3", "ATAT", None),  # This sequence has no quality scores
    ]
    stats = Stats.quality_stats(sequences)
    assert stats['count'] == 2  # Only two sequences have quality scores
    assert stats['min'] == pytest.approx(33.0, rel=1e-4)
    assert stats['max'] == pytest.approx(70.0, rel=1e-4)
    assert stats['mean'] == pytest.approx(51.5, rel=1e-4)
    assert stats['median'] == pytest.approx(51.5, rel=1e-4)
    assert stats['stdev'] == pytest.approx(18.3848, rel=1e-4)

def test_quality_stats_with_no_quality_scores():
    sequences = [
        scripts.Sequence(">seq1", "ATCG", None),
        scripts.Sequence(">seq2", "GCGC", None),
    ]
    with pytest.raises(ValueError):
        scripts.Stats.quality_stats(sequences)

def test_length_stats_with_empty_list():
    sequences = []
    with pytest.raises(ValueError):
        scripts.Stats.length_stats(sequences)

def test_gc_content_stats_with_empty_list():
    sequences = []
    with pytest.raises(ValueError):
        scripts.Stats.gc_content_stats(sequences)