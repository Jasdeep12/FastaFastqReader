import sys
from pathlib import Path
import pytest

TEST_DIR = Path(__file__).resolve().parent
sys.path.append(str(TEST_DIR.parent))
DATA_DIR = TEST_DIR / "test_data"

from sequence import Sequence
from trimmer import SequenceTrimmer

def test_trim_quality_right():
    sequences = [
        Sequence(">seq1", "ATCGAC", "#15A##"),
        Sequence(">seq2", "GCGCAT", "#<<A##"),
    ]
    assert SequenceTrimmer.trim_quality(sequences[0],20,'right').sequence == "ATCG"
    assert SequenceTrimmer.trim_quality(sequences[1],20,'right').sequence == "GCGC"

def test_trim_quality_right():
    sequences = [
        Sequence(">seq1", "ATCGAC", "#15A##"),
        Sequence(">seq2", "GCGCAT", "#<<A##"),
    ]
    assert SequenceTrimmer.trim_quality(sequences[0],20,'left').sequence == "CGAC"
    assert SequenceTrimmer.trim_quality(sequences[1],20,'left').sequence == "CGCAT"