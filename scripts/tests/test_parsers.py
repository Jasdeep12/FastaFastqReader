import sys
from pathlib import Path
import pytest

TEST_DIR = Path(__file__).resolve().parent
sys.path.append(str(TEST_DIR.parent))
DATA_DIR = TEST_DIR / "test_data"

from reader import FileReader
from sequence import Sequence
from fasta_parser import FASTAParser
from fastq_parser import FASTQParser


def test_fastq_parser():
    fastq_file = DATA_DIR / "test.fastq"
    parser = FASTQParser(fastq_file)
    sequences = list(parser.parse())
    assert len(sequences) == 2
    assert sequences[0].header == "seq1"
    assert sequences[0].sequence == "AGCTTAGC"
    assert sequences[0].quality == "#15ADFGH"
    assert sequences[1].header == "seq2"
    assert sequences[1].sequence == "CGATCGAT"
    assert sequences[1].quality == "#<<A4FFF"
def test_fasta_parser():
    fasta_file = DATA_DIR / "test.fasta"
    parser = FASTAParser(fasta_file)
    sequences = list(parser.parse())
    assert len(sequences) == 7
    assert sequences[0].header == "Test_1"
    assert sequences[0].sequence == "TTCCACCGGCCCGGTGCAACTAAAAG"
    assert sequences[1].header == "Test_2"
    assert sequences[1].sequence == "TTTCGCAACGGCGTGATACCATCATC"

