import sys
from pathlib import Path
import pytest
from click.testing import CliRunner


TEST_DIR = Path(__file__).resolve().parent
sys.path.append(str(TEST_DIR.parent))
DATA_DIR = TEST_DIR / "test_data"

from reader import FileReader
from sequence import Sequence
from fasta_parser import FASTAParser
from fastq_parser import FASTQParser
from cli import cli, format_record


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

def test_fasta_roundtrip(tmp_path):
    original = list(FASTAParser(DATA_DIR / "test.fasta").parse())
    out = tmp_path / "rt.fasta"
    out.write_text("".join(format_record(s, "FASTA") for s in original))
    again = list(FASTAParser(out).parse())
    assert [(s.header, s.sequence) for s in again] == \
           [(s.header, s.sequence) for s in original]


def test_fastq_roundtrip(tmp_path):
    original = list(FASTQParser(DATA_DIR / "test.fastq").parse())
    out = tmp_path / "rt.fastq"
    out.write_text("".join(format_record(s, "FASTQ") for s in original))
    again = list(FASTQParser(out).parse())
    assert [(s.header, s.sequence, s.quality) for s in again] == \
           [(s.header, s.sequence, s.quality) for s in original]


def test_quality_stats_on_fasta_is_usage_error():
    r = CliRunner().invoke(cli, ["analyze", "--fasta", str(DATA_DIR / "test.fasta"),
                                 "--stats", "quality"])
    assert r.exit_code == 2
    assert "Error: Error:" not in r.output