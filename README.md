# fastaReader

**A dependency-light Python toolkit and CLI for parsing, summarizing, filtering, quality-trimming, and k-mer profiling FASTA/FASTQ files.**

Built from scratch (no Biopython) to understand the formats end to end: streaming parsers, per-read statistics, Phred quality handling, and a Click command-line interface on top.

## Overview

```text
$ fastareader analyze --fasta tests/test_data/.fasta --stats length --stats gc 
```

```text
--Analyzing FASTA: tests/test_data/test.fasta

Total Sequences: 7

-- Length Statistics --
 Min: 26.0000
 Max: 1000.0000
 Mean: 462.7143
 Median: 559.0000
 Stdev: 347.5705

-- GC Content Statistics --
 Min: 0.4961
 Max: 0.5769
 Mean: 0.5174
 Median: 0.5062

```

---

## Features

| Command   | What it does |
|-----------|--------------|
| `info`    | Quick summary: record count, length range, first record |
| `analyze` | Length, GC content, base composition, and (FASTQ) quality statistics |
| `filter`  | Keep records by length, GC fraction, and (FASTQ) mean quality; optionally write results |
| `trim`    | Quality-trim FASTQ reads from the left, right, or both ends; reports reads trimmed and % bases retained |
| `kmer`    | Most common and rare k-mers, plus k-mer diversity (unique k-mer count) |

All commands accept either `--fasta` or `--fastq` (not both), except `trim`, which is FASTQ-only. Invalid input exits non-zero (usage errors → 2, runtime errors → 1), so commands are safe to chain in shell scripts.

## Usage

### `info`: first look at a file
```bash
python scripts/cli.py info --fasta test_data/test.fasta
```

### `analyze`: summary statistics
```bash
# Defaults to length + GC
python scripts/cli.py analyze --fasta genome.fasta

# Pick statistics explicitly (repeat --stats)
python scripts/cli.py analyze --fastq reads.fastq --stats length --stats composition --stats quality
```
Available statistics: `length`, `gc`, `composition`, `quality` (FASTQ only).

### `filter`: subset records
```bash
python scripts/cli.py filter --fastq reads.fastq \
    --min-length 50 --max-gc 0.6 --min-quality 25 \
    --output filtered.fastq
```
| Option | Meaning |
|--------|---------|
| `--min-length` / `--max-length` | Sequence length bounds (bp); min greater than max is rejected |
| `--min-gc` / `--max-gc` | GC content as a fraction, 0–1 (values outside this range are rejected) |
| `--min-quality` | Minimum mean Phred score (FASTQ only) |
| `--output` | Write passing records in the input's format |

### `trim`: quality trimming
```bash
python scripts/cli.py trim --fastq reads.fastq --quality 20 --side both --output trimmed.fastq
```
`--side` is `left`, `right` (default), or `both`. Bases below `--quality` (default 20) are removed from the chosen end(s).

### `kmer`: k-mer profiling
```bash
python scripts/cli.py kmer --fasta genome.fasta --k 5 --top 10 --rare 5 --diversity
```

---

## Project structure

```text
fastaReader-Project/scripts
├── cli.py            # Click CLI
├── reader.py         # Base file reader
├── sequence.py       # Sequence record (header, sequence, quality, GC, length)
├── fasta_parser.py   # FASTA parser (multi-line records)
├── fastq_parser.py   # FASTQ parser (4-line records)
├── stats.py          # Length / GC / composition / quality statistics
├── filters.py        # Length, GC, and quality filters
├── trimmer.py        # Quality trimming
├── kmer.py           # K-mer counting, rare k-mers, diversity
└── tests/
    ├── test_*.py
    └── test_data/
```

---

## Testing

```bash
pytest -v
```

---

## Design notes and limitations

- **Quality encoding:** assumes Phred+33 (Sanger / Illumina 1.8+). Older Phred+64 files will give wrong scores.
- **Memory:** the parsers are generators, and `info` streams records in a single pass. `analyze`, `filter`, `trim`, and `kmer` load the whole file into memory, which suits test-scale and small files but not multi-GB FASTQ. Streaming `filter` and `trim` is the next planned improvement.
- **Output:** records are written back in their input format with the original header line; FASTA sequences are written unwrapped (one sequence line per record).
- **Empty reads:** `trim` passes zero-length reads through unchanged rather than dropping them.
- **Scope:** single-file, single-end; no paired-end awareness or adapter trimming. Use `fastp` or `cutadapt` for production QC.

## Future work

- Making a proper pyproject.toml so it could installed as an actual package.
- Paired-end awareness
- Adapter trimming
- Support for Phred+64 scores
- Memory fixes
