"""A command line interface to use the FASTA/FASTQ reader using Click"""

from contextlib import nullcontext


import click
from fastq_parser import FASTQParser
from fasta_parser import FASTAParser
from kmer import KmerAnalyzer
from stats import SequenceStats
from filters import SequenceFilter
from trimmer import SequenceTrimmer

# Shared Helpers

def resolve_input(fasta,fastq):
    """Validates --fasta/--fastq and returns the parser, filetype and path"""

    if fasta and fastq:
        raise click.UsageError("Please provide either --fasta or --fastq, not both.")
    if not fasta and not fastq:
        raise click.UsageError("Please provide either --fasta or --fastq, not neither.")

    if fasta:
        return FASTAParser(fasta), "FASTA", fasta
    else:
        return FASTQParser(fastq), "FASTQ", fastq

def check_range(name, low, high):
    """Check given min and max to validate that min is not bigger than max"""
    if low is not None and high is not None and low > high:
        raise click.UsageError(f"--min-{name} ({low}) cannot be higher than --max-{name} ({high})")

def format_record(seq, file_type):
    """Serialize one record in it original form"""
    if file_type == "FASTA":
        return f">{seq.header}\n{seq.sequence}\n"
    else:
        return f"@{seq.header}\n{seq.sequence}\n+\n{seq.quality}\n"

def open_output(output):
    """Opens the output file, or a no-op context if no output was requested."""
    return open(output,"w") if output else nullcontext()

# Commands


@click.group
def cli():
    """FASTA/FASTQ sequence analysis and processing tool"""
    pass


@cli.command()
@click.option('--fasta', type=click.Path(exists=True, dir_okay=False), help='FASTA file to analyze')
@click.option('--fastq', type=click.Path(exists=True, dir_okay=False), help='FASTQ file to analyze')
@click.option('--stats', multiple=True, default=['length','gc'], type=click.Choice(['length','gc','composition','quality']), help='Statistics to calculate (repeatable)')
def analyze(fasta, fastq, stats):
    """Analyze a FASTA or FASTQ file"""
    parser, file_type, path = resolve_input(fasta, fastq)
    if file_type == "FASTA" and 'quality' in stats:
        raise click.UsageError("Quality statistics are only available for FASTQ files")
    
    sequences = list(parser.parse())
    click.echo(f"--Analyzing {file_type}: {path}\n")
    if not sequences:
        raise click.ClickException("No sequences found in the file.")

    click.echo(f"Total Sequences: {len(sequences)}\n")
    
    
    for stat in stats:
        if stat == 'length':
            result = SequenceStats.length_stats(sequences)
            click.echo("-- Length Statistics --")
            click.echo(f" Min: {result['min']:.4f}")
            click.echo(f" Max: {result['max']:.4f}")
            click.echo(f" Mean: {result['mean']:.4f}")
            click.echo(f" Median: {result['median']:.4f}")
            click.echo(f" Stdev: {result['stdev']:.4f}\n")
            
        elif stat == 'gc':
            result = SequenceStats.gc_content_stats(sequences)
            click.echo("-- GC Content Statistics --")
            click.echo(f" Min: {result['min']:.4f}")
            click.echo(f" Max: {result['max']:.4f}")
            click.echo(f" Mean: {result['mean']:.4f}")
            click.echo(f" Median: {result['median']:.4f}\n")
        
        elif stat == 'composition':
            result = SequenceStats.composition_stats(sequences)
            click.echo("-- Composition --")

            for base, freq in sorted(result.items()):
                click.echo(f"{base}: {freq:.4f}")
                
            click.echo()
        
        
        elif stat == 'quality':
            result = SequenceStats.quality_stats(sequences)
            click.echo("-- Quality Statistics --")
            click.echo(f" Min: {result['min_quality']:.4f}")
            click.echo(f" Max: {result['max_quality']:.4f}")
            click.echo(f" Avg: {result['avg_quality']:.4f}")
      

@cli.command(name='filter')
@click.option('--fasta', type=click.Path(exists=True, dir_okay=False), help='FASTA file to filter')
@click.option('--fastq', type=click.Path(exists=True, dir_okay=False), help='FASTQ file to filter')
@click.option('--min-length', type=click.IntRange(min=0), default=None, help='Minimum sequence length')
@click.option('--max-length', type=click.IntRange(min=0), default=None ,help='Maximum sequence length')
@click.option('--min-gc', type=click.FloatRange(0,1), default=None, help='Minimum GC content (0-1)')
@click.option('--max-gc', type=click.FloatRange(0,1), default=None, help='Maximum GC content (0-1)')
@click.option('--min-quality', type=click.IntRange(min=0), default=None, help='Minimum average quality score (FASTQ only)')
@click.option('--output', type=click.Path(writable=True, dir_okay=False), help='Output file (optional)')
def filter_cmd(fastq, fasta, min_length, max_length, min_gc, max_gc, min_quality, output):
    """Filter sequences by length and gc content"""

    parser, file_type, path = resolve_input(fasta,fastq)
    if file_type == "FASTA" and min_quality is not None:
        raise click.UsageError("Minimum quality available only with FASTQ files")   
    check_range("length",min_length,max_length)
    check_range("gc", min_gc, max_gc)

    
    sequences = list(parser.parse())
    if not sequences:
        raise click.ClickException("No sequences found in the file")
        

    click.echo(f"-- Filtering {path} --")

    if min_length is not None or max_length is not None:
        click.echo(f"Filtering by length: min={min_length}, max={max_length}")
        sequences = list(SequenceFilter.by_length(sequences, min_length, max_length))
    if min_gc is not None or max_gc is not None:
        click.echo(f"Filtering by GC content: min={min_gc}, max={max_gc}")
        sequences = list(SequenceFilter.by_gc_content(sequences, min_gc, max_gc))
    if min_quality is not None:
        click.echo(f"Filtering by quality: min_avg_quality={min_quality}")
        sequences = list(SequenceFilter.by_quality(sequences, min_quality))

    count = 0
    with open_output(output) as out:
        for seq in sequences:
            count += 1
            click.echo(f"{seq.header}, {seq.length}bp, GC={seq.gc_content:.2%}")
            if out:
                out.write(format_record(seq,file_type))

    click.echo(f"\n Found {count} sequences matching criteria")
    if output:
        click.echo(f" Written to {output}")
    
 
@cli.command()
@click.option('--fastq', type=click.Path(exists=True, dir_okay=False), required=True, help='FASTQ file to trim')
@click.option('--quality', type=click.IntRange(min=0), default=20, help='Minimum quality threshold')
@click.option('--side', type=click.Choice(['left','right','both']), default='right', help='Which end(s) to trim')
@click.option('--output', type=click.Path(writable=True,dir_okay=False), help='Output file (optional)')
def trim(fastq, quality, side, output):
    """Trim low quality bases from FASTQ sequences"""
    if not fastq:
        raise click.ClickException("Provide a FASTQ file with --fastq")

    processed = trimmed_count = bases_before = bases_after = 0

    click.echo(f"-- Trimming: {fastq} --")
    click.echo(f"Quality Threshold: {quality}, Side: {side}")

    sequences = list(FASTQParser(fastq).parse())

    if not sequences:
        raise click.ClickException("No sequences found in the file")

    with open_output(output) as out:
        for seq in sequences:
            processed += 1
            if seq.length == 0:
                click.echo(f"{seq.header}: empty sequence, passed through unchanged")

                result = seq
            else:
                result = SequenceTrimmer.trim_quality(seq,quality,side)
                bases_before += seq.length
                bases_after += result.length
                if result.length != seq.length:
                    trimmed_count += 1
                    click.echo(f"Trimmed {seq.length-result.length} bases from {seq.header}")
            if out:
                out.write(format_record(result,"FASTQ"))

    click.echo(f"\n Processed {processed} reads; {trimmed_count} were trimmed")
    if bases_before > 0:
        pct = 100 * bases_after/bases_before
        click.echo(f"Total bases: {bases_before} -> {bases_after} ({pct:.1f}% retained)")

    if output:
        click.echo(f"Written to {output}")
    

@cli.command()
@click.option('--fasta', type=click.Path(exists=True, dir_okay=False), help='FASTA file')
@click.option('--fastq', type=click.Path(exists=True, dir_okay=False), help='FASTQ file')
@click.option('--k', type=click.IntRange(min=1), default=3, show_default=True, help='Size of k-mers')
@click.option('--top',type=click.IntRange(min=1), default=10, show_default=True, help='Show top N k-mers')
@click.option('--rare', type=click.IntRange(min=0), default=0, show_default=True, help='Show rare N k-mers (0 = off)')
@click.option('--diversity', is_flag=True, help='Show k-mer diversity')
def kmer(fasta, fastq, k, top, rare, diversity):
    """Analyze k-mers frequencies in sequences"""

    parser, file_type, path = resolve_input(fasta,fastq)

    sequences = list(parser.parse())
     
    if not sequences:
        raise click.ClickException("No sequences found in the file")

    kmers = KmerAnalyzer.count_kmers(sequences, k)

    click.echo(f"Top {top} most common {k}-mers:")
    for kmer, count in kmers.most_common(top):
        click.echo(f"{kmer}: {count}")
    
    if rare > 0:
        click.echo(f"\nRare {k}-mers:")
        for kmer, count in KmerAnalyzer.find_rare_kmers(sequences, k, rare):
            click.echo(f"{kmer}: {count}")

    if diversity:
        diversity = KmerAnalyzer.kmer_diversity(sequences, k)
        click.echo(f"\nK-mer diversity: {diversity} unique {k}-mers")


@cli.command()
@click.option('--fasta', type=click.Path(exists=True, dir_okay=False), help='FASTA file')
@click.option('--fastq', type=click.Path(exists=True, dir_okay=False), help='FASTQ file')
def info(fasta, fastq):
    """Quick overview of a sequence file."""
    parser, file_type, path = resolve_input(fasta, fastq)
 
    count = 0
    first = None
    min_len = max_len = None
    for seq in parser.parse():  
        count += 1
        if first is None:
            first = seq
        min_len = seq.length if min_len is None else min(min_len, seq.length)
        max_len = seq.length if max_len is None else max(max_len, seq.length)
 
    click.echo(f"{file_type} file: {path}")
    click.echo(f"  Sequences: {count}")
    if first is not None:
        click.echo(f"  Length range: {min_len}-{max_len}bp")
        click.echo(f"  First seq: {first.header} ({first.length}bp)")
 


if __name__ == '__main__':
    cli()