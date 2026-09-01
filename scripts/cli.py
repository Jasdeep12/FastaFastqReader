"""A command line interface to use the FASTA/FASTQ reader using Click"""

import Click
from fastq_parser import FASTQParser
from fasta_parser import FASTAParser
from stats import SequenceStats
from filters import SequenceFilters
from trimmer import SequenceTrimmer

@click.group
def cli():
    """FASTA/FASTQ sequence analysis and processing tool"""
    pass:


@cli.command()
@click.option('--fasta', type=click.Path(exists=True), help='FASTA file to analyze')
@click.option('--fastq', type=click.Path(exists=True), help='FASTQ file to analyze')
@click.option('--stats', multiple=True, default=['length','gc'], type=.click.Choice(['length','gc','composition','quality']), help='Statistics to calculate')
def analyze(fasta, fastq, stats):
    """Analyze a FASTA or FASTQ file"""
    
    if not fasta and not fastq:
        click.echo("Error: Provide either --fasta or --fastq", err=True)
    
    
    if fasta:
        click.echo(f" Analyzing FASTA: {fasta}\n")
        parser = fasta_parser(fasta)
        sequences = list(parser.parse())
        file_type = 'FASTA'
    else: 
        click.echo(f" Analyzing FASTQ: {fastq}\n")
        parser = fastq_parser(fastq)
        sequences = list(parser.parse())
        file_type = 'FASTQ'
    
    click.echo(f"Total Sequences: {len(sequences)}\n")
    
    
    for stat in stats:
        if stat = 'length':
            result = SequenceStats.length_stats(sequences)
            click.echo("-- Length Statistics --")
            click.echo(f" Min: {result['min']:.4f}")
            click.echo(f" Max: {result['max']:.4f}")
            click.echo(f" Mean: {result['mean']:.4f}")
            click.echo(f" Median: {result['median']:.4f}")
            click.echo(f" Stdev: {result['stdev']:.4f}\n")
            
        elif stat = 'gc':
            result = SequenceStats.gc_content_stats(sequences)
            click.echo("-- GC Content Statistics --")
            click.echo(f" Min: {result['min']:.4f}")
            click.echo(f" Max: {result['max']:.4f}")
            click.echo(f" Mean: {result['mean']:.4f}")
            click.echo(f" Median: {result['median']:.4f}\n")
        
        elif stat = 'composition':
            result = SequenceStats.composition_stats(sequences)
            click.echo("-- Composition --")

            for base, freq in sorted(result.items()):
                click.echo(f"{base}: {freq:.4f}")
                
            click.echo()
        
        
        elif stat = 'quality':
            result = SequenceStats.composition_stats(sequences)
            click.echo("-- Quality Statistics --")
            click.echo(f" Min: {result['min_quality']:.4f}")
            click.echo(f" Max: {result['max_quality']:.4f}")
            click.echo(f" Avg: {result['avg_quality']:.4f}")
      

@cli.command()
@click.option('--fasta', type=click.Path(exists=True), required=True, help='FASTA file to filter')
@click.option('--min-length', type=int, default=0, help='Minimum sequence length')
@click.option('--max-length', type=int, default=0, help='Maximum sequence length')
@click.option('--min-gc', type=float, default=0.0, help='Minimum GC content (0-1)')
@click.option('--max-gc', type=float, default=0.0, help='Maximum GC content (0-1)')
@click.option('--output', type=click.Path(), help='Output file (optional)')
def filter(fasta, min_length, max_length, min_gc, max_gc, output):
    """Filter sequences by length and gc content"""
    
    click.echo(f"-- Filtering {fasta} --")
    parser = FASTAParser(fasta)
    
    
    sequences = list(parser.parse)
    sequences = SequenceFilters.by_length(sequences, min_length, max_length)
    sequences = SequenceFilters.by_gc_content(sequences, min_gc, max_gc)
    
    count = 0
    for seq in sequences:
        count += 1
        click.echo(f"{seq.identifier}, {seq.length}bp, GC={seq.gc_content:.2%}")
        
        if output:
            
            # TODO Write to file
           
           pass
    
    click.echo(f"\n Found {count} sequences matching criteria")
    
 
@cli.command()
@click.option('--fastq', type=click.Path(exists=True)