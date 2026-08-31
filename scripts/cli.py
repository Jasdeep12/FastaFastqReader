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
      