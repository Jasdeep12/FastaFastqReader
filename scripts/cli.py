"""A command line interface to use the FASTA/FASTQ reader using Click"""


import click
from fastq_parser import FASTQParser
from fasta_parser import FASTAParser
from kmer import KmerAnalyzer
from stats import SequenceStats
from filters import SequenceFilter
from trimmer import SequenceTrimmer

@click.group
def cli():
    """FASTA/FASTQ sequence analysis and processing tool"""
    pass


@cli.command()
@click.option('--fasta', type=click.Path(exists=True), help='FASTA file to analyze')
@click.option('--fastq', type=click.Path(exists=True), help='FASTQ file to analyze')
@click.option('--stats', multiple=True, default=['length','gc'], type=click.Choice(['length','gc','composition','quality']), help='Statistics to calculate')
def analyze(fasta, fastq, stats):
    """Analyze a FASTA or FASTQ file"""
    
    if not fasta and not fastq:
        click.echo("Error: Provide either --fasta or --fastq", err=True)
        return

    if fasta and fastq:
        click.echo("Error: Provide only one of --fasta or --fastq", err=True)
        return

    if fasta and 'quality' in stats:
        click.echo("Error: Quality statistics are only available for FASTQ files", err=True)
        return
    
    
    if fasta:
        click.echo(f" Analyzing FASTA: {fasta}\n")
        parser = FASTAParser(fasta)
        sequences = list(parser.parse())
        file_type = 'FASTA'
    else: 
        click.echo(f" Analyzing FASTQ: {fastq}\n")
        parser = FASTQParser(fastq)
        sequences = list(parser.parse())
        file_type = 'FASTQ'
    
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
      

@cli.command()
@click.option('--fasta', type=click.Path(exists=True), help='FASTA file to filter')
@click.option('--fastq', type=click.Path(exists=True), help='FASTQ file to filter')
@click.option('--min-length', type=int, default=None, help='Minimum sequence length')
@click.option('--max-length', type=int, default=None, help='Maximum sequence length')
@click.option('--min-gc', type=float, default=None, help='Minimum GC content (0-1)')
@click.option('--max-gc', type=float, default=None, help='Maximum GC content (0-1)')
@click.option('--min-quality', type=int, default=None, help='Minimum average quality score (FASTQ only)')
@click.option('--output', type=click.Path(), help='Output file (optional)')
def filter(fastq, fasta, min_length, max_length, min_gc, max_gc, min_quality, output):
    """Filter sequences by length and gc content"""

    if fasta and fastq:
        click.echo("Error: Provide only one of --fasta or --fastq", err=True)
        return
    if fasta and min_quality is not None:
        click.echo("Error: Quality filtering is only available for FASTQ files", err=True)
        return

    if fasta:
        click.echo(f"-- Filtering {fasta} --")
        parser = FASTAParser(fasta)
    elif fastq:
        click.echo(f"-- Filtering {fastq} --")
        parser = FASTQParser(fastq)
    else:
        click.echo("Error: Provide either --fasta or --fastq", err=True)
        return

    
    
    sequences = list(parser.parse())
    if not sequences:
        click.echo("No sequences found in the file", err=True)
        return

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
    for seq in sequences:
        count += 1
        click.echo(f"{seq.identifier}, {seq.length}bp, GC={seq.gc_content:.2%}")
        
    if output:

        with open(output, 'w') as f:
            for seq in sequences:
                if fasta:
                    f.write(f">{seq.identifier}\n{seq.sequence}\n")
                elif fastq:
                    f.write(f"@{seq.identifier}\n{seq.sequence}\n+\n{seq.quality}\n")

    click.echo(f"\n Found {count} sequences matching criteria")
    if output:
        click.echo(f" Written to {output}")
    
 
@cli.command()
@click.option('--fastq', type=click.Path(exists=True), required=True, help='FASTQ file to trim')
@click.option('--quality', type=int, default=20, help='Minimum quality threshold')
@click.option('--side', type=click.Choice(['left','right','both']), default='right', help='Which end to trim')
@click.option('--output', type=click.Path(), help='Output file (optional)')
def trim(fastq, quality, side, output):
    """Trim low quality bases from FASTQ sequences"""

    if not fastq:
        click.echo("Error: Provide a FASTQ file with --fastq", err=True)
        return
    
    click.echo(f"-- Trimming: {fastq} --")
    click.echo(f"Quality Threshold: {quality}, Side: {side}")
    
    parser = FASTQParser(fastq)
    sequences = list(parser.parse())
    
    if not sequences:
        click.echo("No sequences found in the file", err=True)
        return
    
    trimmedCount = 0
    totalBefore = 0
    totalAfter = 0
    trimmedlist = []
    for seq in sequences:
        if seq.length == 0:
            click.echo(f"{seq.identifier}: Sequence is empty, skipping")
            trimmedlist.append(seq) 
            continue
        totalBefore += seq.length
        trimmed =  SequenceTrimmer.trim_quality(seq, quality, side)
        totalAfter += trimmed.length
        trimmedCount += 1
        trimmedlist.append(trimmed)
        
        if seq.length != trimmed.length:
            click.echo(f"{trimmed.identifier}: {seq.length}bp -> {trimmed.length}")
    

    
    click.echo(f"\n Trimmed {trimmedCount} sequences")
    if totalBefore > 0:
        click.echo(f"  Total bases: {totalBefore} -> {totalAfter} ({100*totalAfter/totalBefore:.1f}%)")
    else:
        click.echo("  Total bases: 0 -> 0 (0.0%)")
    
    if output:
        with open(output, 'w') as f:
            for seq in trimmedlist:
                f.write(f"@{seq.identifier}\n{seq.sequence}\n+\n{seq.quality}\n")
        click.echo(f" Written to {output}")
        
        

@cli.command()
@click.option('--fasta', type=click.Path(exists=True), help='FASTA file')
@click.option('--fastq', type=click.Path(exists=True), help='FASTQ file')
@click.option('--k', type=int, default=3, help='Size of k-mers')
@click.option('--top',type=int, default=10, help='Show top N k-mers')
@click.option('--rare', type=int, default=0, help='Show rare N k-mers')
@click.option('--diversity', is_flag=True, help='Show k-mer diversity')
def kmer(fasta, fastq, k, top, rare, diversity):
    """Analyze k-mers frequencies in sequences"""

    if k <= 0:
        click.echo("Error: k must be a positive integer", err=True)
        return
    
    if not fasta and not fastq:
        click.echo("Error: Provide either --fasta or --fastq", err=True)
        return

    if fasta and fastq:
        click.echo("Error: Provide only one of --fasta or --fastq", err=True)
        return

    if fastq:
        click.echo(f"-- Analyzing k-mers in {fastq} --")
        parser = FASTQParser(fastq)
    else:
        click.echo(f"-- Analyzing k-mers in {fasta} --")
        parser = FASTAParser(fasta)

    sequences = list(parser.parse())
     
    if not sequences:
        click.echo("No sequences found in the file", err=True)
        return

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
@click.option('--fasta', type=click.Path(exists=True), help='FASTA file')
@click.option('--fastq', type=click.Path(exists=True), help='FASTQ file')
def info(fasta, fastq):
    """Show quick info about a sequence file"""
    
    if not fasta and not fastq:
        click.echo("Error: Provide either --fasta or --fastq", err=True)
        return

    if fasta and fastq:
        click.echo("Error: Provide only one of --fasta or --fastq", err=True)
        return
    
    if fasta:
        parser = FASTAParser(fasta)
        seqs = list(parser.parse())
        click.echo(f" FASTA File: {fasta}")
        
    elif fastq:
        parser = FASTQParser(fastq)
        seqs = list(parser.parse())
        click.echo(f" FASTQ File: {fastq}")
    
    click.echo(f"   Sequences: {len(seqs)}")
    
    if seqs:
        lengths = [s.length for s in seqs]
        click.echo(f"   Length range: {min(lengths)}-{max(lengths)}bp")
        click.echo(f"   First seq: {seqs[0].identifier} ({seqs[0].length}bp)")


if __name__ == '__main__':
    cli()