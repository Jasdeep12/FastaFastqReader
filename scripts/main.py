"""Example usage of the FASTA/FASTQ reader."""

from fasta_parser import FASTAParser
from fastq_parser import FASTQParser
from filters import SequenceFilter
from stats import SequenceStats



def analyze_fasta(filepath):
    """Parse and analyze a FASTA file"""
    
    parser = FASTAParser(filepath)
    sequences = list(parser.parse())
    
    print(f"Parsing FASTA file @ {filepath}")
    print()
    
    stats = SequenceStats.length_stats(sequences)
    print("-- Length Statistics --")
    print(f"Total Sequences: {stats['count']}")
    print(f"Average Length: {stats['mean']}")
    print(f"Min Length: {stats['min']}")
    print(f"Max Length: {stats['max']}")
    print(f"Median Length: {stats['median']}")
    print(f"Standard Deviation of Length: {stats['stdev']}")
    print()
    
    stats = SequenceStats.gc_content_stats(sequences)
    print("-- GC Content --")
    print(f"Total Sequences: {stats['count']}")
    print(f"Average GC Content: {stats['mean']}")
    print(f"Min GC Content: {stats['min']}")
    print(f"Max GC Content: {stats['max']}")
    print(f"Median GC Content: {stats['median']}")
    print(f"Standard Deviation of GC Content: {stats['stdev']}")
    print()
    
    stats = SequenceStats.composition_stats(sequences)
    print("-- Composition --")
    for key, value in stats.items():
        print(f"{key} : {value}")
    print()
        
        
        
def analyze_fastq(filepath):
    """Parse and analyze a FASTQ file"""
    
    parser = FASTQParser(filepath)
    sequences = list(parser.parse())
    
    print(f"Parsing FASTQ file @ {filepath}")
    print()
    
    stats = SequenceStats.length_stats(sequences)
    print("-- Length Statistics --")
    print(f"Total Sequences: {stats['count']}")
    print(f"Average Length: {stats['mean']}")
    print(f"Min Length: {stats['min']}")
    print(f"Max Length: {stats['max']}")
    print(f"Median Length: {stats['median']}")
    print(f"Standard Deviation of Length: {stats['stdev']}")
    print()
    
    stats = SequenceStats.gc_content_stats(sequences)
    print("-- GC Content --")
    print(f"Total Sequences: {stats['count']}")
    print(f"Average GC Content: {stats['mean']}")
    print(f"Min GC Content: {stats['min']}")
    print(f"Max GC Content: {stats['max']}")
    print(f"Median GC Content: {stats['median']}")
    print(f"Standard Deviation of GC Content: {stats['stdev']}")
    print()
    
    stats = SequenceStats.quality_stats(sequences)
    print("-- Quality Stats --")
    print(f"Total Sequences: {stats['count']}")
    print(f"Average Quality: {stats['avg_quality']}")
    print(f"Min Quality: {stats['min_quality']}")
    print(f"Max Quality: {stats['max_quality']}")
    print()
    
    stats = SequenceStats.composition_stats(sequences)
    print("-- Composition --")
    for key, value in stats.items():
        print(f"{key} : {value}")
    print()
    
def filter_example(filepath, type : str):
    parser = None
    if type.upper() == 'FASTA':
        parser = FASTAParser(filepath)
    else:
        parser = FASTQParser(filepath)
    
    filtered = SequenceFilter.by_length(parser.parse(), min_length=100, max_length=500)
    print("-- Filter Example --")
    for seq in filtered:
        print(seq)
    

def main():
    """Main entry point for examples"""
    example = 'test.fastq'
    
    analyze_fastq(example)
    filter_example(example, example.partition('.')[2])


if __name__ == "__main__":
    main()