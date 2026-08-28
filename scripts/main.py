"""Example usage of the FASTA/FASTQ reader."""

from fasta_parser import FASTAParser
from fastq_parser import FASTQParser
from filters import SequenceFilter
from stats import SequenceStats



def analyze_fasta(filepath):
    """Parse and analyze a FASTA file"""
    
    parser = FASTAParser(filepath)
    sequences = list(parser.parse())
    
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
    for key, value in stats:
        print(f"{key} : {value}")
        
        
def analyze_fastq(filepath):
    """Parse and analyze a FASTQ file"""
    
    parser = FASTAParser(filepath)
    
    
    
    
def filter_example():
    parser = FASTAParser('fastaReader-project/test/test.fasta')\
    
    filtered = SequenceFilter.by_length(parser.parse(), min_length=100, max_length=500)
    for seq in filterd:
        print(seq)
    





def main():
    """Main entry point for examples"""




if __name__ = "__main__":
    main()