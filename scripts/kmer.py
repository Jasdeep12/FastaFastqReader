"""K-mer analysis of a sequence"""
from sequence import Sequence
from collections import Counter



class KmerAnalyzer:
    
    
    
    @staticmethod
    def _sliding_window(sequence, k):
        """Yields all Kmers of length k from a sequence"""
        for i in range(len(sequence) - k + 1):
            yield sequence[i:i + k]
    
    
    @staticmethod
    def count_kmers(sequences: list[Sequence], k : int):
        """Counts k-mers per sequence for all sequences and keeps them in a dictionary"""
        if not sequences:
            raise ValueError("Sequences are empty")
        if k <= 0:
            raise ValueError("k must be a positive integer")
        
        kmers = Counter()
        
        for seq in sequences:
            string = seq.sequence
            kmers.update(KmerAnalyzer._sliding_window(string,k))
        
        return kmers
       
    @staticmethod
    def N_most_common(sequences: list[Sequence], k : int, top_n : int):
        """Return the top N most common kmers"""
        if top_n <= 0:
            raise ValueError("Top n must be a positive integer")
        
        kmers = KmerAnalyzer.count_kmers(sequences,k)
        
        return kmers.most_common(top_n)
        
    @staticmethod
    def find_rare_kmers(sequences, k , count):
        """Returns the least common kmers"""
        if count <= 0:
            raise ValueError("count must be a positive integer")
        kmers = KmerAnalyzer.count_kmers(sequences,k)
        
        return kmers.most_common()[:-count-1:-1]
        
        
    @staticmethod
    def kmer_diversity(sequences, k):
        """Returns the diversity of kmers in a sequence"""
        kmers = KmerAnalyzer.count_kmers(sequences,k)
        
        return len(kmers)
    