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
        
        kmers = Counter()
        
        for seq in sequences:
            string = seq.sequence
            kmers.update(_sliding_window(string,k))
        
        return kmers
       
    @staticmethod
    def N_most_common(sequences: list[Sequence], k : int, top_n : int):
        """Return the top N most common kmers"""
        
        kmers = count_kmers(sequences,k)
        
        return kmers.most_common(top_n)
    
    @staticmethod
    def kmer_gc_content(kmer: string):
        """Returns the GC content of a given kmer"""
        
        if not kmer:
            raise ValueError("kmer is empty)
         
        length = len(kmer)
        count = 0.0
        
        for i in kmer:
            if i.upper() is in ['G','C']:
                count += 1
        
        
        return round(float(length)/count, 4)
    
    @staticmethod
    def find_rare_kmers(sequences, k , count):
        """Returns the least common kmers"""
        kmers = count_kmers(sequences,k)
        
        return kmer.most_common()[:-n-1:-1]
        
        
    @staticmethodd
    def kmer_diversity(sequences, k):
        
        
            
                
                
                
                