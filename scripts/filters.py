"""Options for filtering sequence by specific criteria"""

from typing import Iterator, Callable
from sequence import Sequence
from statistics import mean

class SequenceFilter:
    """Filter sequences based on specified criteria"""
    
    @staticmethod
    def by_length(sequences: Iterator[Sequence], min_length: int = 0, max_length: int = float('inf')) -> Iterator[Sequence]:
        """Filter Sequences based by length range"""
        if min_length is None:
            min_length = 0
        if max_length is None:
            max_length = float('inf')
        
        for seq in sequences:
            if min_length <= seq.length <= max_length:
                yield seq
    
    
    @staticmethod
    def by_gc_content(sequences: Iterator[Sequence], min_gc : float = 0.0, max_gc : float = 1.0) -> Iterator[Sequence]:
        """Filter Sequences based on GC content"""

        if min_gc is None:
            min_gc = 0.0
        if max_gc is None:
            max_gc = 1.0

            
        for seq in sequences:
            if seq.gc_content >= min_gc and seq.gc_content <= max_gc:
                yield seq
                

    @staticmethod
    def by_quality(sequences: Iterator[Sequence], min_avg_quality: int = 0) -> Iterator[Sequence]:
        """Filters Sequences based on quality"""
       
           
        for seq in sequences:
            
            #filters out FASTA sequences
            if not seq.quality:
                continue
                
            qc_score = []
            
            for i in seq.quality:
                qc_score.append(ord(i)-33)
            
            if mean(qc_score) >= min_avg_quality:
                yield seq
    
    
    @staticmethod
    def custom(sequences: Iterator[Sequence], predicate: Callable[[Sequence], bool]) -> Iterator[Sequence]:
        """Filter using a custom predicate function"""
        
        for seq in sequences:
            if predicate(seq):
                yield seq
        