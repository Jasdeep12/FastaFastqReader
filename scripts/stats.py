"""Calculate statistics from sequences."""

from typing import Dict, List, Any
from sequence import Sequence
from statistics import mean, median, stdev

class SequenceStats:
    
    """Calculates various statistics from a collection of sequences."""
    
    
    @staticmethod
    def length_stats(sequences: List[Sequence]) -> Dict[str, float]:
        """Calculate length statistics (min, max, mean, median, stdev)"""
        if not sequences:
            raise ValueError("Sequences list is empty")
        
        lengths = [seq.length for seq in sequences]
         
        return {
            'count': len(lengths),
            'min' : min(lengths),
            'max' : max(lengths),
            'mean' : mean(lengths),
            'median' : median(lengths),
            'stdev' : stdev(lengths) if len(lengths) > 1 else 0,
        }
        
    @staticmethod
    def gc_content_stats(sequences: List[Sequence]) -> Dict[str, float]:
        """Calculate the GC content of each sequence"""
        
        if not sequences:
            raise ValueError("Sequences list is empty")
            
        gc_contents = [seq.gc_content for seq in sequences]
        
        return {
            'count' : len(gc_contents),
            'min' : min(gc_contents),
            'max' : max(gc_contents),
            'mean' : mean(gc_contents),
            'median' : median(gc_contents),
            'stdev' : stdev(gc_contents) if len(gc_contents) > 1 else 0,
        }
        
    @staticmethod
    def quality_stats(sequences: List[Sequence]) -> Dict[str, Any]:
        """Calculates quality score statistics (FASTQ only)"""
        
        if not sequences:
            raise ValueError("Sequences list is empty")
          
          
        qc_stats = []
        
        for seq in sequences:
            
            #filters out FASTA sequences
            if not seq.quality:
                continue
                
            qc_score = []
            
            for i in seq.quality:
                qc_score.append(ord(i)-33)
                
            qc_stats.append(qc_score)
        
        all_scores = [score for scores in qc_stats for score in scores]

        if not all_scores:
            raise ValueError("No quality scores found in the sequences")
        
        return {
            'count' : len(qc_stats),
            'min_quality' : min(all_scores),
            'max_quality' : max(all_scores),
            'avg_quality' : mean(all_scores),       
        }
    
    
    @staticmethod
    def composition_stats(sequences: List[Sequence]) -> Dict[str, Any]:
        """Calculate nucleotide composition percentages across all sequences"""
        
        if not sequences:
            raise ValueError("Sequence list is empty")
        
        total : float = 0.0
        
        counts = {}
        
        for seq in sequences:
            total += len(seq.sequence)
            
            for char in seq.sequence.upper():
                counts[char] = counts.get(char, 0) + 1
        
        
        return {char: round(count/total, 4) for char, count in counts.items()} if total > 0 else {}
                
