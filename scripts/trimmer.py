"""A Sequence trimmer to cut off bad quality reads"""
from sequence import Sequence


class SequenceTrimmer:
    
    
    @staticmethod
    def _trim_right(sequence : Sequence, min_quality = 20):
        """Helper function to trim the right side of the sequence"""
        for index, q in enumerate(reversed(sequence.quality)):
            if (ord(q) - 33) >= min_quality:
                #index = number of bad bases found at end
                if index > 0:
                    sequence.sequence = sequence.sequence[:-index]
                    sequence.quality = sequence.quality[:-index]
                return sequence
                
        #if entire sequence is of low quality
        sequence.quality = ""
        sequence.sequence = ""
        return sequence
        
    
    @staticmethod
    def _trim_left(sequence : Sequence, min_quality = 20):
        """Helper function to trim the left side of the sequence"""
        for index, q in enumerate(sequence.quality):
            if (ord(q)-33) >= min_quality:
                #index = number of bad bases found at start
                sequence.sequence = sequence.sequence[index:]
                sequence.quality = sequence.quality[index:]
                return sequence
        
        #if entire sequence is of low quality
        sequence.quality = ""
        sequence.sequence = ""
        return sequence
    
    
    @staticmethod
    def trim_quality(sequence : Sequence, min_quality = 20, side = 'right'):
        """Trim from specified ends until quality reaches a minimum threshold."""
        
        if not sequence:
            raise ValueError("Sequence is empty")
        if not sequence.quality:
            raise ValueError("Sequence must be a FASTQ sequence")
        if min_quality < 0:
            raise ValueError("Minimum quality cannot be less than 0")
        if side.upper() not in ['RIGHT','BOTH','LEFT']:
            raise ValueError("Side argument must be either: 'RIGHT', 'BOTH', 'LEFT'")
        
        
        if side.upper() == 'RIGHT':
            return SequenceTrimmer._trim_right(sequence, min_quality)
        if side.upper() == 'LEFT':
            return SequenceTrimmer._trim_left(sequence, min_quality)
        if side.upper() == 'BOTH':
            sequence = SequenceTrimmer._trim_left(sequence, min_quality)
            sequence = SequenceTrimmer._trim_right(sequence, min_quality)
            return sequence
        