"""Parser for FASTQ format files"""

from typing import Iterator, Optional
from sequence import Sequence
from reader import FileReader



class FASTQParser:
    
    """Parses FASTQ files and yields sequence objects.
    Handles the 4-line format of FASTQ (Header, sequence, plus, quality)."""
    
    def __init__(self, filepath:str):
        "initialize parser with a FASTQ file path."
        
        self.reader = FileReader(filepath)
        
    def parse(self) -> Iterator[Sequence]:
        
        """Parse FASTQ files and yield Sequence objects.
        FASTQ format : 4 line per record (@header, sequence, + , quality)"""
        
        lines = self.reader.read_lines()
        
        for header_line in lines:
            
            #line 1: header
            if not header_line.startswith('@'):
                raise ValueError(f"Expected @ header, got: {header_line}")
            
            #line 2: Sequence
            try: 
                sequence_line = next(lines)
            except StopIteration:
                raise ValueError("Truncated FASTQ: missing sequence")
            
            #line 3: plus
            try:
                plus_line = next(lines)
            except StopIteration:
                raise ValueError("Truncated FASTQ: missing plus")
                
            if not plus_line.startswith('+'):
                raise ValueError(f"Expected + line, got: {plus_line}")
            
            #line 4: quality
            try:
                quality_line = next(lines)
            except StopIteration:
                raise ValueError("Truncated FASTQ: missing plus")
            
            if len(quality_line) != len(sequence_line):
                raise ValueError(F"Quality Length {len(quality_line)} != {len(sequence_line)}")
                
            yield Sequence(
                header = header_line[1:],
                sequence = sequence_line,
                quality = quality_line
            )
                
                    
                