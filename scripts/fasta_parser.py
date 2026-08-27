"Parser for FASTA format files"
from typing import Iterator, Optional
from sequence import Sequence
from reader import FileReader

class FASTAParser:
    
    "Parses FASTA files and yields the Sequence objects"
    
    def __init__(self, filepath: str):
        "Initialize parser with FASTA file path"
        self.reader = FileReader(filepath)
        
    def parse(self) -> Iterator[Sequence]:
        """Parse FASTA file and yield Sequence objects.
        Headers begin with '>'"""
        
        
        current_header: Optional[str] = None
        current_sequence: list = []
        
        
        for line in self.reader.read_lines():
            if not line: 
                continue
            
            if line.startswith('>'):
                if current_header is not None:
                    yield Sequence(
                        header = current_header,
                        sequence = ''.join(current_sequence)                        
                    )
                current_header = line[1:]
                current_sequence = []
                
            else:
                if current_header is None:
                    raise ValueError(f"Sequence line found before header: {line}")
                current_sequence.append(line.strip())
            
        if current_header is not None:
            yield Sequence(
                header=current_header,
                sequence=''.join(current_sequence)
            )
        
        