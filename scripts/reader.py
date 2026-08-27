"""Low-level file reader supporting both plain and gzip formats."""
import gzip
from pathlib import Path
from typing import Iterator


class FileReader:
    """
    Handles reading from FASTA/FASTQ files.
    Supports both plain text and gzip-compressed files.
    The reader auto-detects file compression based on file extension.
    """
    
    def __init__(self, filepath: str):
        """
        Initialize the reader with a file path.
        Args:
            filepath: Path to the file (can be .gz compressed) 
        Raises:
            FileNotFoundError: If the file doesn't exist
        """
        self.filepath = Path(filepath)
        self._validate_file()
    
    def _validate_file(self) -> None:
        """Check if file exists."""
        if not self.filepath.exists():
            raise FileNotFoundError(f"File not found: {self.filepath}")
    
    def _open_file(self):
        """
        Open file, auto-detecting if it's gzip-compressed.
        Returns:
            An open file object (text mode)
        """
        if self.filepath.suffix == '.gz':
            return gzip.open(self.filepath, 'rt')
        else:
            return open(self.filepath, 'r')
    
    def read_lines(self) -> Iterator[str]:
        """
        Yield lines from the file one at a time.
        Strips trailing newlines from each line.
        Handles both plain text and gzip-compressed files.
        Yields:
            str: Each line without the trailing newline
        """
        with self._open_file() as f:
            for line in f:
                yield line.rstrip('\n')
    
    def read_in_chunks(self, chunk_size: int = 1000) -> Iterator[list]:
        """
        Yield lines in chunks for memory-efficient processing.
        Useful for processing large files without loading everything in memory.
        Args:
            chunk_size: Number of lines per chunk (default: 1000)
        Yields:
            list: A list of strings (lines), up to chunk_size items
        """
        chunk = []
        for line in self.read_lines():
            chunk.append(line)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        
        # Yield remaining lines
        if chunk:
            yield chunk