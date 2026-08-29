from dataclasses import dataclass
from typing import Optional


@dataclass
class Sequence:
    """Represents a biological sequence (DNA, RNA, or protein)."""
    
    header: str
    sequence: str
    quality: Optional[str] = None
    
    @property
    def length(self) -> int:
        """Get sequence length."""
        return len(self.sequence)
    
    @property
    def identifier(self) -> str:
        """Extract sequence identifier (first whitespace-delimited token)."""
        return self.header.split()[0]
    
    @property
    def description(self) -> str:
        """Extract description (everything after first space in header)."""
        parts = self.header.split(None, 1)
        return parts[1] if len(parts) > 1 else ""
    
    @property
    def gc_content(self) -> float:
        """Calculate GC content (0-1)."""
        if not self.sequence:
            return 0.0
        gc_count = self.sequence.upper().count('G') + self.sequence.upper().count('C')
        return gc_count / len(self.sequence)
    
    def __repr__(self) -> str:
        q_info = f", quality={len(self.quality)}bp" if self.quality else ""
        return f"Sequence({self.identifier}, {self.length}bp{q_info})"