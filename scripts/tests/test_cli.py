import sys
from pathlib import Path
import pytest

TEST_DIR = Path(__file__).resolve().parent
sys.path.append(str(TEST_DIR.parent))
DATA_DIR = TEST_DIR / "test_data"

from reader import FileReader
from sequence import Sequence



