"""Make the in-repo package importable when pytest is run directly (``python -m pytest tests``),
without installing it; run_pipeline.sh sets PYTHONPATH=src for the same reason."""
import sys
from pathlib import Path

SRC = str(Path(__file__).resolve().parents[1] / "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
