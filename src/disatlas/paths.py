"""Repository paths. Every script resolves inputs relative to the repository root, never to absolute paths."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW = DATA / "raw"
DERIVED = DATA / "derived"
ZENODO = DATA / "zenodo"          # external Zenodo bundle (not in git)
MANUSCRIPT_FIGURES = ROOT / "figures"


def derived(topic, *parts):
    return DERIVED.joinpath(topic, *parts)
