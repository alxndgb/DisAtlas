"""Repository paths. Every script resolves inputs relative to the repository root, never to absolute paths.

Run mode (env DISATLAS_DATA_MODE): ``warm`` (default) reads processed data (training graph, model, scores)
from the Zenodo bundle in data/zenodo; ``cold`` reads the same files rebuilt from data/raw into build/.
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RAW = DATA / "raw"
DERIVED = DATA / "derived"
ZENODO = DATA / "zenodo"          # external Zenodo bundle (not in git)
BUILD = ROOT / "build"            # cold-mode outputs (not in git)
MANUSCRIPT_FIGURES = ROOT / "figures"
MODE = os.environ.get("DISATLAS_DATA_MODE", "warm")
assert MODE in ("warm", "cold"), f"DISATLAS_DATA_MODE must be warm or cold, got {MODE!r}"
MODEL_NAME = "disatlas_v3.2_ensemble"


def derived(topic, *parts):
    return DERIVED.joinpath(topic, *parts)


def processed(*parts):
    """Processed file: data/zenodo/<parts> in warm mode, build/<parts> in cold mode."""
    return (BUILD if MODE == "cold" else ZENODO).joinpath(*parts)


def model_dir():
    return processed("model", MODEL_NAME)
