"""Shared helpers for the SI biology scripts (clustering, rare/common, FDA, event time, GO release).

Every path is resolved relative to the repository root, so scripts run from any working
directory. Nothing here reads outside the repository tree except ``DISATLAS_MODEL_DIR`` (an
optional override pointing at an unpacked copy of the released model package).
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw"
DERIVED = REPO / "data" / "derived"
ZENODO = REPO / "data" / "zenodo"
ZENODO_SIBIO = ZENODO / "si_bio"
MODEL_NAME = "disatlas_v3.2_ensemble"
SEEDS = (4707, 4708, 4709, 4710, 4711)

try:  # shared blinding helper (written by the data track)
    from disatlas.blinding import blind_target  # type: ignore  # noqa: F401
except Exception:  # pragma: no cover - import guard until the helper lands
    def blind_target(raw_id: str, salt: str = "DisAtlas") -> str:
        """Fallback stable token; replaced by disatlas.blinding at integration."""
        return "T_" + hashlib.sha256((salt + "|" + str(raw_id)).encode()).hexdigest()[:12]


def derived(topic: str) -> Path:
    """data/derived/<topic> (precomputed result tables shipped in git)."""
    p = DERIVED / topic
    p.mkdir(parents=True, exist_ok=True)
    return p


def zenodo_sibio(name: str = "") -> Path:
    """data/zenodo/si_bio/<name> (large inputs and caches shipped on Zenodo)."""
    ZENODO_SIBIO.mkdir(parents=True, exist_ok=True)
    return ZENODO_SIBIO / name if name else ZENODO_SIBIO


def model_dir() -> Path:
    """Released v3.2 model package (5-seed ensemble). Zenodo layout: data/zenodo/model/<MODEL_NAME>."""
    env = os.environ.get("DISATLAS_MODEL_DIR")
    cands = [Path(env)] if env else []
    from disatlas.paths import model_dir as _md
    cands += [_md(), _md().parent]
    for c in cands:
        if (c / "nodes.txt").exists() and (c / "embeddings").is_dir():
            return c
    raise FileNotFoundError(
        "released model not found; unpack the Zenodo model bundle into data/zenodo/model/ "
        "or set DISATLAS_MODEL_DIR")


def save_fig(fig, script_file: str, stem: str | None = None, dpi: int = 300):
    """Save png+pdf next to the calling script; return the two paths."""
    here = Path(script_file).resolve().parent
    stem = stem or Path(script_file).stem
    png, pdf = here / f"{stem}.png", here / f"{stem}.pdf"
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    return png, pdf


def write_numbers(script_file: str, numbers: dict, stem: str | None = None) -> Path:
    """Write the numbers a script reports to <stem>_numbers.json next to it."""
    here = Path(script_file).resolve().parent
    stem = stem or Path(script_file).stem
    out = here / f"{stem}_numbers.json"
    out.write_text(json.dumps(numbers, indent=2, sort_keys=True, default=float) + "\n")
    return out


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
