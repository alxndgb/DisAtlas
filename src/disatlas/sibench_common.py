"""Shared helpers for the SI benchmark and robustness scripts.

Every path is resolved relative to the repository root, so scripts run from any
working directory. Nothing here reads outside the repository tree.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "data" / "raw"
DERIVED = REPO / "data" / "derived"
ZENODO = REPO / "data" / "zenodo"
ZENODO_SIBENCH = ZENODO / "si_bench"

try:  # shared blinding helper (written by the data track)
    from disatlas.blinding import blind_target  # type: ignore  # noqa: F401
except Exception:  # pragma: no cover - import guard until the helper lands
    def blind_target(raw_id: str, salt: str = "DisAtlas") -> str:
        """Fallback stable token; replaced by disatlas.blinding at integration."""
        h = hashlib.sha256((salt + "|" + str(raw_id)).encode()).hexdigest()[:12]
        return "T_" + h


def derived(topic: str) -> Path:
    p = DERIVED / topic
    p.mkdir(parents=True, exist_ok=True)
    return p


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


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
