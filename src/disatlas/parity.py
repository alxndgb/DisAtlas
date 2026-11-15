"""Helpers shared by the figure/table scripts and the parity tests.

Each figure or table script writes, next to its PNG/PDF (or MD/TSV), a ``<stem>_values.tsv`` with two columns
``key`` and ``value``: every number that the figure displays and the manuscript quotes. The parity tests in
``tests/`` compare these files with the manuscript's numbers (``tests/main_text_numbers.tsv``).
"""
import csv
from pathlib import Path


def write_values(path, values):
    """Write an ordered mapping key -> value (numbers or strings) as a two-column TSV."""
    path = Path(path)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["key", "value"])
        for k, v in values.items():
            w.writerow([k, repr(float(v)) if isinstance(v, float) else v])
    return path


def read_values(path):
    with open(path) as fh:
        return {r["key"]: r["value"] for r in csv.DictReader(fh, delimiter="\t")}


def image_similarity(png_a, png_b, size=(600, None)):
    """Mean absolute grey-level difference (0..1) after resizing both images to a common width.
    Used for visual parity against the manuscript image; renderer/font versions change pixels, not content."""
    import numpy as np
    from PIL import Image
    a = Image.open(png_a).convert("L"); b = Image.open(png_b).convert("L")
    w = size[0]; h = int(round(a.height * w / a.width))
    a = np.asarray(a.resize((w, h)), float) / 255; b = np.asarray(b.resize((w, h)), float) / 255
    return float(abs(a - b).mean()), (abs(a.shape[0] / a.shape[1] - b.shape[0] / b.shape[1]))
