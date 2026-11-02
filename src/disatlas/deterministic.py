"""Byte-reproducible vector output for every figure script.

Matplotlib writes the current time into PDF (/CreationDate) and SVG (<dc:date>) files and draws SVG element ids from a
random salt, so a rerun changes every .pdf/.svg even when the drawing is identical.  Importing this module (done by
`import disatlas` and at the top of each figure script) fixes svg.hashsalt and makes the PDF and SVG writers default to
metadata without a date; metadata passed explicitly by a script is kept.  PNG output is already deterministic.
"""
import functools

import matplotlib as mpl
from matplotlib.backends import backend_pdf, backend_svg

PDF_METADATA = {"CreationDate": None}
SVG_METADATA = {"Date": None}
SVG_HASHSALT = "disatlas"


def _wrap(canvas, name, default):
    orig = getattr(canvas, name)
    if getattr(orig, "_disatlas_deterministic", False):
        return

    @functools.wraps(orig)  # keeps the signature, which matplotlib inspects to filter savefig keywords
    def print_vector(self, filename, *args, metadata=None, **kw):
        return orig(self, filename, *args, metadata={**default, **(metadata or {})}, **kw)

    print_vector._disatlas_deterministic = True
    setattr(canvas, name, print_vector)


def install():
    mpl.rcParams["svg.hashsalt"] = SVG_HASHSALT
    _wrap(backend_pdf.FigureCanvasPdf, "print_pdf", PDF_METADATA)
    _wrap(backend_svg.FigureCanvasSVG, "print_svg", SVG_METADATA)


install()
