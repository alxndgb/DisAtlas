"""House figure style shared by every main-text figure (fonts, sizes, panel letters)."""
import os
import matplotlib as mpl

FONT_FAMILY = "DejaVu Sans"
BASE_SIZE = 11.5
TITLE_SIZE = 12
LABEL_SIZE = 12
TICK_SIZE = 11
LEGEND_SIZE = 10.5
ANNOT_SIZE = 9.5                     # small counts / n annotations, colour 0.35
PANEL_LABEL_SIZE = 15
PANEL_LABEL_XY = (-0.10, 1.03)       # axes fraction, va="bottom"
DPI = 300

RC = {
    "font.family": FONT_FAMILY,
    "font.size": BASE_SIZE,
    "axes.titlesize": TITLE_SIZE,
    "axes.titlelocation": "left",
    "axes.labelsize": LABEL_SIZE,
    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "legend.fontsize": LEGEND_SIZE,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": DPI,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
}

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DATA = os.path.join(HERE, "data")


def apply_style():
    """Set the house rcParams (call once, before creating figures)."""
    mpl.rcParams.update(RC)


def panel_label(ax, letter, xy=PANEL_LABEL_XY, **kw):
    """Bold lower-case panel letter at the house position, top-left of ax."""
    return ax.text(xy[0], xy[1], letter.lower(), transform=ax.transAxes,
                   fontweight="bold", fontsize=PANEL_LABEL_SIZE, va="bottom",
                   ha="left", **kw)


def save(fig, stem, **kw):
    """Write <stem>.png (300 dpi) and <stem>.pdf next to the calling script (no date metadata)."""
    out = []
    meta = {"png": {"Software": None}, "pdf": {"CreationDate": None, "Creator": None, "Producer": None}}
    for ext in ("png", "pdf"):
        p = os.path.join(HERE, f"{stem}.{ext}")
        fig.savefig(p, dpi=DPI, metadata=meta[ext], **kw)
        out.append(p)
    return out
