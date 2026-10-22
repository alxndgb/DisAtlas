"""House figure style for the Supplementary Information figures (same settings as figures/main/disatlas_style.py)."""
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
SOURCE_DATA = os.path.join(HERE, "source_data")


def src(name):
    """Path of a shipped source-data table for the SI figures."""
    return os.path.join(SOURCE_DATA, name)


def apply_style():
    """Set the house rcParams (call once, before creating figures)."""
    mpl.rcParams.update(RC)


def panel_label(ax, letter, xy=PANEL_LABEL_XY, **kw):
    """Bold lower-case panel letter at the house position, top-left of ax."""
    return ax.text(xy[0], xy[1], letter.lower(), transform=ax.transAxes,
                   fontweight="bold", fontsize=PANEL_LABEL_SIZE, va="bottom",
                   ha="left", **kw)


# Compact style used by the appendix figures drawn at 7-8 in width (MechDist schematic, IxIDN-Matched funnel).
COMPACT_RC = {
    "font.family": "sans-serif", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6, "xtick.direction": "out",
    "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": True,
    "axes.spines.bottom": True, "axes.grid": False, "legend.frameon": False, "figure.dpi": 200,
    "savefig.dpi": 300, "savefig.bbox": "tight", "axes.titleweight": "normal", "axes.titlelocation": "left",
    "axes.labelweight": "normal", "lines.linewidth": 1.2, "patch.linewidth": 0.6, "pdf.fonttype": 42,
    "ps.fonttype": 42}


def apply_compact_style():
    mpl.rcParams.update(COMPACT_RC)


def compact_panel_letter(ax, letter):
    """Bold lower-case letter outside the top-left corner (compact style)."""
    ax.text(-0.18, 1.02, letter.lower(), transform=ax.transAxes, fontweight="bold",
            fontsize=mpl.rcParams["font.size"] + 1, va="bottom", ha="left")


def muted_focal_pair(focal_color):
    """[focal_color, muted second colour of the default cycle]: the focal-vs-other pair of the compact figures."""
    r, g, b = mpl.colors.to_rgb(mpl.rcParams["axes.prop_cycle"].by_key()["color"][1])
    m = (r + g + b) / 3
    return [focal_color, mpl.colors.to_hex((0.3 * r + 0.7 * m, 0.3 * g + 0.7 * m, 0.3 * b + 0.7 * m))]


def save(fig, stem, **kw):
    """Write <stem>.png (300 dpi) and <stem>.pdf next to the calling script (no date metadata)."""
    out = []
    dpi = kw.pop("dpi", DPI)
    meta = {"png": {"Software": None}, "pdf": {"CreationDate": None, "Creator": None, "Producer": None}}
    for ext in ("png", "pdf"):
        p = os.path.join(HERE, f"{stem}.{ext}")
        fig.savefig(p, dpi=dpi, metadata=meta[ext], **kw)
        out.append(p)
    return out
