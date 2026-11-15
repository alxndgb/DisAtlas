"""House figure style for the DisAtlas manuscript.

Reference: the main-text disease clustering map.
Every main-text figure uses these settings so that fonts, sizes and panel
labels match across the paper.

    from disatlas_style import apply_style, panel_label
    apply_style()
    fig, ax = plt.subplots(...)
    panel_label(ax, "a")
"""
import matplotlib as mpl

FONT_FAMILY = "sans-serif"          # resolves to DejaVu Sans on the cluster
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
}


def apply_style():
    """Set the house rcParams (call once, before creating figures)."""
    mpl.rcParams.update(RC)


def panel_label(ax, letter, xy=PANEL_LABEL_XY, **kw):
    """Bold lower-case panel letter at the house position, top-left of ax."""
    return ax.text(xy[0], xy[1], letter.lower(), transform=ax.transAxes,
                   fontweight="bold", fontsize=PANEL_LABEL_SIZE, va="bottom",
                   ha="left", **kw)


def savefig(fig, path):
    fig.savefig(path, dpi=DPI, bbox_inches="tight")


# ---------------------------------------------------------------------------
# Raster figures without plotting source (pasted composites, PowerPoint
# exports): mask an old label and draw the house label at the same scale.
# Pixel size of the label = PANEL_LABEL_SIZE pt at the figure's print scale:
# px = 15/72 * DPI * (print_width_in / figure_width_in_at_DPI) -> pass
# label_px explicitly, measured from the reference figure.
# ---------------------------------------------------------------------------
def raster_relabel(png_in, png_out, labels, label_px, mask_boxes=(), bg=(255, 255, 255)):
    """labels: [(letter, x_px, y_px_top)]; mask_boxes: [(x0, y0, x1, y1)]."""
    from PIL import Image, ImageDraw, ImageFont
    from matplotlib import font_manager
    im = Image.open(png_in).convert("RGB")
    d = ImageDraw.Draw(im)
    for b in mask_boxes:
        d.rectangle(b, fill=bg)
    fpath = font_manager.findfont(font_manager.FontProperties(family="DejaVu Sans", weight="bold"))
    font = ImageFont.truetype(fpath, int(label_px))
    for letter, x, y in labels:
        d.text((x, y), letter.lower(), fill=(0, 0, 0), font=font)
    im.save(png_out, dpi=(DPI, DPI))
    return png_out
