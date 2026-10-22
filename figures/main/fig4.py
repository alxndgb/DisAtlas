"""Figure 4 (Methods, UPNA-DDAL): two-stage architecture schematic, house style.

Model content is unchanged: 32-d DisAtlas embedding (MechDist + Node2Vec, 5 seeds),
Hadamard product -> 128-unit ReLU -> dropout 0.3 -> sigmoid, 4,353 parameters,
5-seed ensemble mean. Panel b is illustrative and not computed: only the aHUS-PNH edge
corresponds to a reported result (Results case study); the other nodes are placeholders.

Designed at 7.2 in width so that house font sizes stay legible when the figure is
printed at single-column width (~3.5 in; smallest text ~5 pt).

    python figures/main/fig4.py      -> figures/main/fig4_upna_ddal_architecture.{png,pdf}
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle, Circle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from disatlas_style import apply_style, panel_label, save, TITLE_SIZE, ANNOT_SIZE, LEGEND_SIZE, DPI

apply_style()
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["mathtext.fontset"] = "dejavusans"

EDGE = "#6b1f1f"
GREY = "#c8c8c8"
GREEN = "#e3f1d6"
PANEL = "#ebebeb"
NEWFILL = "#fff6e0"
POSCOL = "#e8a33d"
NEGCOL = "#c0504d"
NODECOL = "#8f9fd4"
FS = 10.5          # box text; smallest text in figure is ANNOT_SIZE-ish 10


def group(ax, x, y, w, h, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=2",
                                fc=fc, ec=EDGE, lw=0.9))


def para(ax, cx, cy, w, h, text, fc="white", fs=FS, skew=0.9):
    pts = [(cx - w / 2 - skew, cy - h / 2), (cx + w / 2 - skew, cy - h / 2),
           (cx + w / 2 + skew, cy + h / 2), (cx - w / 2 + skew, cy + h / 2)]
    ax.add_patch(Polygon(pts, closed=True, fc=fc, ec=EDGE, lw=0.8))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs, linespacing=1.25)


def box(ax, cx, cy, w, h, text, fc="white", fs=FS):
    ax.add_patch(Rectangle((cx - w / 2, cy - h / 2), w, h, fc=fc, ec=EDGE, lw=0.8))
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs, linespacing=1.3)


def arrow(ax, p0, p1, lw=1.0):
    ax.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="-|>", color="black", lw=lw,
                mutation_scale=9, shrinkA=0, shrinkB=0))


def seg(ax, pts, lw=1.0):
    xs, ys = zip(*pts); ax.plot(xs, ys, color="black", lw=lw, solid_capstyle="butt")


def build():
    W, HA, HB, GAP = 7.2, 3.55, 2.15, 0.3
    H = HA + HB + GAP + 0.35
    fig = plt.figure(figsize=(W, H))
    # data units: 100 across the width, same scale vertically
    sx = 100 / W
    axa = fig.add_axes([0.0, (HB + GAP) / H, 1.0, HA / H])
    axa.set_xlim(0, 100); axa.set_ylim(0, HA * sx); axa.axis("off")
    axb = fig.add_axes([0.0, 0.0, 1.0, HB / H])
    axb.set_xlim(0, 100); axb.set_ylim(0, HB * sx); axb.axis("off")

    for ax, letter, title in [(axa, "a", "Two-stage architecture"),
                              (axb, "b", "Illustrative output: pair probabilities assembled into a disease network")]:
        panel_label(ax, letter, xy=(0.0, 1.0))
        ax.text(0.04, 1.0, title, transform=ax.transAxes, fontsize=TITLE_SIZE, va="bottom", ha="left")

    # ---- panel a -------------------------------------------------------
    yA, yB, yM = 34, 14, 24
    para(axa, 8.6, yA, 14.6, 8.5, "Disease A\n(e.g., aHUS)", fc=GREY)
    para(axa, 8.6, yB, 14.6, 8.5, "Disease B\n(e.g., PNH)", fc=GREY)
    arrow(axa, (16.5, yA), (19, yA)); arrow(axa, (16.5, yB), (19, yB))

    group(axa, 19, 4, 27, 40, GREEN)
    axa.text(32.5, 46, "Stage 1: unsupervised", ha="center", va="bottom", fontsize=FS)
    box(axa, 32.5, yA, 25.4, 13, "DisAtlas embedding\n(MechDist + Node2Vec,\n32-d, 5 seeds)", fc=NEWFILL, fs=9.8)
    box(axa, 32.5, yB, 25.4, 13, "DisAtlas embedding\n(MechDist + Node2Vec,\n32-d, 5 seeds)", fc=NEWFILL, fs=9.8)
    seg(axa, [(45.2, yA), (48, yA), (48, yB), (45.2, yB)]); arrow(axa, (48, yM), (50.6, yM))

    group(axa, 49.5, 4, 29, 40, PANEL)
    axa.text(64, 46, "Stage 2: supervised", ha="center", va="bottom", fontsize=FS)
    box(axa, 64, yM, 26.8, 36,
        "UPNA-DDAL decoder\n\n"
        r"Hadamard product $e_A \odot e_B$" "\n"
        r"$\downarrow$" "\n128-unit ReLU layer\n" r"$\downarrow$" "\ndropout (0.3)\n"
        r"$\downarrow$" "\nsigmoid output\n\n4,353 parameters;\n5-seed ensemble mean",
        fc=NEWFILL, fs=10)
    arrow(axa, (77.4, yM), (80.6, yM))

    para(axa, 90.0, yM, 17.4, 13, "Association\nprobability\n" r"$p(A,B)\in[0,1]$", fc=GREY, fs=10)

    # ---- panel b -------------------------------------------------------
    group(axb, 1, 1, 98, 26.5, PANEL)
    nodes = {"aHUS": (12, 10), "PNH": (26, 16), "Disease C": (40, 7),
             "Disease D": (53, 15), "Disease E": (66, 6)}
    for a, b in [("PNH", "Disease C"), ("Disease D", "Disease E")]:
        (x1, y1), (x2, y2) = nodes[a], nodes[b]; axb.plot([x1, x2], [y1, y2], color=NEGCOL, lw=2.4, zorder=1)
    (x1, y1), (x2, y2) = nodes["aHUS"], nodes["PNH"]; axb.plot([x1, x2], [y1, y2], color=POSCOL, lw=2.4, zorder=1)
    for name, (x, y) in nodes.items():
        axb.add_patch(Circle((x, y), 2.4, fc=NODECOL, ec=EDGE, lw=0.8, zorder=2))
        axb.text(x, y + 3.3, name, ha="center", va="bottom", fontsize=10, zorder=3)
    axb.plot([74, 79], [18, 18], color=POSCOL, lw=2.4)
    axb.text(80.5, 18, "Positive pairs\n(IxIDN+)", fontsize=LEGEND_SIZE, va="center", linespacing=1.15)
    axb.plot([74, 79], [9, 9], color=NEGCOL, lw=2.4)
    axb.text(80.5, 9, "Negative pairs\n(IxIDN\u2212)", fontsize=LEGEND_SIZE, va="center", linespacing=1.15)
    return fig


if __name__ == "__main__":
    fig = build()
    save(fig, "fig4_upna_ddal_architecture", bbox_inches="tight", pad_inches=0.04)
