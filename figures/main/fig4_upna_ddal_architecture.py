"""Main Figure 4: UPNA-DDAL two-stage architecture schematic (house style).

Panel a names the stages only (dimensions, seeds, Hadamard product, layer size, dropout and parameter count are in
Methods and the SI). Panel b shows six control pairs with their DisAtlas link probabilities (three known related, three
known unrelated), read from data/derived/control_pairs/controls_percentile_rank.tsv (scores of the released five-seed
pair-score matrix; full 17-pair list in the SI) and asserted against the Figure 4 caption values.

Usage: python figures/main/fig4_upna_ddal_architecture.py [OUTDIR] [STEM]
Outputs: UPNA_DDAL_architecture.{png,pdf,svg} next to this script by default.
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle, Circle
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import csv
CTRL = {r["tag"]: float(r["score"]) for r in csv.DictReader(open(os.environ.get("FIG4_CONTROLS", os.path.join(HERE, "..", "..", "data", "derived", "control_pairs", "controls_percentile_rank.tsv"))), delimiter="\t")}
CAPTION = {"aHUS-PNH": "0.955", "Crohn-UC": "0.975", "Psoriasis-PsA": "0.954", "DI-T2D": "0.079", "PKD-PCOS": "0.163", "ALS-SSc": "0.145"}
for t, v in CAPTION.items():
    assert f"{CTRL[t]:.3f}" == v, (t, CTRL[t], v)
from disatlas.figstyle import apply_style, panel_label, TITLE_SIZE, ANNOT_SIZE, LEGEND_SIZE, DPI

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
    W, HA, HB, GAP = 7.2, 3.55, 2.7, 0.3
    H = HA + HB + GAP + 0.35
    fig = plt.figure(figsize=(W, H))
    # data units: 100 across the width, same scale vertically
    sx = 100 / W
    axa = fig.add_axes([0.0, (HB + GAP) / H, 1.0, HA / H])
    axa.set_xlim(0, 100); axa.set_ylim(0, HA * sx); axa.axis("off")
    axb = fig.add_axes([0.0, 0.0, 1.0, HB / H])
    axb.set_xlim(0, 100); axb.set_ylim(0, HB * sx); axb.axis("off")

    for ax, letter, title in [(axa, "a", "Two-stage architecture"),
                              (axb, "b", "Predicted association probabilities for six known pairs")]:
        panel_label(ax, letter, xy=(0.0, 1.0))
        ax.text(0.04, 1.0, title, transform=ax.transAxes, fontsize=TITLE_SIZE, va="bottom", ha="left")

    # ---- panel a -------------------------------------------------------
    yA, yB, yM = 34, 14, 24
    para(axa, 8.6, yA, 14.6, 8.5, "Disease A\n(e.g., aHUS)", fc=GREY)
    para(axa, 8.6, yB, 14.6, 8.5, "Disease B\n(e.g., PNH)", fc=GREY)
    arrow(axa, (16.5, yA), (19, yA)); arrow(axa, (16.5, yB), (19, yB))

    group(axa, 19, 4, 27, 40, GREEN)
    axa.text(32.5, 46, "Stage 1: unsupervised", ha="center", va="bottom", fontsize=FS)
    box(axa, 32.5, yA, 25.4, 13, "DisAtlas embedding\nof disease A", fc=NEWFILL, fs=FS)
    box(axa, 32.5, yB, 25.4, 13, "DisAtlas embedding\nof disease B", fc=NEWFILL, fs=FS)
    seg(axa, [(45.2, yA), (48, yA), (48, yB), (45.2, yB)]); arrow(axa, (48, yM), (50.6, yM))

    group(axa, 49.5, 4, 29, 40, PANEL)
    axa.text(64, 46, "Stage 2: supervised", ha="center", va="bottom", fontsize=FS)
    box(axa, 64, yM, 26.8, 16, "Decoder\n(combines the\ntwo embeddings)", fc=NEWFILL, fs=FS)
    arrow(axa, (77.4, yM), (80.6, yM))

    para(axa, 90.0, yM, 17.4, 13, "Association\nprobability\n" r"$p(A,B)\in[0,1]$", fc=GREY, fs=10)

    # ---- panel b -------------------------------------------------------
    group(axb, 1, 1, 98, 29.5, PANEL)
    pairs = [("aHUS", "PNH", CTRL["aHUS-PNH"], True, False),
             ("Crohn disease", "Ulcerative colitis", CTRL["Crohn-UC"], True, False),
             ("Psoriasis", "Psoriatic arthritis", CTRL["Psoriasis-PsA"], True, False),
             ("Diabetes insipidus", "Type 2 diabetes", CTRL["DI-T2D"], False, False),
             ("Polycystic kidney\ndisease", "Polycystic ovary\nsyndrome", CTRL["PKD-PCOS"], False, False),
             ("Amyotrophic lateral\nsclerosis", "Systemic sclerosis", CTRL["ALS-SSc"], False, False)]
    xs = [9, 25, 41, 57, 73, 89]
    for x, (a, b, sc, rel, miss) in zip(xs, pairs):
        ytop, ybot = 20.5, 10.0
        col = POSCOL if rel else NEGCOL
        axb.plot([x, x], [ytop, ybot], color=col, lw=2.4, ls="-" if rel else (0, (2.2, 1.6)), zorder=1)
        for yy, nm, va, dy in [(ytop, a, "bottom", 2.6), (ybot, b, "top", -2.6)]:
            axb.add_patch(Circle((x, yy), 1.7, fc=NODECOL, ec=EDGE, lw=0.8, zorder=2))
            axb.text(x, yy + dy, nm, ha="center", va=va, fontsize=8.6, linespacing=1.05, zorder=3)
        lab = f"{sc:.3f}" + ("\nmiss" if miss else "")
        axb.text(x + 1.3, (ytop + ybot) / 2, lab, ha="left", va="center", fontsize=8.6,
                 color="black", fontweight="bold" if miss else "normal", linespacing=1.1)
    axb.text(25, 32.0, "Known related (solid)", ha="center", fontsize=LEGEND_SIZE, color=POSCOL, va="bottom")
    axb.text(73, 32.0, "Known unrelated (dashed)", ha="center", fontsize=LEGEND_SIZE, color=NEGCOL, va="bottom")
    return fig


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else HERE
    stem = sys.argv[2] if len(sys.argv) > 2 else "UPNA_DDAL_architecture"
    fig = build()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(out, f"{stem}.{ext}"), dpi=max(DPI, 300),
                    bbox_inches="tight", pad_inches=0.04)
    print("wrote", stem, "{png,svg,pdf} ->", out)
