"""SI figure: approved-repurposing benchmark schematic (moved out of main Figure 2), same style.
    python figS_approved_repurposing.py OUTDIR
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    apply_figure_style()  # figure-style skill helper when run in the analysis kernel
except NameError:
    pass
from matplotlib.patches import FancyBboxPatch
import pandas as pd
# benchmark composition from the shipped pair/split table (data/derived/si_fda/README.md)
_A = pd.read_csv(DERIVED / "si_fda" / "benchmark_pairs_split.tsv", sep="\t")
NP, ND = len(_A), len(set(_A.mondo_a) | set(_A.mondo_b))
NH, NE = int((_A.stratum == "hard").sum()), int((_A.stratum == "easy").sum())
NDEV, NHO = int((_A.split == "dev").sum()), int((_A.split == "holdout").sum())
NCH = int(((_A.stratum == "hard") & (_A.split == "complement_six")).sum())

GREY_F, GREY_E = "#EBEBEB", "#6E6E6E"
YEL_F, YEL_E = "#FFF2C8", "#C9960C"
BLU_F, BLU_E = "#E1EAF5", "#3A62A8"
WHT_F, RED_E = "#FFFFFF", "#7A2A2A"
PINK_F = "#F6E3E3"
PANEL_F, PANEL_E = "#F4F4F4", "#9A9A9A"
MUTED = "#555555"; LINE = "#333333"

SIDE_F, SIDE_E, SIDE_T = "#F7F7F7", "#8C8C8C", "#4A4A4A"

fig, ax = plt.subplots(figsize=(11.8, 3.9))
ax.set_xlim(0, 102); ax.set_ylim(22.5, 47.8); ax.axis("off")
def rbox(x, y, w, h, fc, ec, title, body="", tfs=9.6, bfs=8.2, gap=4.1, ls="-", tc="black"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.9",
                                fc=fc, ec=ec, lw=1.3, ls=ls, zorder=2))
    if body:
        ax.text(x + w / 2, y + h - 1.3, title, ha="center", va="top", fontsize=tfs, weight="bold",
                linespacing=1.2, zorder=3, color=tc)
        ax.text(x + w / 2, y + h - gap, body, ha="center", va="top", fontsize=bfs, linespacing=1.3, zorder=3, color=tc)
    else:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center", fontsize=tfs, weight="bold",
                linespacing=1.25, zorder=3, color=tc)
    return dict(l=(x, y + h / 2), r=(x + w, y + h / 2), t=(x + w / 2, y + h), b=(x + w / 2, y))

def arr(p, q, side=False):
    ax.annotate("", xy=q, xytext=p, zorder=1,
                arrowprops=dict(arrowstyle="-|>", color=SIDE_E if side else LINE, lw=1.1, mutation_scale=10,
                                shrinkA=0, shrinkB=0, ls=(0, (3, 2)) if side else "-"))

def line(xs, ys):
    ax.plot(xs, ys, color=LINE, lw=1.1, zorder=1)

def panel_label(x, y, letter, text):
    ax.text(x, y, letter, fontsize=11, weight="bold", va="center")
    ax.text(x + 2.6, y, text, fontsize=10.5, va="center")


def subpanel(y0, h, title):
    ax.add_patch(FancyBboxPatch((0.5, y0), 69, h, boxstyle="round,pad=0,rounding_size=1.4",
                                fc=PANEL_F, ec=PANEL_E, lw=1.2, zorder=0))
    ax.text(2, y0 + h - 0.8, title, ha="left", va="top", fontsize=9, weight="bold", color="#444444")

subpanel(24.5, 22.5, "Approved repurposing (launched drug-sharing pairs)")
src = rbox(2, 28.5, 15, 13.5, GREY_F, GREY_E, "Shared\nlaunched drug",
           "IxIDN\u207a pairs whose\nshared drug is\nlaunched for both\ndiseases", tfs=8.8, bfs=7.6, gap=5.3)
pc = rbox(19.5, 29.5, 14.5, 11.5, YEL_F, YEL_E, "Positive controls",
          f"{NP:,} pairs over\n{ND:,} diseases\n\nno negative pairs", tfs=8.8, bfs=7.8, gap=3.8)
st = rbox(36.5, 35.5, 16.5, 9, WHT_F, RED_E, "Strata",
          f"hard ({NH:,}): no shared\ngene, HPO or GO link\neasy ({NE:,}): the rest", tfs=8.8, bfs=7.4, gap=3.4)
sl = rbox(36.5, 25.0, 16.5, 9.9, WHT_F, RED_E, "Split (hard pairs)",
          f"dev {NDEV} / holdout {NHO},\ndrug- and mechanism-\ndisjoint; {NCH} named\ncomplement pairs apart", tfs=8.8, bfs=7.2, gap=3.2)
mt = rbox(55.5, 27, 12.5, 16.5, WHT_F, RED_E, "Recovery\nmetric",
          "per-pair\npercentile vs\n20,000 random\ndisease pairs;\nrecovered at\n\u2265 75th", tfs=8.8, bfs=7.4, gap=5.4)
arr(src["r"], pc["l"])
line([34, 35.2, 35.2], [35.25, 35.25, 40]); arr((35.2, 40), (36.5, 40))
line([35.2, 35.2], [35.25, 30]); arr((35.2, 30), (36.5, 30))
line([53, 54.2], [40, 40]); line([53, 54.2], [30, 30]); line([54.2, 54.2], [30, 40])
arr((54.2, 35.25), (55.5, 35.25))


t1 = rbox(74, 29.75, 23, 11, GREY_F, GREY_E, "Task 1\nDisease\u2013disease\nassociation", tfs=9)
line([68, 71], [35.25, 35.25]); arr((71, 35.25), (74, 35.25))
ax.text(85.5, 26.6, "used for Task 1 only;\nTask 2 uses IxIDN-Matched only", ha="center", va="center",
        fontsize=7.6, color=MUTED, style="italic")
ax.text(35, 23.4, "benchmark pairs never enter the graph (checked automatically)", ha="center", va="center",
        fontsize=8.2, color=MUTED, style="italic")
import os, sys
OUT = str(HERE)
fig.savefig(os.path.join(OUT, "SI_approved_repurposing_schematic.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(OUT, "SI_approved_repurposing_schematic.pdf"), bbox_inches="tight", facecolor="white")
