"""Main Figure 2 (three panels): (a) graph and embedding construction, (b) retrospective IxIDN-Matched benchmark
feeding Task 1 and Task 2 (embeddings from (a) drawn into both tasks), (c) simulated prospective event-time test.
Box text for (c) follows Methods "Event-time drug-centred evaluation" and the Fig 7 legend.
Pure schematic (no data inputs). Usage: python figures/main/fig2_disatlas_pipeline.py [OUTDIR]  (default: next to this script)
"""
import sys as _s, pathlib as _p; _s.path.insert(0, str(_p.Path(__file__).resolve().parents[2] / "src")); import disatlas.deterministic  # noqa: E402,F401  reproducible PDF/SVG bytes
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

GREY_F, GREY_E = "#EBEBEB", "#6E6E6E"
YEL_F, YEL_E = "#FFF2C8", "#C9960C"
BLU_F, BLU_E = "#E1EAF5", "#3A62A8"
WHT_F, RED_E = "#FFFFFF", "#7A2A2A"
PINK_F = "#F6E3E3"
PANEL_F, PANEL_E = "#F4F4F4", "#9A9A9A"
MUTED = "#555555"; LINE = "#333333"

SIDE_F, SIDE_E, SIDE_T = "#F7F7F7", "#3C3C3C", "#1A1A1A"   # dark text/borders; solid lines (panel title marks the simulation)
BODY_T = "#262626"; NOTE_T = "#333333"

fig, ax = plt.subplots(figsize=(11.8, 12.3))
ax.set_xlim(0, 102); ax.set_ylim(22.5, 136); ax.axis("off")

def rbox(x, y, w, h, fc, ec, title, body="", tfs=9.6, bfs=8.2, gap=4.1, ls="-", tc="black", bc=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.9",
                                fc=fc, ec=ec, lw=1.3, ls=ls, zorder=2))
    if body:
        ax.text(x + w / 2, y + h - 1.3, title, ha="center", va="top", fontsize=tfs, weight="bold",
                linespacing=1.2, zorder=3, color=tc)
        ax.text(x + w / 2, y + h - gap, body, ha="center", va="top", fontsize=bfs, linespacing=1.3, zorder=3, color=bc or tc)
    else:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center", fontsize=tfs, weight="bold",
                linespacing=1.25, zorder=3, color=tc)
    return dict(l=(x, y + h / 2), r=(x + w, y + h / 2), t=(x + w / 2, y + h), b=(x + w / 2, y))

def arr(p, q, side=False):
    ax.annotate("", xy=q, xytext=p, zorder=1,
                arrowprops=dict(arrowstyle="-|>", color=SIDE_E if side else LINE, lw=1.1, mutation_scale=10,
                                shrinkA=0, shrinkB=0, ls="-"))

def line(xs, ys):
    ax.plot(xs, ys, color=LINE, lw=1.1, zorder=1)

def panel_label(x, y, letter, text):
    ax.text(x, y, letter, fontsize=11, weight="bold", va="center")
    ax.text(x + 2.6, y, text, fontsize=10.5, va="center")

# ================= Panel a: from evidence to embeddings (offset O) =================
O = 56
panel_label(0.3, 78.5 + O, "a", "Biological evidence is combined into a mechanism-aware disease graph and embedded")

s  = rbox(0.5, 64 + O, 14, 7.5, GREY_F, GREY_E, "Disease\nsymptoms")
g  = rbox(0.5, 53 + O, 14, 7.5, GREY_F, GREY_E, "Disease genes")
lt = rbox(0.5, 39.5 + O, 14, 11, GREY_F, GREY_E, "Literature", "disease\u2013gene\nstatements, dated\nby publication year")

br = rbox(20, 60 + O, 24, 12.5, GREY_F, GREY_E, "Broad similarity",
          "four views: gene function (GO),\ncell-type expression (HPA),\nprotein domains (Pfam) and\nphenotype (HPO); mean rank",
          bfs=7.6, gap=3.6)
sp = rbox(20, 44.5 + O, 24, 12.5, YEL_F, YEL_E, "Specific links",
          "disease\u2013gene and entity\nlinks from curated sources\nand the literature")
mn = rbox(49, 60 + O, 21, 12.5, YEL_F, YEL_E, "Mechanism nodes",
          "shared pathways and\ndruggable target genes,\ne.g. complement, TNF,\nJAK, IL-1, FcRn")
mb = rbox(49, 44.5 + O, 21, 12.5, YEL_F, YEL_E, "Mechanism balancing",
          "sparsely studied mechanisms\nup-weighted toward the\nbest-curated one (\u2264 10\u00d7)")
gr = rbox(76, 62 + O, 25.5, 10.5, BLU_F, BLU_E, "Disease graph\n(full, all years)",
          "all evidence pooled", tfs=9.3, gap=6.4)
nv = rbox(76, 43 + O, 25.5, 15, BLU_F, BLU_E, "Weighted node2vec",
          "random walks follow edge\nweights; weights from evidence\nstrength and mechanism\nbalancing \u2014 no time term")
em = rbox(76, 29 + O, 25.5, 10, BLU_F, BLU_E, "DisAtlas embeddings",
          "one vector per disease;\nmean of 5 seeds")

arr(s["r"], (20, 67.75 + O))
arr((14.5, 58.5 + O), (20, 62.5 + O)); arr((14.5, 55 + O), (20, 52 + O))
arr(lt["r"], (20, 48 + O))
ax.text(17.2, 68.6 + O, "HPO", ha="center", va="bottom", fontsize=7.5, color=MUTED, style="italic")
ax.text(17.2, 61.0 + O, "GO, HPA,\nPfam", ha="center", va="bottom", fontsize=6.6, color=MUTED, style="italic", linespacing=1.05)
arr(sp["r"], mb["l"]); arr(mn["b"], mb["t"])
line([44, 46.5, 46.5], [66.25 + O, 66.25 + O, 75 + O]); line([46.5, 73], [75 + O, 75 + O]); line([73, 73], [75 + O, 69.5 + O])
arr((73, 69.5 + O), (76, 69.5 + O))
line([70, 73, 73], [50.75 + O, 50.75 + O, 65 + O]); arr((73, 65 + O), (76, 65 + O))
arr(gr["b"], nv["t"]); arr(nv["b"], em["t"])
ax.text(32, 41.5 + O, "no drug or clinical-trial information enters the graph", ha="center", va="top",
        fontsize=8, color=MUTED, style="italic")

# ================= Panel b: retrospective benchmark =================
panel_label(0.3, 82.5, "b", "Retrospective benchmark built from clinical-trial outcomes")
DY = 9.5
def subpanel(y0, h, title):
    ax.add_patch(FancyBboxPatch((0.5, y0), 69, h, boxstyle="round,pad=0,rounding_size=1.4",
                                fc=PANEL_F, ec=PANEL_E, lw=1.2, zorder=0))
    ax.text(2, y0 + h - 0.8, title, ha="left", va="top", fontsize=9, weight="bold", color="#444444")
subpanel(48.5 + DY, 22, "Trial co-investigation (IxIDN)")
tr = rbox(2, 54 + DY, 9.5, 9.5, GREY_F, GREY_E, "Drug\u2013\ndisease\ntrials", tfs=8.8)
ix = rbox(13.5, 53 + DY, 11.5, 11.5, YEL_F, YEL_E, "IxIDN", "trial graph:\ndiseases linked\nby shared drugs", bfs=7.8)
po = rbox(27.5, 59.5 + DY, 20, 8.5, WHT_F, RED_E, "Positive pairs (IxIDN\u207a)",
          "trialled with the same drug;\ndated by first joint trial", tfs=8.8, bfs=7.8, gap=3.7)
ne = rbox(27.5, 50 + DY, 20, 8.5, PINK_F, RED_E, "Negative pairs (IxIDN\u207b)",
          "\u2265 3 links apart, no shared-\nmechanism drug; matched", tfs=8.8, bfs=7.8, gap=3.7)
bm = rbox(50.5, 50 + DY, 17.5, 18, WHT_F, RED_E, "IxIDN-Matched\nbenchmark",
          "disease-disjoint\ntune / test split\n\nmetric: AUROC\n\ntime-censored",
          tfs=8.8, bfs=7.8, gap=6.4)
arr(tr["r"], (13.5, 58.75 + DY))
line([25, 26.3, 26.3], [58.75 + DY, 58.75 + DY, 63.75 + DY]); arr((26.3, 63.75 + DY), (27.5, 63.75 + DY))
line([26.3, 26.3], [58.75 + DY, 54.25 + DY]); arr((26.3, 54.25 + DY), (27.5, 54.25 + DY))
line([47.5, 49.2], [63.75 + DY, 63.75 + DY]); line([47.5, 49.2], [54.25 + DY, 54.25 + DY]); line([49.2, 49.2], [54.25 + DY, 63.75 + DY])
arr((49.2, 59 + DY), (50.5, 59 + DY))
t1 = rbox(74, 52.5 + DY, 23, 11, GREY_F, GREY_E, "Task 1\nDisease\u2013disease\nassociation", tfs=9)
t2 = rbox(74, 35.5 + DY, 23, 10, GREY_F, GREY_E, "Task 2\nDisease clustering", tfs=9)
line([68, 71.5], [59 + DY, 59 + DY]); arr((71.5, 59 + DY), (74, 59 + DY))
line([71.5, 71.5], [59 + DY, 41 + DY]); arr((71.5, 41 + DY), (74, 41 + DY))
ax.text(45, 46.6 + DY, "benchmark pairs never enter the graph (checked automatically)", ha="center", va="center",
        fontsize=8.2, color=MUTED, style="italic")

# full-graph embeddings from (a) feed both tasks (blue = embedding flow); routed along the right edge
def barr(p, q):
    ax.annotate("", xy=q, xytext=p, zorder=1, arrowprops=dict(arrowstyle="-|>", color=BLU_E, lw=1.3, mutation_scale=11, shrinkA=0, shrinkB=0))
# one trunk from the embeddings box down the right margin, forking into the right side of Task 1 and Task 2
XT = 99.6
ax.plot([XT, XT], [29 + O, t2["r"][1]], color=BLU_E, lw=1.3, zorder=1, solid_capstyle="butt")
ax.plot([XT], [t1["r"][1]], marker="o", ms=3.2, color=BLU_E, zorder=1)
barr((XT, t1["r"][1]), t1["r"])
barr((XT, t2["r"][1]), t2["r"])
# ================= Panel c: prospective (event-time) arm =================
panel_label(0.3, 41.5, "c", "Simulated prospective test, year by year")
yb, hb = 25.0, 13.0
kw = dict(tfs=8.8, bfs=7.5, ls="-", tc=SIDE_T, bc=BODY_T)
d1 = rbox(0.5, yb, 15, hb, SIDE_F, SIDE_E, "Dated edges", "each graph edge in (a)\ndated by its first\nevidence", gap=3.9, **kw)
d2 = rbox(18.5, yb, 16, hb, SIDE_F, SIDE_E, "Yearly graph\nsnapshots", "evidence up to the\nend of year y\u22121", gap=6.6, **kw)
d3 = rbox(37.5, yb, 16, hb, SIDE_F, SIDE_E, "Snapshot\nembeddings", "DisAtlas retrained\non each snapshot", gap=6.6, **kw)
d4 = rbox(56.5, yb, 20, hb, SIDE_F, SIDE_E, "Score each drug",
          "candidates: diseases not yet\ntrialled with the drug; score =\nhighest similarity to the\ndrug\u2019s earlier trial diseases", gap=3.9, **kw)
d5 = rbox(79.5, yb, 22, hb, SIDE_F, SIDE_E, "Rank its new trials",
          "percentile of each disease the\ndrug first enters in year y,\namong diseases of similar\ntrial activity (and degree),\nvs weighted Adamic\u2013Adar", gap=3.4, bfs=7.2, tfs=8.8, ls="-", tc=SIDE_T, bc=BODY_T)
for p, q in [(d1, d2), (d2, d3), (d3, d4), (d4, d5)]:
    arr(p["r"], q["l"], side=True)
# trial table behind IxIDN (panel b) supplies the dated trial starts used to score panel c; dark arrow, lands on d4 only
XI = ix["b"][0]; YR = 47.0
line([XI, XI], [ix["b"][1], YR]); line([XI, d4["t"][0]], [YR, YR]); arr((d4["t"][0], YR), d4["t"])
ax.text((XI + d4["t"][0]) / 2, YR + 0.5, "dated trial starts (scoring only)", ha="center", va="bottom", fontsize=7.8, color=NOTE_T, style="italic")

import os, sys
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(OUT, "disatlas_pipeline.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(OUT, "disatlas_pipeline.pdf"), bbox_inches="tight", facecolor="white")
