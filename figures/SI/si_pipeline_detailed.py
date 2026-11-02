"""SI figure: detailed pipeline schematic (si_pipeline_detailed).

Layout and design text follow Methods; the counts drawn in the boxes are cross-checked against
data/derived/pipeline_counts/ by tests/test_sibench_parity.py."""
import sys as _s, pathlib as _p; _s.path.insert(0, str(_p.Path(__file__).resolve().parents[2] / "src")); import disatlas.deterministic  # noqa: E402,F401  reproducible PDF/SVG bytes
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    apply_figure_style()  # figure-style skill helper when run in the analysis kernel
except NameError:
    pass
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

GRN_F, GRN_E, GRN_T = "#E6F1E3", "#3E7D3A", "#3E7D3A"


fig, ax = plt.subplots(figsize=(12.6, 20.82))
ax.set_xlim(0, 102); ax.set_ylim(2, 207); ax.axis("off")

def rbox(x, y, w, h, fc, ec, title, body="", tfs=9.2, bfs=7.6, gap=3.9, ls="-", tc="black", bc=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.9",
                                fc=fc, ec=ec, lw=1.3, ls=ls, zorder=2))
    if body:
        ax.text(x + w / 2, y + h - 1.2, title, ha="center", va="top", fontsize=tfs, weight="bold",
                linespacing=1.2, zorder=3, color=tc)
        ax.text(x + w / 2, y + h - gap, body, ha="center", va="top", fontsize=bfs, linespacing=1.3, zorder=3, color=bc or tc)
    else:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center", fontsize=tfs, weight="bold",
                linespacing=1.25, zorder=3, color=tc)
    return dict(l=(x, y + h / 2), r=(x + w, y + h / 2), t=(x + w / 2, y + h), b=(x + w / 2, y), x=x, y=y, w=w, h=h)

def arr(p, q, side=False, color=None):
    ax.annotate("", xy=q, xytext=p, zorder=1,
                arrowprops=dict(arrowstyle="-|>", color=color or (SIDE_E if side else LINE), lw=1.1, mutation_scale=10,
                                shrinkA=0, shrinkB=0, ls="-"))

def line(xs, ys, color=LINE):
    ax.plot(xs, ys, color=color, lw=1.1, zorder=1)

def panel_label(x, y, letter, text):
    ax.text(x, y, letter, fontsize=11, weight="bold", va="center")
    ax.text(x + 2.6, y, text, fontsize=10.5, va="center")

# ================= Panel a =================
panel_label(0.3, 205, "a", "A broad and a mechanism-specific channel are merged into one disease graph and embedded")
ax.text(0.5, 200.5, "Broad channel", fontsize=9.2, weight="bold", color=MUTED, va="center")
gn = rbox(0.5, 180.5, 17, 17, GREY_F, GREY_E, "Disease genes", "OpenTargets\nOrphanet, OMIM\nGenCC\nAZ BIKG literature", tfs=9.2, bfs=7.6, gap=4.2)
sy = rbox(0.5, 165.5, 17, 11.5, GREY_F, GREY_E, "Disease symptoms", "HPO annotations\nOMIM", tfs=8.8, bfs=7.6, gap=4.2)
CH = [("GO gene-function similarity", "96% of pairs", GREY_F, GREY_E, 190),
      ("HPA cell-type expression", "91% of pairs", YEL_F, YEL_E, 182.5),
      ("Pfam protein domains", "94% of pairs", YEL_F, YEL_E, 175),
      ("HPO phenotype", "35% of pairs", GREY_F, GREY_E, 167.5)]
chb = [rbox(22, yy, 27, 6.3, fc, ec, t, b, tfs=8.6, bfs=7.6, gap=3.4) for t, b, fc, ec, yy in CH]
rk = rbox(53, 166, 19, 31.5, YEL_F, YEL_E, "Rank-mean fusion",
          "each pair is ranked\nin every channel;\nranks are averaged;\nmissing channel\nranks last; best-ranked\npairs kept, at least\n3 per disease\n\n331,234 edges\n13,020 diseases\n(after merging)", tfs=9.2, bfs=7.6, gap=4.6)
for b_ in chb[:3]:
    arr((17.5, 189), b_["l"])
arr(sy["r"], chb[3]["l"])
for b_ in chb:
    arr(b_["r"], (53, 181.75))
G = rbox(83.5, 180.5, 18, 17, BLU_F, BLU_E, "MechDist graph",
         "13,021 diseases\n498,155 weighted\nedges after merging\n(23,086 in both\nchannels)", tfs=9.2, bfs=7.6, gap=4.4)
arr(rk["r"], (83.5, 189))
RM = rbox(83.5, 165, 18, 10.5, WHT_F, RED_E, "Benchmark pairs\nremoved", "10,052 pairs \u2192 488,103\ntraining edges", tfs=8.4, bfs=7.2, gap=6.2)
NV = rbox(83.5, 149, 18, 11.5, BLU_F, BLU_E, "node2vec embedding", "weighted, 32 dimensions\n= DisAtlas embeddings", tfs=8.8, bfs=7.4, gap=4.2)
DC = rbox(83.5, 132, 18, 12, BLU_F, BLU_E, "Pair decoder", "MLP on the element-wise\nproduct, 5-seed ensemble;\nscores any disease pair", tfs=8.8, bfs=7.2, gap=3.9)
arr(G["b"], RM["t"]); arr(RM["b"], NV["t"]); arr(NV["b"], DC["t"])

ax.text(0.5, 161.5, "Mechanism-specific channel", fontsize=9.2, weight="bold", color=GRN_T, va="center")
ME = rbox(0.5, 136, 35, 22.5, GRN_F, GRN_E, "Disease\u2013mechanism evidence", "", tfs=9.2)
ax.texts[-1].set_va("top"); ax.texts[-1].set_y(136 + 22.5 - 1.2)
rows = [("BIKG literature", "gene, pathway, process and\nphenotype links (text-mined)"),
        ("OpenTargets", "genetic association, somatic\nmutation, animal model"),
        ("CTD", "marker and mechanism genes"),
        ("MetaBase", "curated pathway maps"),
        ("Autoantibodies", "autoantibody target genes")]
yy = 154.2
for k_, v_ in rows:
    ax.text(1.6, yy, k_, fontsize=7.6, weight="bold", va="top", zorder=3)
    ax.text(13.3, yy, v_, fontsize=7.3, va="top", linespacing=1.15, zorder=3)
    yy -= 3.3 if "\n" in v_ else 2.3
ax.plot([1.6, 34.4], [yy + 0.6, yy + 0.6], color=GRN_E, lw=0.7, ls=(0, (2, 2)), zorder=3)
yy -= 0.4
ax.text(1.6, yy, "Reactome", fontsize=7.6, weight="bold", style="italic", va="top", zorder=3)
ax.text(13.3, yy, "pathway gene sets\n(disease genes by enrichment)", fontsize=7.3, style="italic", va="top", linespacing=1.15, zorder=3)

def bip(x, y, w, h, rows_, indent=()):
    b_ = rbox(x, y, w, h, GRN_F, GRN_E, "Bipartite graph", "", tfs=9.2)
    ax.texts[-1].set_va("top"); ax.texts[-1].set_y(y + h - 1.2)
    ax.text(x + 0.8, y + h - 5.2, "diseases", fontsize=6.8, style="italic", color=MUTED, va="center", zorder=3)
    ax.text(x + w - 0.8, y + h - 5.2, "entities", fontsize=6.8, style="italic", color=MUTED, va="center", ha="right", zorder=3)
    yy_ = y + h - 8.4
    for i, (d_, lab, e_) in enumerate(rows_):
        ax.text(x + 5.0, yy_, d_, fontsize=7.6, ha="right", va="center", zorder=3)
        ax.text(x + (7.2 if i in indent else 5.9), yy_, lab, fontsize=7.6, va="center", zorder=3)
        ax.text(x + w - 0.8, yy_, e_, fontsize=7.6, ha="right", va="center", zorder=3)
        yy_ -= 2.45
    return b_
# icon: diseases (blue) linked to entities (green)
def icon(cx, cy):
    L = [(cx - 2.2, cy + 0.9 - i * 0.9) for i in range(3)]; Rr = [(cx + 2.2, cy + 1.2 - i * 0.8) for i in range(4)]
    for p_ in L:
        for q_ in Rr:
            ax.plot([p_[0], q_[0]], [p_[1], q_[1]], color=GRN_E, lw=0.5, zorder=3)
    for p_ in L: ax.plot(*p_, "o", ms=3.4, color="#3A62A8", zorder=4)
    for q_ in Rr: ax.plot(*q_, "s", ms=3.2, color=GRN_E, zorder=4)
B1 = bip(38.5, 136, 25, 22.5, [("3,677", "Mechanisms", "357"), ("3,839", "Pathways", "1,642"), ("3,757", "Reactome", "1,377"),
                                ("1,245", "literature", "265"), ("2,437", "Biological processes", "2,234")], indent=(2, 3))
icon(51, 153.3)
SE = rbox(66.5, 128, 14.5, 26.5, GRN_F, GRN_E, "Specific edges",
          "diseases scored\nby shared\nentities; link\nweight = source\nfamilies + ln(1 +\ndated documents);\n50 nearest\nper disease\n\n190,007 edges\n5,194 diseases\n(after merging)", tfs=9.0, bfs=7.2, gap=4.2)
PE = rbox(0.5, 121.5, 35, 10.5, GRN_F, GRN_E, "Disease\u2013pathology evidence", "", tfs=9.2)
ax.texts[-1].set_va("top"); ax.texts[-1].set_y(121.5 + 10.5 - 1.2)
ax.text(1.6, 128.0, "BIKG pathology", fontsize=7.6, weight="bold", va="top", zorder=3)
ax.text(13.3, 128.0, "anatomy, cell type and\nphenotype (text-mined)", fontsize=7.3, va="top", linespacing=1.15, zorder=3)
B2 = bip(38.5, 115, 25, 17.5, [("2,448", "Anatomical sites", "2,819"), ("1,406", "Cell types", "464"), ("2,249", "Phenotype terms", "1,088")])
arr(ME["r"], B1["l"]); arr(PE["r"], (38.5, PE["r"][1]))
arr(B1["r"], (66.5, B1["r"][1]))
line([63.5, 73.75], [123.75, 123.75]); arr((73.75, 123.75), SE["b"])
line([81, 82.3, 82.3], [141.25, 141.25, 184]); arr((82.3, 184), (83.5, 184))
ax.text(0.5, 119.4, "BIKG = AstraZeneca Biological Insight Knowledge Graph;\nCTD = Comparative Toxicogenomics Database",
        fontsize=7.3, style="italic", color=MUTED, va="top", linespacing=1.2)
ax.text(50, 112.4, "no drug or clinical-trial information enters the graph", ha="center", va="center",
        fontsize=8, color=MUTED, style="italic")

# ================= Panel b =================
panel_label(0.3, 105, "b", "Retrospective benchmark built from clinical-trial outcomes")
ax.add_patch(FancyBboxPatch((0.5, 54.5), 71, 47, boxstyle="round,pad=0,rounding_size=1.4", fc=PANEL_F, ec=PANEL_E, lw=1.2, zorder=0))
ax.text(2, 100.7, "Trial co-investigation (IxIDN)", ha="left", va="top", fontsize=9, weight="bold", color="#444444")
kb = dict(tfs=8.6, bfs=7.0, gap=3.6)
tr = rbox(2, 58, 16.5, 39, GREY_F, GREY_E, "Drug\u2013disease trials",
          "all trial records:\n419,240 trials,\n38,389 drugs,\n1,593 diseases\n\nremoved in turn:\n38 non-disease\nindications\n2,266 drugs with no\nspecific mechanism\nof action\n52,056 trials of\n> 2 drugs\n\nkept: 299,122 trials,\n32,023 drugs,\n1,437 diseases", **kb)
ix = rbox(20.5, 58, 12, 39, YEL_F, YEL_E, "IxIDN",
          "trial graph:\ndiseases linked\nwhen the same\ndrug was trialled\nfor both\n\n1,379 diseases\n43,112 links", tfs=8.6, bfs=7.2, gap=3.6)
po = rbox(35, 78.5, 19.5, 18.5, WHT_F, RED_E, "Positive pairs (IxIDN\u207a)",
          "linked in IxIDN: 27,415 candidate\npairs in the 985-disease\nevaluation universe; 4,302 with\na matched negative were checked\nby a language model; 2,182 kept\n(1,285 strong, 897 plausible)", **kb)
ne = rbox(35, 58, 19.5, 18.5, PINK_F, RED_E, "Negative pairs (IxIDN\u207b)",
          "\u2265 3 IxIDN steps apart, no shared\ntrial drug: 155,375 candidates;\nmatched to a positive on trial-\ngraph connectivity, organ system,\nclinical area, name and ontology\nsimilarity, gene coverage;\n2,182 matched negatives", **kb)
bm = rbox(57, 58, 13.5, 39, WHT_F, RED_E, "IxIDN-Matched\nbenchmark",
          "4,364 pairs\n(2,182 IxIDN\u207a,\n2,182 IxIDN\u207b)\n938 diseases\n\ndisease-disjoint\nsplit by MONDO-ID\nhash: 461 tune,\n477 test diseases\n\ntest split: 1,094\npairs with both\ndiseases in the\ntest half (491\nIxIDN\u207a, 603 IxIDN\u207b)\n\nmetric: AUROC\ntime-censored", tfs=8.6, bfs=7.0, gap=6.4)
arr(tr["r"], ix["l"])
line([32.5, 33.8], [77.5, 77.5]); line([33.8, 33.8], [67.25, 87.75]); arr((33.8, 87.75), po["l"]); arr((33.8, 67.25), ne["l"])
line([54.5, 55.8], [87.75, 87.75]); line([54.5, 55.8], [67.25, 67.25]); line([55.8, 55.8], [67.25, 87.75]); arr((55.8, 77.5), bm["l"])
t1 = rbox(75.5, 78, 22, 14.5, GREY_F, GREY_E, "Task 1\nDisease\u2013disease\nassociation", "decoder scores pairs;\nAUROC on IxIDN-Matched", tfs=8.8, bfs=7.2, gap=10.4)
t2 = rbox(75.5, 64, 22, 12.5, GREY_F, GREY_E, "Task 2\nDisease clustering", "13,021 diseases,\nK = 55 clusters", tfs=8.8, bfs=7.2, gap=7.4)
line([70.5, 73.2], [77.5, 77.5]); line([73.2, 73.2], [85.25, 70.25]); arr((73.2, 85.25), (75.5, 85.25)); arr((73.2, 70.25), (75.5, 70.25))
ax.text(52, 52.6, "benchmark pairs never enter the graph (checked automatically)", ha="center", va="center",
        fontsize=8.2, color=MUTED, style="italic")
def barr(p, q):
    ax.annotate("", xy=q, xytext=p, zorder=1, arrowprops=dict(arrowstyle="-|>", color=BLU_E, lw=1.3, mutation_scale=11, shrinkA=0, shrinkB=0))
XT = 99.6
ax.plot([DC["b"][0], DC["b"][0], XT], [DC["b"][1], 109.5, 109.5], color=BLU_E, lw=1.3, zorder=1, solid_capstyle="butt")
ax.plot([XT, XT], [109.5, t2["r"][1]], color=BLU_E, lw=1.3, zorder=1, solid_capstyle="butt")
ax.text(DC["b"][0] - 0.8, 119, "DisAtlas\nembeddings", ha="right", va="center", fontsize=8, color=BLU_E, style="italic")
ax.plot([XT], [t1["r"][1]], marker="o", ms=3.2, color=BLU_E, zorder=1)
barr((XT, t1["r"][1]), t1["r"]); barr((XT, t2["r"][1]), t2["r"])

# ================= Panel c =================
panel_label(0.3, 46, "c", "Simulated prospective test, year by year")
kw = dict(tfs=8.6, bfs=7.2, ls="-", tc=SIDE_T, bc=BODY_T)
yb, hb = 24, 15
d1 = rbox(0.5, yb, 16, hb, SIDE_F, SIDE_E, "Dated edges", "each graph edge in (a)\ndated by its first\nevidence", gap=3.6, **kw)
d2 = rbox(19, yb, 17, hb, SIDE_F, SIDE_E, "Yearly graph\nsnapshots", "for each decision year y,\n1996\u20132021: evidence\nup to the end of\nyear y\u22121 only", gap=6.6, **kw)
d3 = rbox(38.5, yb, 16.5, hb, SIDE_F, SIDE_E, "Snapshot\nembeddings", "DisAtlas retrained\non each snapshot", gap=6.6, **kw)
d4 = rbox(57.5, yb, 20, hb, SIDE_F, SIDE_E, "Score each drug",
          "candidates: diseases in that\nyear\u2019s model not yet trialled\nwith the drug (median 2,241\nin the test pool); score =\nmaximum similarity to the\ndrug\u2019s earlier trial diseases", gap=3.6, **kw)
d5 = rbox(80, yb, 21.5, hb, SIDE_F, SIDE_E, "Rank its new trials",
          "event = a drug\u2019s first trial in\na new disease in year y;\npercentile of that disease\namong comparison candidates", gap=3.6, **kw)
for p, q in [(d1, d2), (d2, d3), (d3, d4), (d4, d5)]:
    arr(p["r"], q["l"], side=True)
yr, hr = 4, 14
e1 = rbox(0.5, yr, 31, hr, SIDE_F, SIDE_E, "Events",
          "74,345 new drug\u2013disease trials in\n8,064 drugs; 35,793 trials in\n6,496 drugs among test-pool\ndiseases (not used to tune DisAtlas)", gap=3.6, **kw)
e2 = rbox(34.5, yr, 33, hr, SIDE_F, SIDE_E, "Comparison sets (co-primary)",
          "candidates of the same prior trial\nactivity (popularity), or of the same\nactivity and similar weighted degree\nin that year\u2019s graph; groups need\n\u2265 20 candidates (adjacent groups merged)", gap=3.6, **kw)
e3 = rbox(70.5, yr, 31, hr, SIDE_F, SIDE_E, "Endpoints vs weighted Adamic\u2013Adar",
          "median over drugs of each drug\u2019s\nmedian event percentile; share of\ndrugs with best new trial in top 10%;\npaired bootstrap over drugs, 2,000\nreplicates", gap=3.6, **kw)
line([16, 90.75], [21, 21], color=SIDE_E)
for e in (e1, e2, e3):
    arr((e["x"] + e["w"] / 2, 21), e["t"], side=True)
ax.plot([90.75, 90.75], [yb, 21], color=SIDE_E, lw=1.1, zorder=1)
XI = ix["b"][0]; YR = 49.8
line([XI, XI], [ix["b"][1], YR]); line([XI, d4["t"][0]], [YR, YR]); arr((d4["t"][0], YR), d4["t"])
ax.text(52, YR - 0.6, "dated trial starts (scoring only)", ha="center", va="top", fontsize=7.8, color=NOTE_T, style="italic")

import os, sys
OUT = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(OUT, "si_pipeline_detailed.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(OUT, "si_pipeline_detailed.pdf"), bbox_inches="tight", facecolor="white")
