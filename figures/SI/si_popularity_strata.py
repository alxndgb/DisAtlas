#!/usr/bin/env python
"""SI figure: popularity-controlled drug-centred holdout (pre-registered).

Reads the summary / boundary / decision tables, draws a 3-panel figure, writes
the caption + one Results sentence with every number formatted from the tables,
and writes a numbers-check table (every number used, with its source row).

Usage: python figures/SI/si_popularity_strata.py  (inputs: data/derived/popularity_strata/)
"""
import sys as _s, pathlib as _p; _s.path.insert(0, str(_p.Path(__file__).resolve().parents[2] / "src")); import disatlas.deterministic  # noqa: E402,F401  reproducible PDF/SVG bytes
import json
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

from pathlib import Path
HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]
import sys as _sys; _sys.path.insert(0, str(REPO / "src"))
DV = REPO / "data" / "derived" / "popularity_strata"
SUMMARY, BOUNDS = DV / "summary_popularity_strata.tsv", DV / "strata_boundaries.tsv"
DECISION, DTABLE = DV / "decision_popularity_strata.json", DV / "decision_table_all_cells.tsv"
OUT = HERE
SRC_S, SRC_B, SRC_J, SRC_D = (Path(p).name for p in (SUMMARY, BOUNDS, DECISION, DTABLE))

S = pd.read_csv(SUMMARY, sep="\t")
BND = pd.read_csv(BOUNDS, sep="\t")
DEC = json.load(open(DECISION))
DT = pd.read_csv(DTABLE, sep="\t")
# sign convention for every figure: positive = DisAtlas better
_m = S.endpoint.str.startswith("med_")
S.loc[_m, "delta"] = -S.loc[_m, "delta"]
S.loc[_m, ["delta_lo", "delta_hi"]] = -S.loc[_m, ["delta_hi", "delta_lo"]].values

PRIMARY = (2010, "test")
CELLS = [(2005, "test"), (2005, "all"), (2010, "test"), (2010, "all"), (2015, "test"), (2015, "all")]
COMPARATORS = [  # (method id, display label, family)
    ("adamic_adar_w", "Weighted Adamic–Adar", "graph"),
    ("adamic_adar", "Adamic–Adar", "graph"),
    ("degree", "Degree", "graph"),
    ("strength", "Strength (weighted degree)", "graph"),
    ("gobp_simgic", "GO-BP semantic similarity", "annot"),
    ("gene_jaccard", "Shared genes (Jaccard)", "annot"),
    ("hpo_simgic", "HPO phenotype similarity", "annot"),
]
CUT_COL = {2005: "#E69F00", 2010: "#0072B2", 2015: "#009E73"}  # Okabe–Ito
BAND = "#DCE9F5"
GREY = "#555555"

# ------------------------------------------------------------------ numbers log
LOG = []


def row(t, pool, endpoint, method):
    r = S[(S.t == t) & (S.pool == pool) & (S.endpoint == endpoint) & (S.method == method)]
    assert len(r) == 1, (t, pool, endpoint, method, len(r))
    return r.index[0], r.iloc[0]


def take(use, t, pool, endpoint, method, col):
    idx, r = row(t, pool, endpoint, method)
    v = r[col]
    assert pd.notna(v), (t, pool, endpoint, method, col)
    LOG.append(dict(figure="SI", used_in=use, value=v, source=SRC_S, source_row=int(idx) + 2,  # 1-based incl. header
                    t=t, pool=pool, endpoint=endpoint, method=method, column=col))
    return v


def take_b(use, t, stratum, col):
    r = BND[(BND.t == t) & (BND.stratum == stratum)]
    assert len(r) == 1
    v = r.iloc[0][col]
    LOG.append(dict(figure="SI", used_in=use, value=v, source=SRC_B, source_row=int(r.index[0]) + 2,
                    t=t, pool="", endpoint=f"stratum {stratum}", method="", column=col))
    return v


def ci3(use, t, pool, endpoint, method, prefix="delta"):
    lo, hi = (prefix + "_lo", prefix + "_hi") if prefix == "delta" else ("lo", "hi")
    return tuple(take(use, t, pool, endpoint, method, c) for c in (prefix, lo, hi))


# ------------------------------------------------------------------ style
mpl.rcParams.update({
    "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "legend.fontsize": 6,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.spines.top": False,
    "axes.spines.right": False, "xtick.direction": "out", "ytick.direction": "out",
    "legend.frameon": False, "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

MM = 1 / 25.4
fig = plt.figure(figsize=(180 * MM, 200 * MM))
axA = fig.add_axes([0.215, 0.12, 0.265, 0.81])
axB = fig.add_axes([0.725, 0.715, 0.26, 0.215])
axC = fig.add_axes([0.725, 0.42, 0.255, 0.17])
axD = fig.add_axes([0.725, 0.12, 0.255, 0.19])
TITLE = {}  # axes -> (letter x, title text)


def marker_kw(t, pool, primary):
    c = CUT_COL[t]
    return dict(marker="o" if pool == "test" else "s", ms=4.2 if primary else 3.2,
                mfc=c if pool == "test" else "white", mec="black" if primary else c,
                mew=1.0 if primary else 0.8, color=c, lw=1.6 if primary else 0.9, zorder=3)


def forest_point(ax, y, est, lo, hi, kw):
    ax.plot([lo, hi], [y, y], color=kw["color"], lw=kw["lw"], solid_capstyle="butt", zorder=2)
    ax.plot([est], [y], ls="none", **{k: v for k, v in kw.items() if k not in ("color", "lw")})


# ------------------------------------------------------------------ panel a
step, gap = 1.0, 1.6
y = 0.0
yticks, ylabels, group_spans = [], [], []
for mid, label, fam in COMPARATORS:
    ys = []
    for (t, pool) in CELLS:
        prim = (t, pool) == PRIMARY
        d, lo, hi = ci3(f"a: {label}, {t} {pool}", t, pool, "med_strat", mid)
        if prim:
            axA.axhspan(y - 0.5, y + 0.5, color=BAND, zorder=0, lw=0)
        forest_point(axA, y, d, lo, hi, marker_kw(t, pool, prim))
        ys.append(y)
        y += step
    yticks.append(np.mean(ys)); ylabels.append(label); group_spans.append((ys[0], ys[-1], fam))
    y += gap
axA.axvline(0, color="black", lw=0.7, zorder=1)
axA.set_yticks(yticks, ylabels)
axA.tick_params(axis="y", length=0)
axA.set_ylim(y - gap + 0.2, -1.4)
axA.set_xlabel("DisAtlas advantage over reference score\n(within-stratum median percentile, points)")
axA.margins(x=0.04)
xl = axA.get_xlim(); axA.set_xlim(min(xl[0], -1.5), xl[1])
# family brackets (left, outside tick labels)
for fam, txt in (("graph", "Disease-graph scores"), ("annot", "Annotation-based")):
    spans = [s for s in group_spans if s[2] == fam]
    y0, y1 = spans[0][0] - 0.4, spans[-1][1] + 0.4
    axA.annotate("", xy=(-0.70, y0), xycoords=("axes fraction", "data"), xytext=(-0.70, y1),
                 textcoords=("axes fraction", "data"),
                 arrowprops=dict(arrowstyle="-", color=GREY, lw=0.7), annotation_clip=False)
    axA.text(-0.735, (y0 + y1) / 2, txt, transform=mpl.transforms.blended_transform_factory(axA.transAxes, axA.transData),
             rotation=90, ha="center", va="center", fontsize=6, color=GREY)
TITLE[axA] = (0.01, "Within popularity strata, DisAtlas ranks later-trialled indications\nhigher than every reference score, in every cell")
axA.text(0.98, 0.005, "DisAtlas better →", transform=axA.transAxes, fontsize=6, color=GREY, ha="right", va="bottom")

# legend (cutoff colour, pool marker, primary band)
h = [Line2D([], [], ls="none", marker="o", mfc=CUT_COL[t], mec=CUT_COL[t], ms=3.5, label=f"Evidence up to {t}") for t in CUT_COL]
h += [Line2D([], [], ls="none", marker="o", mfc=GREY, mec=GREY, ms=3.5, label="Test-pool diseases"),
      Line2D([], [], ls="none", marker="s", mfc="white", mec=GREY, ms=3.5, label="All diseases"),
      Patch(facecolor=BAND, edgecolor="black", lw=0.6, label="Primary analysis (2010, test pool)")]
fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, -0.01), ncol=3, fontsize=6, handletextpad=0.4, columnspacing=2.0)

# ------------------------------------------------------------------ panel b
t, pool = PRIMARY
kw = marker_kw(t, pool, True)
bt, bl = [], []
for k in range(5):
    ep = f"med_strat_k{k}"
    d, lo, hi = ci3(f"b: stratum k{k}", t, pool, ep, "adamic_adar_w")
    nd = take(f"b: stratum k{k} n drugs", t, pool, ep, "adamic_adar_w", "n_drugs")
    nl = take(f"b: stratum k{k} n later", t, pool, ep, "adamic_adar_w", "n_later")
    pmin = take_b(f"b: stratum k{k} popularity min", t, k, "pop_min")
    pmax = take_b(f"b: stratum k{k} popularity max", t, k, "pop_max")
    forest_point(axB, k, d, lo, hi, kw)
    rng = "0" if pmin == pmax else f"{int(pmin)}–{int(pmax):,}"
    bt.append(k); bl.append(f"Popularity {rng}\n{int(nd):,} drugs, {int(nl):,} later indications")
axB.axvline(0, color="black", lw=0.7, zorder=1)
axB.set_yticks(bt, bl); axB.tick_params(axis="y", length=0)
axB.set_ylim(4.9, -0.5)
axB.margins(x=0.06)
xl = axB.get_xlim(); axB.set_xlim(min(xl[0], -0.8), xl[1])
axB.xaxis.set_major_locator(mpl.ticker.MultipleLocator(4))
axB.set_xlabel("DisAtlas advantage over weighted\nAdamic–Adar (points)")
TITLE[axB] = (0.53, "The advantage over weighted Adamic–Adar holds in every\npopularity stratum (2010, test pool)")
axB.text(0.99, 0.01, "DisAtlas better →", transform=axB.transAxes, fontsize=6, color=GREY, ha="right", va="bottom")

# ------------------------------------------------------------------ panel c
ct, cl = [], []
for i, (t, pool) in enumerate(CELLS):
    prim = (t, pool) == PRIMARY
    d, lo, hi = ci3(f"c: fused vs popularity, {t} {pool}", t, pool, "med_fused", "popularity")
    if prim:
        axC.axhspan(i - 0.5, i + 0.5, color=BAND, zorder=0, lw=0)
    forest_point(axC, i, d, lo, hi, marker_kw(t, pool, prim))
    ct.append(i); cl.append(f"{t}, {'test pool' if pool == 'test' else 'all diseases'}")
axC.axvline(0, color="black", lw=0.7, zorder=1)
axC.set_yticks(ct, cl); axC.tick_params(axis="y", length=0)
axC.set_ylim(6.3, -0.6)
axC.margins(x=0.08)
axC.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
xl = axC.get_xlim(); axC.set_xlim(min(xl[0], -0.37), xl[1])
axC.set_xlabel("DisAtlas + popularity vs\npopularity alone (points)")
TITLE[axC] = (0.53, "Combining DisAtlas with popularity does not rank higher\nthan popularity alone in the test pool")
axC.text(0.015, 6.2, "combination better →", transform=axC.transData, fontsize=6, color=GREY, ha="left", va="bottom")
axC.text(-0.015, 6.2, "← popularity better", transform=axC.transData, fontsize=6, color=GREY, ha="right", va="bottom")

# ------------------------------------------------------------------ panel d (top-10% precision, %)
PREC = [("popularity", "Trial popularity"), ("adamic_adar_w", "Weighted\nAdamic–Adar"), ("disatlas_ens5", "DisAtlas")]
y = 0.0; dt, dl = [], []
for mid, label in PREC:
    ys = []
    for (t, pool) in CELLS:
        prim = (t, pool) == PRIMARY
        v, lo, hi = ci3(f"d: mean top-10% precision, {label.replace(chr(10), ' ')}, {t} {pool}", t, pool, "prec10_mean", mid, prefix="value")
        if prim:
            axD.axhspan(y - 0.5, y + 0.5, color=BAND, zorder=0, lw=0)
        forest_point(axD, y, 100 * v, 100 * lo, 100 * hi, marker_kw(t, pool, prim))
        ys.append(y); y += 1.0
    dt.append(np.mean(ys)); dl.append(label); y += 1.2
axD.set_yticks(dt, dl); axD.tick_params(axis="y", length=0)
for tl, (mid, _) in zip(axD.get_yticklabels(), PREC):
    if mid == "disatlas_ens5":
        tl.set_fontweight("bold")
axD.set_ylim(y - 1.2 - 0.4, -0.6)
axD.set_xlim(0, None)
axD.set_xlabel("Top-10% precision (%; higher = better)")
TITLE[axD] = (0.53, "Trial popularity has the highest top-10% precision in\nevery cell; DisAtlas and weighted Adamic–Adar are close")

for ax, L in ((axA, "a"), (axB, "b"), (axC, "c"), (axD, "d")):
    x0, txt = TITLE[ax]
    top = ax.get_position().y1
    fig.text(x0, top + 0.012, L, fontsize=9, fontweight="bold", va="bottom", ha="left")
    fig.text(x0 + 0.025, top + 0.012, txt, fontsize=7, va="bottom", ha="left")

# claim checks behind titles
assert (S[(S.endpoint == "med_strat") & S.method.isin([m for m, _, _ in COMPARATORS])].delta_lo > 0).all()  # a: every cell, every comparator
assert (S[(S.endpoint.str.startswith("med_strat_k")) & (S.method == "adamic_adar_w") & (S.t == 2010) & (S.pool == "test")].delta_lo > 0).all()  # b
_p = S[S.endpoint == "prec10_mean"].pivot_table(index=["t", "pool"], columns="method", values="value")
assert (_p.idxmax(axis=1) == "popularity").all()  # d
_f = S[(S.endpoint == "med_fused") & (S.method == "popularity") & (S.pool == "test")]
assert (_f.delta_hi > 0).sum() == 1 and (_f.delta < 0).all()  # c: test pool never above popularity; one interval spans zero
pd.DataFrame(LOG).to_csv(OUT / "si_popularity_strata_numbers.tsv", sep="\t", index=False)
# ------------------------------------------------------------------ save + bbox check
for ext in ("png", "pdf"):
    fig.savefig(OUT / f"si_popularity_strata.{ext}", bbox_inches="tight")
r = fig.canvas.get_renderer()
texts = [(t_, t_.get_window_extent(r)) for t_ in fig.findobj(mpl.text.Text) if t_.get_text().strip() and t_.get_visible()]
ov = [(a.get_text(), b.get_text()) for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
fb = fig.bbox
outside = [t_.get_text() for t_, bb in texts if bb.x0 < fb.x0 - 1 or bb.y0 < fb.y0 - 1 or bb.x1 > fb.x1 + 1 or bb.y1 > fb.y1 + 1]
print("text overlaps:", ov)
print("outside figure:", outside)
