#!/usr/bin/env python
"""SI figure: rare-common neighbour coverage in the DisAtlas embedding, current model.

Re-render of FIGURES/SI_rare_common_neighbour_coverage.png with the same layout, from the current-model outputs of
the original analysis code (rc coverage), out/rc_coverage:
  a  fraction of rare diseases with >= 1 common disease among their k nearest neighbours, k = 1..50, observed vs three
     label-permutation nulls (1000 permutations; mean and 95% band): labels shuffled, shuffled within annotation-depth
     bins, shuffled within connectivity (training-graph degree) bins  [coverage_by_k_summary.tsv]
  b  same for common diseases with >= 1 rare neighbour  [coverage_by_k_summary.tsv]
  c  distribution of the number of common diseases among the 10 nearest neighbours of each rare disease
     [per_disease_first_other_class.tsv.gz, column n_other_class_in_10]
Outputs: next to this script: SI_rare_common_neighbour_coverage.{png,pdf} and *_values.tsv
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import os, sys, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter

from disatlas.disatlas_style import apply_style, panel_label

D = str(DERIVED / "si_rare_common")
SUMF, PERF = (f"{D}/coverage_by_k_summary.tsv", f"{D}/per_disease_first_other_class.tsv.gz")
OUTDIR = str(HERE)
NAME = "SI_rare_common_neighbour_coverage"

S = pd.read_csv(SUMF, sep="\t")
P = pd.read_csv(PERF, sep="\t")
rare = P[P.is_rare == 1]; common = P[P.is_rare == 0]
# consistency: observed curve at k = 10 equals the per-disease table
for cls, sub in (("rare", rare), ("common", common)):
    obs10 = S[(S.model == "observed") & (S.anchor_class == cls) & (S.k == 10)]["mean"].item()
    assert abs(obs10 - (sub.first_other_class_rank <= 10).mean()) < 1e-9, cls

BLUE = "#2160a8"; GREY = "0.6"
NULLS = [("null_annotation", ":", "chance, matched on annotation depth"),
         ("null_degree", "--", "chance, matched on connectivity"),
         ("null_uniform", "-", "chance (labels shuffled)")]
apply_style()
fig, ax = plt.subplots(1, 3, figsize=(10.6, 4.5), gridspec_kw=dict(width_ratios=[1, 1, 0.82], wspace=0.55))
VAL = []
for a, cls, other, ttl in ((ax[0], "rare", "common", "Rare diseases"), (ax[1], "common", "rare", "Common diseases")):
    for m, ls, _ in NULLS:
        g = S[(S.model == m) & (S.anchor_class == cls)].sort_values("k")
        a.fill_between(g.k, g.lo, g.hi, color="0.85", lw=0, zorder=1)
        a.plot(g.k, g["mean"], ls=ls, color=GREY, lw=1.6, zorder=2)
        r10 = g[g.k == 10].iloc[0]
        VAL.append((ttl, m, 10, r10["mean"], r10.lo, r10.hi, r10.p_lower))
    o = S[(S.model == "observed") & (S.anchor_class == cls)].sort_values("k")
    a.plot(o.k, o["mean"], color=BLUE, lw=2.6, zorder=3)
    y10 = o[o.k == 10]["mean"].item()
    a.plot(10, y10, "o", color=BLUE, ms=6, zorder=4)
    xt, yt = (17, y10 - 0.14) if cls == "rare" else (17, y10 - 0.29)
    a.annotate(f"{y10 * 100:.0f}% at k = 10", (10, y10), (xt, yt), color=BLUE, fontsize=9.5, va="center",
               arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.9, shrinkA=0, shrinkB=3))
    for k in (1, 5, 10, 20, 50):
        VAL.append((ttl, "observed", k, o[o.k == k]["mean"].item(), np.nan, np.nan, np.nan))
    a.set_xlim(0, 51); a.set_ylim(0, 1.03); a.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    a.set_xlabel("Nearest neighbors considered (k)")
    a.set_ylabel(f"{cls.capitalize()} diseases with ≥1 {other}\ndisease among neighbors")
    a.set_title(ttl, loc="left"); a.tick_params(labelsize=8.5)

c = ax[2]
frac = rare.n_other_class_in_10.value_counts(normalize=True).reindex(range(11), fill_value=0)
c.bar(frac.index, frac.values, width=0.75, color=BLUE)
c.annotate(f"none: {frac[0] * 100:.0f}%", (0.4, frac[0] - 0.005), (2.5, frac[0] - 0.025), color=BLUE, fontsize=9.5,
           va="center", arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.9, shrinkA=0, shrinkB=0))
for n_, v in frac.items():
    VAL.append(("Rare diseases, k = 10", f"{n_} common among 10 nearest", 10, v, np.nan, np.nan, np.nan))
c.set_xticks(range(0, 11, 2)); c.set_xlim(-0.7, 10.7)
c.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
c.set_xlabel("Common diseases among\nthe 10 nearest neighbors", fontsize=11); c.set_ylabel("Rare diseases")
c.set_title("Rare diseases, k = 10", loc="left"); c.tick_params(labelsize=8.5)
for a_, l in zip(ax, "abc"):
    panel_label(a_, l, xy=(-0.18 if l != "c" else -0.28, 1.03))

h = [Line2D([], [], color=BLUE, lw=2.6, label="DisAtlas")] + \
    [Line2D([], [], color=GREY, ls=ls, lw=1.6, label=lab) for _, ls, lab in NULLS]
fig.legend(handles=h, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.47, -0.06), fontsize=8.5)

os.makedirs(OUTDIR, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(f"{OUTDIR}/{NAME}.{ext}", dpi=300, bbox_inches="tight")
V = pd.DataFrame(VAL, columns=["panel", "series", "k", "value", "lo", "hi", "p_lower"]).round(4)
med = float(np.median(rare.first_other_class_rank))
extra = pd.DataFrame([("text", "median rank of first common neighbor (rare)", np.nan, med, np.nan, np.nan, np.nan),
                      ("text", "n rare / n common", np.nan, len(rare), np.nan, len(common), np.nan)], columns=V.columns)
pd.concat([V, extra]).to_csv(f"{OUTDIR}/{NAME}_values.tsv", sep="\t", index=False)
print(f"{OUTDIR}/{NAME}.png")
