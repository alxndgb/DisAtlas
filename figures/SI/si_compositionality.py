"""SI figure preview: compositionality of DisAtlas space (current release, five seeds).
Every value is read from data/derived/si_compositionality/ (written by si_compositionality_recompute.py):
  systematic_mondo_multiparent.tsv  per-disease ranks, 2,535 MONDO diseases with >= 2 embedded is_a parents (all drawn)
  systematic_summary.tsv            top-1% shares with Wilson 95% CIs
  examples.tsv / examples_null_draws.tsv.gz   both worked examples and all 1,000 null draws each
Panel a: share of diseases ranked in the nearest 1% to each query, 95% CI. Panel b: cumulative distribution of the
disease's percentile rank for each query (all diseases). The Supplementary Figure has panels a and b only.
Usage: python figures/SI/si_compositionality.py [--with-panel-c]
  default         panels a and b (the Supplementary Figure); writes SI_compositionality.{png,pdf}
  --with-panel-c  also draws both worked examples against their random-combination nulls; writes SI_compositionality_abc.*
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import sys, os, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from disatlas import figstyle_si as fs
EV, OUT = DERIVED / "si_compositionality", HERE; NOPC = "--with-panel-c" not in sys.argv[1:]; os.makedirs(OUT, exist_ok=True); fs.apply_figure_style()
R = pd.read_csv(f"{EV}/systematic_mondo_multiparent.tsv", sep="\t")
S = pd.read_csv(f"{EV}/systematic_summary.tsv", sep="\t").set_index("query")
E = pd.read_csv(f"{EV}/examples.tsv", sep="\t"); Nd = pd.read_csv(f"{EV}/examples_null_draws.tsv.gz", sep="\t")
n = len(R); assert n == int(S.n.iloc[0]) and len(E) == 2 and len(Nd) == 2000
Q = [("Random parents", "null_random_parents", "Random diseases", "#999999"),
     ("One real parent + random others", "null_one_real_rest_random", "One real parent + random diseases", "#E69F00"),
     ("Barycentre of real parents", "bary", "Average of all real parents", "#0072B2"),
     ("Best single real parent", "best_single_parent", "Closest single real parent", "#009E73")]
for q, col, _, _ in Q:  # summary table and per-disease table agree
    assert int((R[f"pct_{col}"] <= .01).sum()) == int(S.loc[q, "top1pct_k"])
if NOPC:
    fig = plt.figure(figsize=(7.2, 2.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=.42, bottom=.30)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1])
else:
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1, 1.25], hspace=.62, wspace=.42)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1])
    gc = gs[1, :].subgridspec(1, 2, wspace=.22); c1 = fig.add_subplot(gc[0]); c2 = fig.add_subplot(gc[1], sharey=c1)
y = np.arange(len(Q))[::-1]
for k, (q, col, lab, cc) in enumerate(Q):
    s = S.loc[q]
    a.errorbar(100 * s.top1pct, y[k], xerr=[[100 * (s.top1pct - s.top1pct_lo)], [100 * (s.top1pct_hi - s.top1pct)]], fmt="o", ms=4, color=cc, lw=1, capsize=0)
    a.text(100 * s.top1pct_hi + 2, y[k], f"{100 * s.top1pct:.1f}%", va="center", fontsize=6.5, color=cc)
a.set_yticks(y); a.set_yticklabels([q[2] for q in Q]); a.set_xlim(0, 70)
a.set_xlabel(f"Diseases ranked in the nearest 1% (%)"); a.set_title(f"a  Query built from a disease's MONDO parents (n = {n:,})", loc="left", x=-1.05)
grid = np.logspace(-4, 0, 400)
for q, col, lab, cc in Q:
    v = np.sort(R[f"pct_{col}"].to_numpy()); b.plot(grid, 100 * np.searchsorted(v, grid, side="right") / n, color=cc, lw=1.2, label=lab)
b.axvline(.01, color="0.4", ls="--", lw=.7); b.set_xscale("log"); b.set_xlim(1e-4, 1); b.set_ylim(0, 100)
b.set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1]); b.set_xticklabels(["0.01%", "0.1%", "1%", "10%", "100%"])
b.set_xlabel("Rank of the disease among all diseases (percentile, nearest first)"); b.set_ylabel("Diseases (cumulative %)")
b.set_title("b  All diseases, each query", loc="left")
if not NOPC:
    for ax, (_, e) in zip([c1, c2], E.iterrows()):
        d = Nd[Nd.example == e.example].null_pct.to_numpy(); bins = np.logspace(-4, 0, 41)
        ax.hist(np.clip(d, 1e-4, 1), bins=bins, color="0.8", edgecolor="white", lw=.3)
        ax.axvline(e.pct, color="#0072B2", lw=1.4); ax.set_xscale("log"); ax.set_xlim(1e-4, 1)
        ax.set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1]); ax.set_xticklabels(["0.01%", "0.1%", "1%", "10%", "100%"])
        ax.text(e.pct * .8, .97, f"rank {int(e['rank']):,} of {int(e.out_of):,}\np = {e.p:.3f}", transform=ax.get_xaxis_transform(), ha="right", va="top", fontsize=6.5, color="#0072B2")
        ax.set_xlabel("Rank of the target disease (percentile)")
        ax.set_title(f"{e.target[0].upper() + e.target[1:]} ({int(e.n_components)} components)", loc="left", fontsize=7)
    c1.set_ylabel("Random combinations (of 1,000)"); plt.setp(c2.get_yticklabels(), visible=False)
    fig.text(c1.get_position().x0 - .07, c1.get_position().y1 + .065, "c  The two worked examples against random combinations of as many diseases", fontsize=8)
h, l = b.get_legend_handles_labels(); fig.legend(h[::-1], l[::-1], loc="upper center", bbox_to_anchor=(.5, .14 if NOPC else .03), ncol=4, handlelength=1.6, columnspacing=1.2, fontsize=6.5)
J = json.load(open(f"{EV}/compositionality_summary.json"))
vals = dict(n=n, top1={q: dict(k=int(S.loc[q, "top1pct_k"]), share=float(S.loc[q, "top1pct"]), lo=float(S.loc[q, "top1pct_lo"]), hi=float(S.loc[q, "top1pct_hi"]),
            median_pct=float(S.loc[q, "median_pct"])) for q, *_ in Q},
            bary_better_than_best_single=J["share_bary_better_than_best_single"], bary_ties=J["share_bary_ties_best_single"],
            bary_worse=J["share_bary_worse_than_best_single"], examples=E.to_dict("records"))
vals["png_bytes"] = fs.save(fig, OUT, "SI_compositionality" + ("" if NOPC else "_abc")); vals["overlaps"] = fs.overlap_check(fig)
json.dump(vals, open(f"{OUT}/SI_compositionality_preview{'' if NOPC else '_abc'}_values.json", "w"), indent=1, default=float); print(json.dumps(vals, default=float)[:1500])
