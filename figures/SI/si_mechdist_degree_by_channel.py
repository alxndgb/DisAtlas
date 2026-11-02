"""Supplementary Figure (fig:SI_degree_by_channel): number of partner diseases in each graph of the MechDist network.

  a-c  Share of diseases with degree >= k (CCDF) for rare, common and all diseases in the broad channel, the specific
       channel and the combined graph; dotted lines at each group's median (over diseases with >= 1 edge where the
       median is 0); dashed black = discrete power law fitted to the tail k >= xmin,
       P(K >= k | K >= xmin) = zeta(alpha, k) / zeta(alpha, xmin), scaled by n_tail / N; tick at xmin.
  d-f  p(k) over diseases with >= 1 edge, all diseases pooled: logarithmic bins (about 10 per decade, integer edges,
       density, arithmetic bin centres) and each single degree (faint grey); dashed = fitted pmf
       k^-alpha / zeta(alpha, xmin) scaled by n_tail / N_nonzero; dotted = xmin.  The fit uses the unbinned degrees.
Binning and fit lines follow the powerlaw package (Alstott et al. 2014; statistics.pdf with linear_bins=False and
Power_Law.plot_ccdf / plot_pdf evaluated at the unique tail degrees), re-implemented with numpy/scipy so that the
figure needs no extra dependency.  xmin, sigma, the bootstrap goodness-of-fit p and the likelihood-ratio tests are
read from the shipped fit table; alpha is refitted by discrete maximum likelihood at that xmin as a check (|diff| < 5e-4).

Inputs (data/derived/):
  degree_distributions/degree_per_disease.tsv   degree of each of the 13,021 diseases in the three graphs, and class
  degree_distributions/degree_summary.tsv       per graph and class: median / max degree, isolated diseases, top-1% share
  degree_distributions/powerlaw_summary.tsv     power-law tail fits (fit == "free") and likelihood ratios
  training_graph/training_graph_stats.json      statistics of the combined (training) graph
Run:  python figures/SI/si_mechdist_degree_by_channel.py
Writes si_mechdist_degree_by_channel.png/.pdf and si_mechdist_degree_by_channel_values.json next to this script.
"""
import sys as _s, pathlib as _p; _s.path.insert(0, str(_p.Path(__file__).resolve().parents[2] / "src")); import disatlas.deterministic  # noqa: E402,F401  reproducible PDF/SVG bytes
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.special import zeta
from scipy.optimize import minimize_scalar
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from disatlas.disatlas_style import apply_style, panel_label, ANNOT_SIZE, DPI

HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]
DD = REPO / "data" / "derived" / "degree_distributions"
RARE_COL, COMMON_COL, ALL_COL, FIT_COL, BIN_COL = "#d1495b", "#2a6ebb", "0.25", "black", "#2a78d6"
NAMES = {"broad": "Broad channel", "specific": "Specific channel", "combined": "Combined graph"}
GRAPHS = ["broad", "specific", "combined"]
ALTS = {"lognormal": "lognormal", "exponential": "exponential", "stretched_exponential": "stretched exp.",
        "truncated_power_law": "truncated power law"}
E_STATED = {"broad": 331234, "specific": 190007, "combined": 498155}

P = pd.read_csv(DD / "degree_per_disease.tsv", sep="\t")
S = pd.read_csv(DD / "degree_summary.tsv", sep="\t")
PS = pd.read_csv(DD / "powerlaw_summary.tsv", sep="\t"); PS = PS[PS.fit == "free"].set_index("graph")
TG = json.load(open(REPO / "data" / "derived" / "training_graph" / "training_graph_stats.json"))
rare = (P.disease_class == "rare").to_numpy(); common = (P.disease_class == "common").to_numpy()
deg = {g: P[f"deg_{g}"].to_numpy() for g in GRAPHS}; n = len(P)
for g in GRAPHS:
    assert deg[g].sum() // 2 == E_STATED[g], (g, deg[g].sum() // 2)
dc = deg["combined"]          # panel c is the training graph of Methods
assert (n, int(dc.sum() // 2), int(dc.min()), int(dc.max()), int((dc == 0).sum())) == \
    (TG["n_nodes"], TG["n_edges"], TG["min_degree"], TG["max_degree"], TG["n_isolated"])
assert abs(dc.mean() - TG["mean_degree"]) < 1e-9 and float(np.median(dc)) == TG["median_degree"]


def ccdf(x):
    k = np.unique(x[x > 0]); xs = np.sort(x)
    return k, 1 - np.searchsorted(xs, k, side="left") / len(x)


def log_binned_pdf(x):
    """As powerlaw.statistics.pdf(x, linear_bins=False): ceil(10 x decades) log-spaced edges, floored, density."""
    x = x[x > 0]; lo, hi = np.log10(x.min()), np.log10(x.max())
    bins = np.logspace(lo, hi, num=int(np.ceil((hi - lo) * 10)))
    bins[:-1] = np.floor(bins[:-1]); bins[-1] = np.ceil(bins[-1])
    return np.histogram(x, np.unique(bins), density=True)[::-1]


def fit_alpha(x, xmin):
    """Discrete power-law maximum-likelihood alpha for the tail x >= xmin, alpha in [1, 10] as in the shipped fits."""
    t = x[x >= xmin].astype(float); sl = np.log(t).sum(); m = len(t)
    r = minimize_scalar(lambda a: m * np.log(zeta(a, xmin)) + a * sl, bounds=(1.0001, 10.0), method="bounded",
                        options=dict(xatol=1e-10))
    return float(r.x), m


def pl_ccdf(k, a, xmin):
    return zeta(a, k) / zeta(a, xmin)


def pl_pdf(k, a, xmin):
    return k ** -a / zeta(a, xmin)


def pfmt(p):
    if p < 1e-15: return "< 10⁻¹⁵"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"{m}×10" + str(int(e)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{p:.2f}" if p >= 0.01 else f"{p:.3f}"


FIT = {}
for g in GRAPHS:
    x = deg[g][deg[g] > 0]; r = PS.loc[g]; xm = float(r.xmin)
    a_ml, m = fit_alpha(x, xm); a = float(r.alpha)   # draw the shipped fit; the refit checks it
    assert m == int(r.n_tail) and abs(a_ml - a) < 5e-4, (g, a_ml, a, m, int(r.n_tail))
    best = min(ALTS, key=lambda q: r[f"LR_{q}_R"])
    FIT[g] = dict(alpha=a, sigma=float(r.sigma), xmin=xm, n_tail=m, p_gof=float(r.p_gof), n_boot=int(r.n_boot),
                  n_boot_D_ge_obs=int(r.n_boot_D_ge_obs), best=best, best_R=float(r[f"LR_{best}_R"]),
                  best_p=float(r[f"LR_{best}_p"]))

apply_style()
fig, AX = plt.subplots(2, 3, figsize=(15, 12.2), sharex=True, sharey="row", gridspec_kw=dict(hspace=0.55, wspace=0.10))
axs, bxs = AX[0], AX[1]
xmax = max(int(deg[g].max()) for g in GRAPHS)
TXT = dict(fontsize=ANNOT_SIZE, color="0.3", va="top", ha="left", linespacing=1.35)
V = {"graphs": {}}
# ---------------- top row: CCDF by group + fitted power-law CCDF ----------------
for ax, lab, g in zip(axs, "abc", GRAPHS):
    d = deg[g]
    for grp, msk, col, lw, z in (("all", np.ones(n, bool), ALL_COL, 1.3, 3), ("common", common, COMMON_COL, 2.0, 2),
                                 ("rare", rare, RARE_COL, 2.0, 1)):
        k, y = ccdf(d[msk])
        ax.step(k, y, where="post", color=col, lw=lw, ls=(0, (4, 2)) if grp == "all" else "-", zorder=z)
        med = float(np.median(d[msk]))
        if med < 1: med = float(np.median(d[msk][d[msk] > 0]))
        if grp != "all": ax.axvline(med, color=col, lw=0.9, ls=":", zorder=0)
    f = FIT[g]; a_, xm = f["alpha"], f["xmin"]
    kt = np.unique(d[d >= xm]).astype(float)
    ax.plot(kt, pl_ccdf(kt, a_, xm) * f["n_tail"] / n, color=FIT_COL, lw=1.6, ls="--", zorder=5)
    ax.plot([xm, xm], [0.5 / n, 1.4 / n], color=FIT_COL, lw=1.0, zorder=5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.8, xmax * 1.4); ax.set_ylim(0.5 / n, 1.3)
    ax.grid(True, which="major", color="0.9", lw=0.6); ax.set_axisbelow(True)
    ax.set_title(f"{NAMES[g]}: {E_STATED[g]:,} edges", fontsize=12)
    ax.set_xlabel("Degree k (number of partner diseases)"); ax.xaxis.set_tick_params(labelbottom=True)
    s = S[S.graph == g].set_index("group")

    def f0(grp):
        k_, v = int(s.loc[grp, "n_isolated"]), s.loc[grp, "pct_isolated"]
        return "0" if k_ == 0 else (f"{v:.1f}%" if v >= 1 else f"{k_:,}")
    if s.loc["all", "median_degree"] < 1:
        med_l = "median, diseases with ≥ 1 edge: " + " / ".join(f"{s.loc[q, 'median_degree_connected']:.0f}" for q in ("all", "rare", "common"))
    else:
        med_l = "median degree: " + " / ".join(f"{s.loc[q, 'median_degree']:.0f}" for q in ("all", "rare", "common"))
    txt = ["all / rare / common", med_l,
           f"degree 0: {f0('all')} / {f0('rare')} / {f0('common')}" + ("" if s.loc["all", "pct_isolated"] >= 1 else " diseases"),
           f"max degree {int(s.loc['all', 'max_degree']):,}",
           f"top 1% of diseases hold {s.loc['all', 'pct_endpoints_top1pct']:.1f}% of edge ends"]
    ax.text(0.0, -0.17, "\n".join(txt), transform=ax.transAxes, **TXT)
    panel_label(ax, lab)
    V["graphs"][g] = dict(edges=E_STATED[g], n_isolated=int((d == 0).sum()), max_degree=int(d.max()), **f)
axs[0].set_ylabel("Share of diseases with degree ≥ k")
axs[0].legend(handles=[Line2D([], [], color=RARE_COL, lw=2, label=f"Rare ({rare.sum():,})"),
                       Line2D([], [], color=COMMON_COL, lw=2, label=f"Common ({common.sum():,})"),
                       Line2D([], [], color=ALL_COL, lw=1.3, ls=(0, (4, 2)), label=f"All ({n:,})"),
                       Line2D([], [], color="0.4", lw=0.9, ls=":", label="Group median (b: ≥ 1 edge)"),
                       Line2D([], [], color=FIT_COL, lw=1.6, ls=(0, (6, 3)), label="Power-law fit, k ≥ xmin")],
              loc="lower left", frameon=False, fontsize=9.5)
# ---------------- bottom row: log-binned p(k) + fitted power-law pmf ----------------
for bx, lab, g in zip(bxs, "def", GRAPHS):
    d = deg[g]; x = d[d > 0]; nz = len(x); n0 = int((d == 0).sum())
    kv, cnt = np.unique(x, return_counts=True)
    bx.scatter(kv, cnt / nz, s=4, color="0.55", alpha=0.12, edgecolors="none", zorder=1)
    edges, hist = log_binned_pdf(x)
    cen = (edges[1:] + edges[:-1]) / 2; hp = hist.astype(float).copy(); hp[hp == 0] = np.nan
    bx.plot(cen, hp, color=BIN_COL, lw=0.8, marker="o", ms=5.5, mec="white", mew=0.5, zorder=4)
    f = FIT[g]; a_, xm = f["alpha"], f["xmin"]
    kt = np.unique(x[x >= xm]).astype(float)
    bx.plot(kt, pl_pdf(kt, a_, xm) * f["n_tail"] / nz, color=FIT_COL, lw=1.6, ls="--", zorder=5)
    bx.axvline(xm, color=FIT_COL, lw=0.9, ls=":", zorder=1)
    bx.set_xscale("log"); bx.set_yscale("log")
    bx.grid(True, which="major", color="0.9", lw=0.6); bx.set_axisbelow(True)
    bx.set_title(f"{NAMES[g]}: all diseases pooled", fontsize=12)
    bx.set_xlabel("Degree k (number of partner diseases)")
    pg = f"< {1 / f['n_boot']:.3f}" if f["n_boot_D_ge_obs"] == 0 else f"= {f['p_gof']:.3f}"
    bx.text(0.97, 0.97, f"α = {a_:.2f} ± {f['sigma']:.2f}\nxmin = {xm:.0f}\np(GoF) {pg}", transform=bx.transAxes,
            fontsize=ANNOT_SIZE, color="0.15", va="top", ha="right", linespacing=1.3,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.85", lw=0.6))
    bp = pfmt(f["best_p"])
    txt = [f"tail k ≥ xmin: {f['n_tail']:,} of {nz:,} diseases",
           f"best alternative: {ALTS[f['best']]}",
           f"   R = {f['best_R']:.1f}, p {bp if bp.startswith('<') else '= ' + bp}",
           f"not shown: {n0:,} disease{'s' if n0 != 1 else ''} with degree 0"]
    bx.text(0.0, -0.17, "\n".join(txt), transform=bx.transAxes, **TXT)
    panel_label(bx, lab)
    V["graphs"][g].update(n_nonzero=nz, log_bin_edges=[float(e) for e in edges],
                          density_integral=float(np.sum(hist * np.diff(edges))))
bxs[0].set_ylabel("p(k)")
bxs[0].set_ylim(1e-7, 0.3)
bxs[0].legend(handles=[Line2D([], [], color=BIN_COL, marker="o", ms=5.5, mec="white", lw=0.8, label="Log-binned p(k)"),
                       Line2D([], [], color="0.55", marker="o", ls="", ms=2.5, alpha=0.3, label="Unbinned p(k)"),
                       Line2D([], [], color=FIT_COL, lw=1.6, ls="--", label="Power-law fit, k ≥ xmin"),
                       Line2D([], [], color=FIT_COL, lw=0.9, ls=":", label="xmin")],
              loc="lower left", frameon=False, fontsize=9.5)
V.update(n_diseases=n, n_rare=int(rare.sum()), n_common=int(common.sum()), combined_mean_degree=float(dc.mean()),
         combined_min_degree=int(dc.min()), combined_max_degree=int(dc.max()),
         combined_avg_clustering=float(TG["avg_clustering"]))
for ext in ("png", "pdf"):
    fig.savefig(HERE / f"si_mechdist_degree_by_channel.{ext}", dpi=DPI, bbox_inches="tight", facecolor="white")
json.dump(V, open(HERE / "si_mechdist_degree_by_channel_values.json", "w"), indent=1, sort_keys=True)
print("si_mechdist_degree_by_channel:", {g: (round(V["graphs"][g]["alpha"], 4), V["graphs"][g]["xmin"]) for g in GRAPHS})
