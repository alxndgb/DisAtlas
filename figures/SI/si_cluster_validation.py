"""SI figure: extended cluster validation at K = 55 on the 13,021-disease index (#fig:SI_cluster_validation).
Panels: a mean silhouette of k-medoids on D = 1 - p over K across 100 subsamples of 80% of diseases (band 2.5-97.5%;
        dashed = best mean minus one SD; K = 55 marked); b per-disease mean distance to own vs nearest other cluster
        (all 13,021 diseases, K = 55); c agreement with the 19 MONDO body-system classes (observed vs shuffled labels and
        random pairs; adjusted mutual information).
Inputs: k_selection_boot.tsv.gz (select_k_bootstrap_lfl.py on the 13,021 matrix), within_between_K55_13021.tsv.gz and
        sifu_cluster_13021_values.json (sifu_cluster_13021.py). Same panel design as the manuscript figure.
Usage: si_cluster_validation_k55.py <k_selection_boot.tsv.gz> <within_between.tsv.gz> <values.json> <outdir>
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import sys, os, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from disatlas import figstyle_si as fs
T = DERIVED / "si_cluster"; BOOT, WBF, VJ, OUT = T / "k_selection_boot.tsv.gz", T / "within_between_K55_13021.tsv.gz", T / "sifu_cluster_13021_values.json", HERE; fs.apply_figure_style(); K_SEL = 55
B = pd.read_csv(BOOT, sep="\t"); B = B[B.kind == "boot"]
g = B.groupby("K").silhouette
S = pd.DataFrame(dict(mean=g.mean(), sd=g.std(), lo=g.quantile(.025), hi=g.quantile(.975), n=g.size())).reset_index()
assert (S.n == 100).all(), S.n.unique()
best = S.loc[S["mean"].idxmax()]; thr = best["mean"] - best["sd"]; k_one_se = int(S.loc[S["mean"] >= thr, "K"].min())
WB = pd.read_csv(WBF, sep="\t"); assert len(WB) == 13021 and WB.cluster.nunique() == K_SEL
V = json.load(open(VJ)); ext, intl = V["m3_external"], V["internal"]
med_a, med_b = float(np.median(WB.within)), float(np.median(WB.nearest_other))
assert abs(med_a - intl["median_within"]) < 1e-3 and abs(med_b - intl["median_nearest_other"]) < 1e-3
BLUE, GREY, ORANGE = "#0072B2", "0.55", "#E69F00"
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), gridspec_kw=dict(width_ratios=[1, 1, .95], wspace=.75))
ax = axs[0]
ax.fill_between(S.K, S.lo, S.hi, color=BLUE, alpha=.18, lw=0); ax.plot(S.K, S["mean"], "-o", color=BLUE, ms=2.5, lw=1)
ax.axhline(thr, color="0.4", ls="--", lw=.7); ax.text(81, thr + .013, "best mean − 1 SD", fontsize=6, color="0.35", ha="right")
ks = S[S.K == K_SEL].iloc[0]; ax.plot([K_SEL], [ks["mean"]], "o", ms=6, mfc="none", mec="black", mew=1)
ax.annotate(f"K = {K_SEL}", (K_SEL, ks["mean"]), xytext=(K_SEL + 2, ks["mean"] - .045), ha="center", fontsize=6,
            arrowprops=dict(arrowstyle="-", color="0.3", lw=.6))
ax.set_xlabel("Number of clusters K"); ax.set_ylabel("Mean silhouette"); ax.set_xticks(range(10, 81, 10))
ax.set_title("a  Silhouette over K", loc="left")
ax = axs[1]; bins = np.linspace(0, 1, 51)
ax.hist(WB.within, bins=bins, color=BLUE, alpha=.75, label="own cluster")
ax.hist(WB.nearest_other, bins=bins, color=ORANGE, alpha=.6, label="nearest other")
ym = ax.get_ylim()[1] * 1.45; ax.set_ylim(0, ym)
for v, col in [(med_a, BLUE), (med_b, "#9a5a1c")]: ax.axvline(v, ymax=.82, color=col, ls="--", lw=.8)
ax.text(med_a - .02, ym * .66, f"median\n{med_a:.2f}", ha="right", fontsize=6, color=BLUE)
ax.text(med_b + .02, ym * .66, f"median\n{med_b:.2f}", ha="left", fontsize=6, color="#9a5a1c")
ax.set_xlabel("Mean distance (1 − p)"); ax.set_ylabel("Diseases"); ax.set_xlim(0, 1)
ax.legend(loc="upper left", fontsize=6, ncol=2, columnspacing=.8, handlelength=1, frameon=True, framealpha=1, edgecolor="none", facecolor="white")
ax.set_title("b  Within vs between clusters", loc="left")
ax = axs[2]
rows = [("Same-cluster\npairs sharing\na body system", ext["pair_share_within_cluster"], ext["pair_share_within_cluster_null_mean"], ext["pair_share_random"]),
        ("Adjusted\nmutual\ninformation", ext["ami"], ext["ami_null_max"], None)]
for y, (lab, obs, null, rnd) in zip([1, 0], rows):
    ax.plot([null, obs], [y, y], color="0.8", lw=1.5, zorder=1)
    ax.plot([null], [y], "o", ms=5, mfc="white", mec=GREY, mew=1, zorder=2); ax.plot([obs], [y], "o", ms=5, color=BLUE, zorder=3)
    ax.text(obs, y + .14, f"{obs:.2f}", ha="center", va="bottom", fontsize=6, color=BLUE)
    ax.text(max(null, .02), y - .14, f"{null:.2f}" if null >= .01 else f"max {null:.3f}", ha="center" if null >= .01 else "left", va="top", fontsize=6, color="0.35")
    if rnd is not None: ax.plot([rnd], [y], "|", ms=9, color="0.2", mew=1, zorder=2)
ax.set_yticks([1, 0]); ax.set_yticklabels([r[0] for r in rows]); ax.set_ylim(-.7, 1.9); ax.set_xlim(-.03, .8)
ax.set_xlabel("Share of pairs / AMI")
h = [Line2D([], [], marker="o", ls="", color=BLUE, ms=4, label="observed"),
     Line2D([], [], marker="o", ls="", mfc="white", mec=GREY, ms=4, label="shuffled labels"),
     Line2D([], [], marker="|", ls="", color="0.2", ms=7, label="random pairs")]
ax.legend(handles=h, loc="upper right", fontsize=6, handletextpad=.2, borderaxespad=0)
ax.set_title("c  MONDO body systems", loc="left")
N = dict(png_bytes=fs.save(fig, OUT, "extended_cluster_validation"), overlaps=fs.overlap_check(fig), best_K=int(best.K),
         best_mean=best["mean"], best_sd=best["sd"], k_one_se=k_one_se, k55_mean=ks["mean"], k55_lo=ks.lo, k55_hi=ks.hi,
         within_1sd=S.loc[S["mean"] >= thr, "K"].tolist(), med_within=med_a, med_between=med_b)
json.dump(N, open(f"{OUT}/extended_cluster_validation_values.json", "w"), indent=1, default=float); print(N)
