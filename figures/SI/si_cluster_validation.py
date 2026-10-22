"""Supplementary Figure SI_cluster_validation: extended cluster validation.

Appendix section: Cluster validation (#fig:SI_cluster_validation).
  a  mean silhouette of k-medoids clustering on the DisAtlas distance (1 - p) over K, across 100 subsamples
     of 80% of diseases (band: 2.5-97.5% range); grey shading = plateau K 45-80; dashed line = best mean
     minus one SD (one-SE rule), K = 50 marked.
  b  per-disease mean distance to the own cluster vs to the nearest other cluster (the two silhouette terms),
     all 13,019 clustered diseases at K = 50; dashed lines = medians.
  c  concordance with the 19 MONDO body-system classes: share of same-cluster pairs sharing a class
     (observed vs size-preserving label-permutation null, and vs random pairs), and adjusted mutual
     information (observed vs shuffled-label null maximum).
Inputs: figures/SI/source_data/cluster_k_selection_bootstrap.tsv.gz (a),
        figures/SI/source_data/cluster_within_between_k50.tsv.gz (b),
        figures/SI/source_data/cluster_mondo_concordance.json (c, and the panel-b medians check)
Output: figures/SI/extended_cluster_validation.{png,pdf}
Run:    python figures/SI/si_cluster_validation.py
"""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _style import apply_style, panel_label, save, src, ANNOT_SIZE, LEGEND_SIZE

K_SEL = 50
PLATEAU = (45, 80)
BLUE, GREY, ORANGE = "#1f4e79", "0.55", "#d9822b"


def main():
    apply_style()
    B = pd.read_csv(src("cluster_k_selection_bootstrap.tsv.gz"), sep="\t")
    B = B[B.kind == "boot"]
    g = B.groupby("K").silhouette
    S = pd.DataFrame(dict(mean=g.mean(), sd=g.std(), lo=g.quantile(0.025), hi=g.quantile(0.975),
                          n=g.size())).reset_index()
    assert (S.n == 100).all(), S.n.unique()
    best = S.loc[S["mean"].idxmax()]
    thr = best["mean"] - best["sd"]
    k_one_se = int(S.loc[S["mean"] >= thr, "K"].min())
    assert k_one_se == K_SEL, (k_one_se, best.K)
    WB = pd.read_csv(src("cluster_within_between_k50.tsv.gz"), sep="\t")
    J = json.load(open(src("cluster_mondo_concordance.json")))
    ext, intl = J["mondo_body_system_concordance"], J["within_between_medians"]
    assert len(WB) == 13019 and ext["n_categories"] == 19
    med_a, med_b = float(np.median(WB.within)), float(np.median(WB.nearest_other))
    assert abs(med_a - intl["median_within_distance"]) < 2e-3 and abs(med_b - intl["median_nearest_other_cluster_distance"]) < 2e-3

    fig, axs = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw=dict(width_ratios=[1.0, 1.0, 0.95], wspace=0.55))

    # a: silhouette over K
    ax = axs[0]
    ax.axvspan(*PLATEAU, color="0.92", zorder=0, lw=0)
    ax.fill_between(S.K, S.lo, S.hi, color=BLUE, alpha=0.18, lw=0)
    ax.plot(S.K, S["mean"], "-o", color=BLUE, ms=4.5, lw=1.8)
    ax.axhline(thr, color="0.4", ls="--", lw=1)
    ksel = S[S.K == K_SEL].iloc[0]
    ax.plot([K_SEL], [ksel["mean"]], "o", ms=10, mfc="none", mec="black", mew=1.6)
    ax.annotate(f"K = {K_SEL}", (K_SEL, ksel["mean"]), xytext=(K_SEL - 4, ksel["mean"] - 0.06),
                ha="center", fontsize=ANNOT_SIZE + 1, arrowprops=dict(arrowstyle="-", color="0.3", lw=0.8))
    ax.text(62.5, S.lo.min() + 0.005, "plateau", ha="center", va="bottom", fontsize=ANNOT_SIZE, color="0.35")
    ax.text(11, thr + 0.004, "best mean − 1 SD", ha="left", va="bottom", fontsize=ANNOT_SIZE, color="0.35")
    ax.set_xlabel("Number of clusters K"); ax.set_ylabel("Mean silhouette")
    ax.set_xticks(range(10, 81, 10)); ax.set_xlim(7, 83)
    ax.set_title("Silhouette over K (100 subsamples)", loc="left")
    panel_label(ax, "a")

    # b: within vs nearest-other cluster distance
    ax = axs[1]
    bins = np.linspace(0, 1, 51)
    ax.hist(WB.within, bins=bins, color=BLUE, alpha=0.75, label="own cluster")
    ax.hist(WB.nearest_other, bins=bins, color=ORANGE, alpha=0.6, label="nearest other cluster")
    for v, col in [(med_a, BLUE), (med_b, ORANGE)]:
        ax.axvline(v, color=col, ls="--", lw=1.2)
    ymax = ax.get_ylim()[1] * 1.32; ax.set_ylim(0, ymax)
    ax.text(med_a - 0.015, ymax * 0.79, f"median {med_a:.2f}", ha="right", fontsize=ANNOT_SIZE, color=BLUE)
    ax.text(med_b + 0.015, ymax * 0.79, f"median {med_b:.2f}", ha="left", fontsize=ANNOT_SIZE, color="#9a5a1c")
    ax.set_xlabel(f"Mean distance (1 − p) per disease (n = {len(WB):,})"); ax.set_ylabel("Diseases")
    ax.set_xlim(0, 1)
    ax.legend(loc="upper left", frameon=False, ncol=2, title="mean distance to", title_fontsize=LEGEND_SIZE, alignment="left", columnspacing=1.2)
    ax.set_title(f"Distance within and between clusters", loc="left")
    panel_label(ax, "b")

    # c: MONDO concordance
    ax = axs[2]
    rows = [("Same-cluster\npairs sharing\na body-system\nclass", ext["pair_share_within_cluster"],
             ext["pair_share_within_cluster_null_mean"], ext["pair_share_within_cluster_null_95"], ext["pair_share_random"]),
            ("Adjusted\nmutual\ninformation", ext["ami"], ext["ami_null_max"], None, None)]
    ys = [1, 0]
    for y, (lab_, obs, null, ci, rnd) in zip(ys, rows):
        ax.plot([null, obs], [y, y], color="0.8", lw=2, zorder=1)
        ax.plot([null], [y], "o", ms=9, mfc="white", mec=GREY, mew=1.6, zorder=2)
        ax.plot([obs], [y], "o", ms=9, color=BLUE, zorder=3)
        ax.text(obs, y + 0.13, f"{obs:.2f}", ha="center", va="bottom", fontsize=ANNOT_SIZE + 1, color=BLUE)
        nl = f"{null:.2f}" if null >= 0.01 else f"max {null:.4f}"
        ax.text(max(null, 0.03), y - 0.13, nl, ha="center" if null >= 0.01 else "left", va="top", fontsize=ANNOT_SIZE + 1, color="0.35")
        if rnd is not None:
            ax.plot([rnd], [y], "|", ms=14, color="0.2", mew=1.4, zorder=2)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows])
    ax.set_ylim(-0.6, 1.6); ax.set_xlim(-0.03, 0.8)
    ax.set_xlabel("Share of pairs / AMI")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=BLUE, ms=8, label="observed clusters"),
         Line2D([], [], marker="o", ls="", mfc="white", mec=GREY, mew=1.5, ms=8, label="shuffled cluster labels"),
         Line2D([], [], marker="|", ls="", color="0.2", ms=12, mew=1.4, label="random disease pairs")]
    ax.legend(handles=h, loc="center right", bbox_to_anchor=(1.0, 0.5), frameon=False)
    ax.set_title("Agreement with MONDO body systems", loc="left")
    panel_label(ax, "c", xy=(-0.36, 1.03))

    print(save(fig, "extended_cluster_validation", bbox_inches="tight"))
    print(dict(best_K=int(best.K), best_mean=round(best["mean"], 4), sd=round(best["sd"], 4), k_one_se=k_one_se,
               k50_mean=round(ksel["mean"], 4), med_within=round(med_a, 3), med_between=round(med_b, 3),
               share_obs=round(ext["pair_share_within_cluster"], 3), share_null=round(ext["pair_share_within_cluster_null_mean"], 3),
               share_rand=round(ext["pair_share_random"], 3), ami=round(ext["ami"], 3), ami_null_max=round(ext["ami_null_max"], 4),
               plateau_within_1sd=S.loc[S["mean"] >= thr, "K"].tolist()))


if __name__ == "__main__":
    main()
