"""Supplementary Figure SI_placement_confidence: per-disease placement confidence, exclusion sensitivity.

Appendix section: Per-disease placement confidence (#sec:SI_placement_confidence).
Panels: a, all 4,364 IxIDN-Matched pairs; b, the 1,094-pair held-out test split. Each row is the change in
AUROC when one group of lower-confidence diseases is excluded (point, 95% disease-cluster bootstrap CI), against
the 95% range when the same number of diseases is excluded at random (grey bar). Nothing is retrained.

Input:  figures/SI/source_data/placement_exclusion_sensitivity.tsv
Output: figures/SI/SI_placement_confidence.{png,pdf}
Run:    python figures/SI/si_placement_confidence.py
"""
import os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import apply_style, panel_label, save, src, ANNOT_SIZE

BLUE = "#1f4e79"
PANELS = [("all", "All IxIDN-Matched pairs (n = 4,364)"),
          ("test", "Held-out test split (n = 1,094)")]
ROWS = [("unstable", "Unstable placement"),
        ("ill_defined", "Ill-defined disease term"),
        ("unstable_or_ill_defined", "Unstable or ill-defined"),
        ("lowest_confidence_tier", "Lowest-confidence tier")]


def check(T):
    """Appendix text: 19-46% of pairs removed, AUROC change -0.012 to +0.016, every CI spans 0, all within null."""
    assert (round(100 - T.pct_pairs_kept.max()), round(100 - T.pct_pairs_kept.min())) == (19, 46)
    assert (round(T.delta.min(), 3), round(T.delta.max(), 3)) == (-0.012, 0.016)
    assert ((T.delta_lo <= 0) & (T.delta_hi >= 0)).all()
    assert ((T.delta >= T.null_delta_lo) & (T.delta <= T.null_delta_hi)).all()
    assert set(T.n_pairs_total[T.benchmark == "all"]) == {4364} and set(T.n_pairs_total[T.benchmark == "test"]) == {1094}


def main():
    apply_style()
    T = pd.read_csv(src("placement_exclusion_sensitivity.tsv"), sep="\t")
    check(T)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True, gridspec_kw=dict(wspace=0.32))
    for k, (ax, (bench, title)) in enumerate(zip(axes, PANELS)):
        S = T[T.benchmark == bench].set_index("exclusion")
        for i, (key, _) in enumerate(ROWS):
            r = S.loc[key]; y = len(ROWS) - 1 - i
            ax.barh(y, r.null_delta_hi - r.null_delta_lo, left=r.null_delta_lo, height=0.55, color="0.88", zorder=1)
            ax.errorbar(r.delta, y, xerr=[[r.delta - r.delta_lo], [r.delta_hi - r.delta]], fmt="o", color=BLUE,
                        ms=7, capsize=3, lw=1.6, zorder=3)
            ax.text(0.085, y, f"{r.pct_pairs_kept:.0f}%", va="center", ha="right", fontsize=ANNOT_SIZE, color="0.35")
        ax.axvline(0, color="0.3", ls="--", lw=1)
        ax.set_xlim(-0.06, 0.09); ax.set_xticks([-0.04, 0, 0.04]); ax.set_xticklabels(["\u22120.04", "0", "+0.04"])
        ax.set_ylim(-0.6, len(ROWS) - 0.4)
        ax.set_yticks(range(len(ROWS))); ax.set_yticklabels([r[1] for r in ROWS][::-1])
        ax.set_xlabel("Change in AUROC after exclusion")
        ax.set_title(title, loc="left")
        ax.text(0.085, len(ROWS) - 0.45, "pairs kept", ha="right", va="bottom", fontsize=ANNOT_SIZE, color="0.35")
        panel_label(ax, "ab"[k], xy=(-0.10 if k else -0.42, 1.03))
    h = [Line2D([], [], marker="o", color=BLUE, ms=7, label="observed change (95% CI)"),
         Patch(color="0.88", label="random exclusion of the same number of diseases (95% range)")]
    fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.55, -0.07))
    print(save(fig, "SI_placement_confidence", bbox_inches="tight"))


if __name__ == "__main__":
    main()
