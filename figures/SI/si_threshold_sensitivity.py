"""Supplementary Figure SI_threshold_sensitivity: robustness to the GO/HPO distance cutoff.

Appendix section: Robustness to the cutoff value (#sec:SI_threshold_sweep).
  a  AUROC on all 4,364 IxIDN-Matched pairs (solid) and on the 1,094 test-split pairs (dashed) for networks rebuilt
     at cutoffs 0.10-0.30 with DisAtlas retrained from scratch (three-seed ensembles); open marker at 0.10 = the
     4,104 pairs scorable at that cutoff; star = the five-seed DisAtlas ensemble at the chosen cutoff 0.20.
  b  diseases retained in the network at each cutoff.
Rendered with matplotlib defaults (this figure predates the house style; kept to match the manuscript image).
Inputs: figures/SI/source_data/threshold_sweep.tsv, figures/SI/source_data/threshold_reference_disatlas_ensemble.tsv
Output: figures/SI/threshold_sensitivity.{png,pdf}
Run:    python figures/SI/si_threshold_sensitivity.py
"""
import os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import save, src

t = pd.read_csv(src("threshold_sweep.tsv"), sep="\t")
REF = pd.read_csv(src("threshold_reference_disatlas_ensemble.tsv"), sep="\t").set_index("pairs").auroc
ref_all, ref_tt = float(REF["all"]), float(REF["test"])

# appendix text values (3 d.p.): 0.15-0.30 retain 12,848-13,051 diseases, 13-31% edges change, AUROC 0.753-0.764
# (test 0.766-0.784); 0.10 loses 0.035 on the pairs every cutoff scores; sweep 0.764 at 0.20 vs ensemble 0.765.
s = t[t.cutoff >= 0.15]
assert (s.diseases.min(), s.diseases.max()) == (12848, 13051)
assert (round(100 * (1 - s["edge_jaccard_vs_cutoff_0.20"][s.cutoff != 0.2]).min()),
        round(100 * (1 - s["edge_jaccard_vs_cutoff_0.20"]).max())) == (13, 31)
assert (round(s.auroc_ens.min(), 3), round(s.auroc_ens.max(), 3)) == (0.753, 0.764)
assert (round(s.auroc_ens_test_test.min(), 3), round(s.auroc_ens_test_test.max(), 3)) == (0.766, 0.784)
c = t.set_index("cutoff").auroc_ens_common_pairs
assert round(c[0.2] - c[0.1], 3) == 0.035 and round(ref_all, 3) == 0.765 and round(ref_tt, 3) == 0.782
dark, light, grey, shade, acc = "#1f4e9c", "#6a9fd6", "#9a9a9a", "#f7ead8", "#c0392b"
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1.5, 1]))
for ax in (a, b):
    ax.axvspan(0.175, 0.225, color=shade, alpha=.6, lw=0); ax.spines[["top", "right"]].set_visible(False)
    ax.set_xticks(t.cutoff); ax.set_xlabel("GO/HPO distance cutoff", fontsize=8); ax.tick_params(labelsize=6.5)
x = t.cutoff
a.plot(x, t.auroc_ens_test_test, "--s", color=light, ms=5, lw=1.4)
a.plot(x, t.auroc_ens, "-o", color=dark, ms=5, lw=1.8)
a.plot([0.10], [t.auroc_ens.iloc[0]], "o", mfc="white", mec=dark, ms=5, mew=1.4)
a.plot([0.20], [ref_all], "*", color=acc, ms=10, zorder=5)
a.annotate(f"DisAtlas\n(5 seeds) {ref_all:.3f}", (0.20, ref_all), (0.222, 0.792), fontsize=6.5, color=acc,
           arrowprops=dict(arrowstyle="-", color=acc, lw=.6))
a.text(0.307, t.auroc_ens_test_test.iloc[-1] + 0.002, "test-split\npairs", color=light, fontsize=7, va="center")
a.text(0.307, t.auroc_ens.iloc[-1] - 0.002, "all 4,364\npairs", color=dark, fontsize=7, va="center")
a.annotate("0.10: 260 pairs\nunscorable", (0.10, t.auroc_ens.iloc[0]), (0.106, 0.703), fontsize=6, color="#666666",
           arrowprops=dict(arrowstyle="-", color="#666666", lw=.6))
a.text(0.20, 0.806, "chosen (0.20)", color="#d9822b", fontsize=6.5, ha="center")
a.set_ylim(0.695, 0.81); a.set_xlim(0.085, 0.36); a.set_ylabel("IxIDN-Matched AUROC", fontsize=8)
a.set_title("Retrained at each cutoff (3 seeds)", fontsize=8, loc="left")
b.plot(x, t.diseases / 1000, "-o", color="#555555", ms=5, lw=1.6); b.set_ylim(10.6, 13.4)
b.set_ylabel("Diseases in network (thousands)", fontsize=8); b.set_title("Only 0.10 drops diseases", fontsize=8, loc="left")
i20 = int((t.cutoff - 0.20).abs().idxmin())
b.annotate(f"−{int(t.diseases.iloc[i20]-t.diseases.iloc[0]):,} diseases", (0.10, t.diseases.iloc[0]/1000), (0.135, 11.2), fontsize=6.5,
           color="#555555", arrowprops=dict(arrowstyle="-", color="#555555", lw=.6))
for ax, l in ((a, "a"), (b, "b")): ax.text(-0.14, 1.04, l, transform=ax.transAxes, fontsize=9, fontweight="bold")
fig.tight_layout(); print(save(fig, "threshold_sensitivity"), ref_all, ref_tt)
