"""Supplementary Figure SI_uncertainty: seed and sampling uncertainty on IxIDN-Matched.

Appendix section: Uncertainty Quantification (#sec:SI_uncertainty), Figure #fig:SI_uncertainty.
  a  SD of AUROC across five DisAtlas training seeds (blue) and across a disease-cluster bootstrap (B = 2,000, orange),
     on all 4,364 IxIDN-Matched pairs and on the 1,094 test-split pairs.
  b  AUROC and AUPRC of the five-seed DisAtlas ensemble with 95% disease-cluster bootstrap CIs; grey ticks = no skill.
Uses its own plot style (kept to match the manuscript image), not the SI house style.
Inputs: figures/SI/source_data/uncertainty_seed_vs_sample.tsv, figures/SI/source_data/uncertainty_disatlas_ci.tsv
Output: figures/SI/uncertainty_quantification.{png,pdf}
Run:    python figures/SI/si_uncertainty.py
"""
import os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _style import src, save  # noqa: E402

S = pd.read_csv(src("uncertainty_seed_vs_sample.tsv"), sep="\t").set_index(["stratum", "metric"])
C = pd.read_csv(src("uncertainty_disatlas_ci.tsv"), sep="\t").set_index(["stratum", "metric"])
ST = [("all_pairs", "All pairs\n(n = 4,364)"), ("test_split", "Test split\n(n = 1,094)")]
NOSKILL_AUPRC = {"all_pairs": 0.500, "test_split": 0.449}
BLUE, ORANGE, INK = "#3B6EA8", "#D9822B", "#222222"

# ---- value check: appendix #sec:SI_uncertainty ----
_sd = S.loc[("all_pairs", "auroc")]
assert (round(_sd.seed_sd, 4), round(_sd.cluster_bootstrap_sd, 4)) == (0.0015, 0.0203)
assert round(_sd.cluster_bootstrap_sd / _sd.seed_sd) == 13
assert round(C.loc[("all_pairs", "auroc"), "value"], 3) == 0.765 and round(C.loc[("test_split", "auroc"), "value"], 3) == 0.782
print("appendix values: OK")

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.3})

# a: where the uncertainty comes from
x = range(len(ST)); w = 0.36
seed = [S.loc[(s, "auroc"), "seed_sd"] for s, _ in ST]
samp = [S.loc[(s, "auroc"), "cluster_bootstrap_sd"] for s, _ in ST]
a.bar([i - w / 2 for i in x], seed, w, color=BLUE, label="Retraining with a new seed (5 seeds)")
a.bar([i + w / 2 for i in x], samp, w, color=ORANGE, label="Resampling the diseases evaluated")
for i in x:
    a.text(i + w / 2, samp[i] + 0.0006, f"{samp[i]/seed[i]:.0f}× larger", ha="center", va="bottom", fontsize=8, color=INK)
a.set_xticks(list(x)); a.set_xticklabels([l for _, l in ST], fontsize=8)
a.set_ylabel("SD of AUROC"); a.set_ylim(0, max(samp) * 1.45)
a.legend(frameon=False, fontsize=8, loc="upper left")
a.set_title("Test-set sampling, not the seed,\ndrives the uncertainty", loc="left", fontsize=9)

# b: DisAtlas performance with 95% CI
mets = [("auroc", "AUROC"), ("auprc", "AUPRC")]
for j, (s, lab) in enumerate(ST):
    for k, (m, mlab) in enumerate(mets):
        y = j * 2.6 + k
        r = C.loc[(s, m)]
        b.plot([r.cluster_ci_lo, r.cluster_ci_hi], [y, y], color=BLUE, lw=2, solid_capstyle="butt")
        b.plot(r.value, y, "o", color=BLUE, ms=6)
        b.text(r.cluster_ci_hi + 0.01, y, f"{r.value:.2f} ({r.cluster_ci_lo:.2f}–{r.cluster_ci_hi:.2f})", va="center", fontsize=8)
        ns = 0.5 if m == "auroc" else NOSKILL_AUPRC[s]
        b.plot([ns, ns], [y - 0.35, y + 0.35], color="0.5", lw=1.2)
yt = [j * 2.6 + k for j in range(len(ST)) for k in range(len(mets))]
b.set_yticks(yt); b.set_yticklabels([m for _ in ST for _, m in mets], fontsize=8)
for j, (s, lab) in enumerate(ST):
    b.text(0.41, j * 2.6 - 0.75, lab.replace("\n", " "), ha="left", va="center", fontsize=8, fontweight="bold")
b.text(0.505, -0.45, "no skill", ha="left", va="center", fontsize=7, color="0.4")
b.set_xlim(0.40, 1.0); b.set_ylim(len(ST) * 2.6 - 1.3, -1.1)
b.set_xlabel("Score")
b.set_title("DisAtlas on IxIDN-Matched\n(95% CI, disease-cluster bootstrap)", loc="left", fontsize=9)
b.tick_params(axis="y", length=0)

for ax, L in [(a, "A"), (b, "B")]:
    ax.text(-0.02, 1.08, L, transform=ax.transAxes, fontweight="bold", fontsize=11, ha="right")
print(save(fig, "uncertainty_quantification", bbox_inches="tight"))
