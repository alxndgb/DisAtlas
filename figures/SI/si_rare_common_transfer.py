"""Supplementary Figure SI_rare_common_stratified: rare/common transfer.

Appendix section: Rare-common stratified evaluation (#fig:SI_rare_common_transfer).
  a  test-split pairs: AUROC by pair type (rare-rare, rare-common, common-common, all) for DisAtlas (95% CI)
     and the HPO / GO similarity baselines; positives/negatives printed under each stratum. A stratum with
     fewer than MIN_N pairs in either class gets no AUROC ("n too small").
  b  anchored analysis, all pairs: four anchor/partner groups under one per-anchor rule (an anchor qualifies
     with >= 1 positive and >= 1 negative partner). A plotted group with sparse negatives (< SPARSE_NEG
     negative pairs) is drawn hollow.
Inputs: figures/SI/source_data/rare_common_strata_auroc.tsv, figures/SI/source_data/rare_common_anchored_groups.tsv
Output: figures/SI/SI_rare_common_transfer.{png,pdf} and SI_rare_common_transfer_values.tsv
Run:    python figures/SI/si_rare_common_transfer.py
"""
import os, sys, json
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _style import apply_style, panel_label, save, src, ANNOT_SIZE

STEM = "SI_rare_common_transfer"
T = pd.read_csv(src("rare_common_strata_auroc.tsv"), sep="\t").rename(columns={
    "hpo_resnik_bma_auroc": "HPO Resnik-BMA | auroc", "go_resnik_genebma_auroc": "GO Resnik gene-BMA | auroc"})
T["scope"] = T.scope.map({"all": "All 4364 pairs", "test": "Held-out test-test pairs"})

DA_COL = "#d9822b"
BL = [("HPO Resnik-BMA", "HPO similarity", "#7fa6cf", "s"), ("GO Resnik gene-BMA", "GO similarity", "#4c9a6a", "^")]
ST = [("rare-rare", "Rare–\nrare"), ("rare-common", "Rare–\ncommon"), ("common-common", "Common–\ncommon"), ("all", "All\npairs")]
MIN_N = 5  # minimum pairs per class to plot an AUROC
SC = [("All 4364 pairs", "All IxIDN-Matched pairs"), ("Held-out test-test pairs", "Test-split pairs")]
DIR = [("rare_anchor_rank_common_partners | all", "Rare\nanchor,\ncommon\npartners"),
       ("common_anchor_rank_rare_partners | all", "Common\nanchor,\nrare\npartners"),
       ("rare_anchor_rank_rare_partners | all", "Rare\nanchor,\nrare\npartners"),
       ("common_anchor_rank_common_partners | all", "Common\nanchor,\ncommon\npartners")]
MIN_ANCH = 5     # minimum qualifying anchors to plot an anchored group
SPARSE_NEG = 50  # plotted groups with fewer negative pairs than this are drawn hollow
A4 = pd.read_csv(src("rare_common_anchored_groups.tsv"), sep="\t").set_index("group")

# n annotations below: positives/negatives per stratum
# ---- value check: every plotted number equals the appendix value (3 d.p.) ----
EXPECTED = {
    ("All 4364 pairs", "rare-rare"): (0.722, 0.672, 0.772, 0.575, 0.601, 1284, 1168),
    ("All 4364 pairs", "rare-common"): (0.758, 0.668, 0.835, 0.534, 0.609, 253, 988),
    ("All 4364 pairs", "all"): (0.765, 0.729, 0.799, 0.544, 0.645, 2182, 2182),
    ("Held-out test-test pairs", "rare-rare"): (0.743, 0.656, 0.817, 0.610, 0.605, 301, 321),
    ("Held-out test-test pairs", "rare-common"): (0.799, 0.656, 0.920, 0.574, 0.505, 68, 279),
    ("Held-out test-test pairs", "all"): (0.782, 0.721, 0.841, 0.576, 0.625, 491, 603),
}
for (sc, st), e in EXPECTED.items():
    r = T[(T.scope == sc) & (T.stratum == st)].iloc[0]
    got = (r.auroc, r.ci_lo, r.ci_hi, r["HPO Resnik-BMA | auroc"], r["GO Resnik gene-BMA | auroc"], r.n_pos, r.n_neg)
    assert all(abs(g - x) < 6e-4 for g, x in zip(got, e)), (sc, st, got, e)
# test-split common-common stratum: counts only, no AUROC (too few negatives)
cc = T[(T.scope == "Held-out test-test pairs") & (T.stratum == "common-common")].iloc[0]
assert (int(cc.n_pos), int(cc.n_neg)) == (122, 3) and pd.isna(cc.auroc), cc[["n_pos", "n_neg", "auroc"]].tolist()
EXP_ANCH = {DIR[0][0]: (0.833, 0.717, 0.935, 0.516, 0.627, 31, 120),
            DIR[1][0]: (0.767, 0.706, 0.822, 0.557, 0.627, 104, 934)}
for k, e in EXP_ANCH.items():
    for d in (A4.loc[k.split(" | ")[0]],):
        got = (d["mean_auroc"], d["ci_lo"], d["ci_hi"], d["mean_auroc_hpo"], d["mean_auroc_go"], d["n_anchors"], d["n_pairs"])
        assert all(abs(g - x) < 6e-4 for g, x in zip(got, e)), (k, got, e)

apply_style()
fig, ax = plt.subplots(1, 2, figsize=(14.0, 6.4), gridspec_kw=dict(width_ratios=[1, 1], wspace=0.22))
MS, CAP, LW, S = 10, 4, 2.0, 60


def da_point(a, x, y, lo, hi, hollow=False):
    a.errorbar(x, y, yerr=[[y - lo], [hi - y]], fmt="o", color=DA_COL, ms=MS, capsize=CAP, lw=LW, zorder=3,
               mfc="white" if hollow else DA_COL, mec=DA_COL, mew=LW)


def bl_point(a, x, y, col, m, hollow=False):
    a.scatter(x, y, marker=m, s=S, zorder=2, facecolors="white" if hollow else col, edgecolors=col,
              linewidths=1.6 if hollow else 0)


PLOTTED = []
a = ax[0]
Rt = T[T.scope == "Held-out test-test pairs"].set_index("stratum")
for i, (st, lab) in enumerate(ST):
    r = Rt.loc[st]
    if min(r.n_pos, r.n_neg) < MIN_N:
        a.text(i + 0.1, 0.62, "n too small", ha="center", va="center", fontsize=ANNOT_SIZE, color="0.35", style="italic")
        PLOTTED.append(("a", lab.replace("\n", " ").replace("– ", "–"), None, None, None, None, None, int(r.n_pos), int(r.n_neg)))
    else:
        da_point(a, i, r.auroc, r.ci_lo, r.ci_hi)
        for k, (c_, _, col, m) in enumerate(BL):
            bl_point(a, i + 0.14 + 0.08 * k, r[f"{c_} | auroc"], col, m)
        PLOTTED.append(("a", lab.replace("\n", " ").replace("– ", "–"), r.auroc, r.ci_lo, r.ci_hi, r["HPO Resnik-BMA | auroc"],
                        r["GO Resnik gene-BMA | auroc"], int(r.n_pos), int(r.n_neg)))
    a.text(i + 0.1, 0.355, f"{int(r.n_pos)}/{int(r.n_neg)}", ha="center", va="bottom", fontsize=ANNOT_SIZE, color="0.35")
a.set_xticks(range(len(ST))); a.set_xticklabels([l.replace("\n", "\n", 1) for _, l in ST]); a.set_xlim(-0.45, len(ST) - 0.35)
a.set_title("Test-split pairs", loc="left")
a.set_ylabel("AUROC (DisAtlas: 95% CI)")

c = ax[1]
HOLLOW_B = []
for i, (k, lab) in enumerate(DIR):
    d = A4.loc[k.split(" | ")[0]]
    nm = lab.replace("\n", " ")
    if d.n_anchors < MIN_ANCH:
        c.text(i, 0.62, "too few\nnegatives", ha="center", va="center", fontsize=ANNOT_SIZE, color="0.35", style="italic")
        PLOTTED.append(("b", nm, None, None, None, None, None, int(d.n_anchors), int(d.n_pairs)))
    else:
        hol = d.n_neg < SPARSE_NEG
        if hol: HOLLOW_B.append(int(d.n_neg))
        da_point(c, i, d.mean_auroc, d.ci_lo, d.ci_hi, hollow=hol)
        for k2, (key, col, m) in enumerate([("mean_auroc_hpo", BL[0][2], "s"), ("mean_auroc_go", BL[1][2], "^")]):
            bl_point(c, i + 0.14 + 0.08 * k2, d[key], col, m, hollow=hol)
        PLOTTED.append(("b", nm, d.mean_auroc, d.ci_lo, d.ci_hi, d.mean_auroc_hpo, d.mean_auroc_go,
                        int(d.n_anchors), int(d.n_pairs)))
    c.text(i + 0.07, 0.355, f"{int(d.n_anchors)} anchors\n{int(d.n_pairs)} pairs", ha="center", va="bottom",
           fontsize=ANNOT_SIZE, color="0.35")
c.set_xticks(range(len(DIR))); c.set_xticklabels([l for _, l in DIR]); c.set_xlim(-0.45, len(DIR) - 0.35)
c.set_ylabel("Mean per-anchor AUROC"); c.set_title("Anchored, all pairs", loc="left")

for a_, letter in zip(ax, "ab"):
    a_.axhline(0.5, color="0.8", ls="--", lw=1.0, zorder=1)
    a_.set_ylim(0.33, 1.0)
    panel_label(a_, letter)

h = [Line2D([], [], color=DA_COL, marker="o", ms=MS, ls="", label="DisAtlas")] + \
    [Line2D([], [], color=col, marker=m, ms=8, ls="", label=lab) for _, lab, col, m in BL]
assert len(HOLLOW_B) <= 1
if HOLLOW_B:
    h.append(Line2D([], [], color="0.35", marker="o", mfc="white", mec="0.35", mew=1.6, ms=MS, ls="",
                    label=f"hollow = few negatives (n = {HOLLOW_B[0]})"))
fig.legend(handles=h, frameon=False, ncol=len(h), loc="upper center", bbox_to_anchor=(0.5, -0.04))
print(save(fig, STEM, bbox_inches="tight"))
TAB = os.path.join(HERE, STEM + "_values.tsv")
pd.DataFrame(PLOTTED, columns=["panel", "stratum", "disatlas_auroc", "ci_lo", "ci_hi", "hpo_auroc", "go_auroc",
                               "n_pos_or_anchors", "n_neg_or_pairs"]).round(3).to_csv(TAB, sep="\t", index=False)
print(TAB)
