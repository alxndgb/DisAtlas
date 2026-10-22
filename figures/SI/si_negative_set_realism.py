"""Supplementary Figure SI_negrealism: negative-set realism and hard-negative stress test.

Appendix section: Negative-set realism (#sec:SI_negrealism), Figure #fig:SI_negrealism.
(a) ORDO max-LCA distance for random disease pairs and for IxIDN-Matched negatives / positives;
(b) DisAtlas AUROC of all positives against each hard-negative subset, and within shared-class strata.
The DisAtlas scores are unchanged; only the negative pairs scored are filtered.
Inputs:  figures/SI/source_data/negative_realism_ordo_distance.tsv.gz
         figures/SI/source_data/negative_realism_hard_negative_auroc.tsv
Output:  figures/SI/negative_set_realism.{png,pdf}
Run:     python figures/SI/si_negative_set_realism.py
"""
import os, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import save, SOURCE_DATA as SRC  # noqa: E402

L = pd.read_csv(os.path.join(SRC, "negative_realism_ordo_distance.tsv.gz"), sep="\t")
T = pd.read_csv(os.path.join(SRC, "negative_realism_hard_negative_auroc.tsv"), sep="\t")

SETS = [("random_disease_pairs", "Random disease\npairs", "#d9d9d9"),
        ("ixidn_matched_negatives", "IxIDN-Matched\nnegatives", "#d9822b"),
        ("ixidn_matched_positives", "IxIDN-Matched\npositives", "#3b6ea8")]
SUB = [("All negatives (reference)", "All negatives"),
       ("Share an organ system with each other", "Same organ system"),
       ("Share an ORDO clinical area", "Same ORDO clinical area"),
       ("MONDO distance <= positive median (ontology-close)", "Negatives close in MONDO"),
       ("ORDO distance <= positive median (ontology-close)", "Negatives close in ORDO"),
       ("GO Resnik gene-BMA >= positive median (function-similar)", "Negatives similar by GO function"),
       ("HPO Resnik-BMA >= positive median (phenotype-similar)", "Negatives similar by HPO phenotype"),
       ("Close or similar on >= 2 of the 4 axes above", "Negatives close or similar on ≥2 axes"),
       ("Both classes share an organ system", "Within stratum: organ system"),
       ("Both classes share an ORDO clinical area", "Within stratum: ORDO area")]
ALL, TEST = "All 4,364 pairs", "Test-split pairs"

# ---- appendix values (#sec:SI_negrealism) -----------------------------------------------------
med = L.groupby("set").ordo_dist.median()
assert med["ixidn_matched_negatives"] == 8 and med["ixidn_matched_positives"] == 6
ra = T[T.scope == ALL].set_index("subset")
assert round(ra.loc["Share an organ system with each other", "pct_neg_kept"]) == 37
for s in ["Share an organ system with each other", "MONDO distance <= positive median (ontology-close)"]:
    assert f"{ra.loc[s, 'auroc_disatlas']:.2f}" == "0.72", s
for s in ["ORDO distance <= positive median (ontology-close)", "GO Resnik gene-BMA >= positive median (function-similar)"]:
    assert f"{ra.loc[s, 'auroc_disatlas']:.2f}" in ("0.66", "0.67"), s
h = ra.loc["HPO Resnik-BMA >= positive median (phenotype-similar)"]
assert round(h.pct_neg_kept) == 7 and f"{h.auroc_disatlas:.2f}" == "0.54"
o = ra.loc["Both classes share an organ system"]
assert (f"{o.auroc_disatlas:.3f}", f"{o.auroc_hpo_resnik_bma:.3f}", f"{o.auroc_go_resnik_gene_bma:.3f}") == ("0.780", "0.585", "0.641")
print("appendix values: OK")

# ---- figure -----------------------------------------------------------------------------------
plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False})
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.4), gridspec_kw=dict(width_ratios=[1, 1.45], wspace=1.25))
data = [L.loc[L.set == s, "ordo_dist"].values for s, _, _ in SETS]
bp = a.boxplot(data, vert=False, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color="k"))
for p, (_, _, col) in zip(bp["boxes"], SETS):
    p.set_facecolor(col)
a.set_yticks(range(1, len(SETS) + 1)); a.set_yticklabels([f"{l}\n(n = {len(d):,})" for (_, l, _), d in zip(SETS, data)]); a.invert_yaxis()
a.set_xlabel("ORDO distance (max-LCA)"); a.set_title("a", loc="left", fontweight="bold", x=-0.55)
y = np.arange(len(SUB))
for sc, col, off, lab in [(ALL, "k", -0.17, "All pairs"), (TEST, "#d9822b", 0.17, "Test-split pairs")]:
    r = T[T.scope == sc].set_index("subset").loc[[s for s, _ in SUB]]
    b.errorbar(r.auroc_disatlas, y + off, xerr=[r.auroc_disatlas - r.ci_lo, r.ci_hi - r.auroc_disatlas], fmt="o", ms=3.5,
               color=col, capsize=1.5, lw=0.9, label=lab)
rs = ra.loc[[s for s, _ in SUB]]
for i, (s, _) in enumerate(SUB):
    t = (f"{int(rs.n_pos.iloc[i]):,} / {int(rs.n_neg.iloc[i]):,} pairs" if s.startswith("Both")
         else f"{int(rs.n_neg.iloc[i]):,} neg ({rs.pct_neg_kept.iloc[i]:.0f}%)")
    b.text(1.005, i, t, transform=b.get_yaxis_transform(), fontsize=6, va="center")
b.axvline(0.5, color="0.6", ls="--", lw=0.8); b.axhline(7.5, color="0.85", lw=0.8)
b.set_yticks(y); b.set_yticklabels([l for _, l in SUB]); b.invert_yaxis(); b.set_xlim(0.28, 0.92)
b.set_xlabel("AUROC, DisAtlas (95% CI)"); b.legend(frameon=False, loc="lower left", fontsize=6.5)
b.set_title("b", loc="left", fontweight="bold", x=-0.85)
print(save(fig, "negative_set_realism", bbox_inches="tight"))
