"""SI figure: benchmark AUROC after excluding diseases by placement-confidence flag (#fig:SI_placement_confidence).
All 20 rows of the exclusion table are drawn (4 benchmarks x 5 rules); none is selected.
Input: disatlas/v3_revised_manuscript/si_figures/ixidn_si/v3.2.0/dq_exclusion_sensitivity.tsv
       (revision/disease_data_quality_audit/code/s_exclusion_sensitivity.py; disease-cluster bootstrap CIs).
Panel a: change in AUROC after exclusion (kept minus full), 95% CI; grey bar = 95% range when the same number of
diseases is removed at random (null columns of the table). Panel b: AUROC of the excluded pairs vs the kept pairs.
Usage: si_placement_confidence.py <dq_exclusion_sensitivity.tsv> <outdir>
"""
import sys, os, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]
import sys as _sys; _sys.path.insert(0, str(REPO / "src"))
import disatlas.figstyle_si as fs
SRC, OUT = str(REPO / "data" / "derived" / "placement_confidence" / "exclusion_sensitivity.tsv"), str(HERE); fs.apply_figure_style()
T = pd.read_csv(SRC, sep="\t"); assert len(T) == 20, len(T)
BEN = [("IxIDN-Matched 4,364 (all pools)", "IxIDN-Matched, all 4,364 pairs"),
       ("IxIDN-Matched test-test (held out)", "IxIDN-Matched, test split"),
       ("test_dm (held-out MechDist test split)", "Held-out link-prediction test split"),
       ("IxIDN-910 (tune split)", "Tuning split, 906 pairs")]
RULE = [("C5 unstable", "Unstable placement"), ("C6 ill-defined", "Ill-defined term"), ("C5 or C6", "Unstable or ill-defined"),
        ("Critical tier", "Critical tier"), ("Critical (evidence-only)", "Critical tier, evidence flags only")]
COL = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
crit = T[(T.benchmark == BEN[0][0]) & (T.exclusion == "Critical tier")].iloc[0]
assert crit.n_dis_excluded == 105 and crit.n_pairs_total - crit.n_pairs_kept == 1140
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.9), gridspec_kw=dict(width_ratios=[1.35, 1], wspace=.16))
y0 = np.arange(len(RULE))[::-1].astype(float); off = np.array([.27, .09, -.09, -.27])
for j, (bk, bl) in enumerate(BEN):
    S = T[T.benchmark == bk].set_index("exclusion").loc[[r for r, _ in RULE]]
    yy = y0 + off[j]
    for k in range(len(S)): a.plot([S.null_delta_lo.iloc[k], S.null_delta_hi.iloc[k]], [yy[k]] * 2, color="0.82", lw=3.2, solid_capstyle="butt", zorder=1)
    a.errorbar(S.delta, yy, xerr=[S.delta - S.delta_lo, S.delta_hi - S.delta], fmt="o", ms=3.5, color=COL[j], lw=.9, capsize=0, label=bl)
    b.scatter(S.auroc_kept, yy, s=14, color=COL[j], zorder=3)
    b.scatter(S.auroc_touching, yy, s=16, facecolor="white", edgecolor=COL[j], lw=.9, zorder=3)
    for k in range(len(S)): b.plot([S.auroc_kept.iloc[k], S.auroc_touching.iloc[k]], [yy[k]] * 2, color=COL[j], lw=.6, alpha=.6, zorder=2)
a.axvline(0, color="0.4", ls="--", lw=.7); a.set_xlim(-.055, .105); a.set_xticks(np.arange(-.04, .101, .02))
a.set_yticks(y0); a.set_yticklabels([l for _, l in RULE]); a.set_ylim(-.6, len(RULE) - .4)
a.set_xlabel("Change in AUROC after exclusion (95% CI)"); a.set_title("a  AUROC change when flagged diseases are removed")
b.set_yticks(y0); b.set_yticklabels([]); b.set_ylim(a.get_ylim()); b.set_xlabel("AUROC"); b.set_xlim(.74, .93); b.set_xticks(np.arange(.75, .926, .05))
b.set_title("b  Pairs touching flagged diseases (open) vs kept (filled)")
from matplotlib.lines import Line2D
h, l = a.get_legend_handles_labels(); h.append(Line2D([], [], color="0.82", lw=3.2)); l.append("Random removal of as many diseases (95% range)")
fig.legend(h, l, loc="lower center", ncol=2, bbox_to_anchor=(.5, -.17), handlelength=2.2, handletextpad=.6, columnspacing=1.6)
N = dict(png_bytes=fs.save(fig, OUT, "si_placement_confidence"), overlaps=fs.overlap_check(fig),
         critical_all=dict(n_excluded=int(crit.n_dis_excluded), pairs_excluded=int(crit.n_pairs_total - crit.n_pairs_kept),
                           auroc_full=crit.auroc_full, auroc_kept=crit.auroc_kept, delta=crit.delta, lo=crit.delta_lo, hi=crit.delta_hi,
                           auroc_touching=crit.auroc_touching, null_lo=crit.null_delta_lo, null_hi=crit.null_delta_hi, pct_null_ge_obs=crit.pct_null_ge_obs))
json.dump(N, open(f"{OUT}/si_placement_confidence_values.json", "w"), indent=1, default=float); print(N)
