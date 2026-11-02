#!/usr/bin/env python
"""Supplementary figure: seed vs sampling uncertainty on IxIDN-Matched (current five-seed DisAtlas ensemble).

Layout matches the earlier FIGURES/uncertainty_quantification.png: (A) SD of AUROC across the five training seeds vs
across the disease-cluster bootstrap (B = 2,000), all pairs and test split, annotated with the ratio; (B) AUROC/AUPRC of
the ensemble with 95% disease-cluster bootstrap CIs and the no-skill line (0.5 for AUROC, positive fraction for AUPRC).
Inputs (current model only):
  SD   : disatlas/v3_revised_manuscript/si_figures/ixidn_si/v3.2.0/R3M5_seed_vs_sample_v3.2.0.tsv  (seed_sd, sample_sd_cluster)
  CI   : revision/v32_rewrite/headline_swap/out/core/headline_core_values.tsv (B = 2,000 disease-cluster bootstrap; n, n_pos)
         cross-checked against R3M5_ci_pair_vs_cluster_v3.2.0.tsv (cluster_ci_lo/hi), must agree to 1e-3.
usage: python figures/SI/si_uncertainty.py
"""
import sys as _s, pathlib as _p; _s.path.insert(0, str(_p.Path(__file__).resolve().parents[2] / "src")); import disatlas.deterministic  # noqa: E402,F401  reproducible PDF/SVG bytes
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]
import sys as _sys; _sys.path.insert(0, str(REPO / "src"))
DV = REPO / "data" / "derived" / "uncertainty"
SVS = str(DV / "seed_vs_sample_sd.tsv"); CIC = str(DV / "ci_pair_vs_cluster.tsv")
HC = str(DV / "headline_core_values.tsv")
OUT = str(HERE / "si_uncertainty.png")
STR = {"all": ("all_4364", "All pairs"), "test": ("test_test_heldout", "Test split")}

sv = pd.read_csv(SVS, sep="\t").query("metric == 'auroc'").set_index("stratum")
hc = pd.read_csv(HC, sep="\t").iloc[:, :8]
hc.columns = ["metric", "stratum", "n", "n_pos", "value", "lo", "hi", "B"]
hc = hc[hc.metric.isin(["AUROC p", "AUPRC p"]) & hc.stratum.isin(["all", "test"])].copy()
hc["metric"] = hc.metric.str.split().str[0]
assert len(hc) == 4 and (hc.B.astype(int) == 2000).all(), hc
cc = pd.read_csv(CIC, sep="\t").query("system == 'dis2vec_ensemble'").set_index(["stratum", "metric"])
for r in hc.itertuples():
    k = (STR[r.stratum][0], r.metric.lower()); c = cc.loc[k]
    assert abs(c.value - r.value) < 1e-3 and abs(c.cluster_ci_lo - r.lo) < 1e-3 and abs(c.cluster_ci_hi - r.hi) < 1e-3, (k, r)
N = {s: int(hc[hc.stratum == s].n.iloc[0]) for s in STR}; NP = {s: int(hc[hc.stratum == s].n_pos.iloc[0]) for s in STR}

rows = []
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
BLUE, ORANGE, GREY = "#3a6aa8", "#d9822b", "#7a7a7a"
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.8, 3.9), dpi=300, gridspec_kw=dict(width_ratios=[1, 1.2], wspace=0.32))

# ---- A
x = np.arange(2); w = 0.36
seed = [sv.loc[STR[s][0], "seed_sd"] for s in STR]; samp = [sv.loc[STR[s][0], "sample_sd_cluster"] for s in STR]
a1.bar(x - w / 2, seed, w, color=BLUE, label="Retraining with a new seed (5 seeds)")
a1.bar(x + w / 2, samp, w, color=ORANGE, label="Resampling the diseases evaluated")
for i, s in enumerate(STR):
    ratio = samp[i] / seed[i]
    a1.text(x[i] + w / 2, samp[i] + 0.0007, "%d\u00d7 larger" % round(ratio), ha="center", va="bottom", fontsize=8)
    rows += [("A", STR[s][1], "seed SD of AUROC", seed[i], SVS + ":seed_sd"),
             ("A", STR[s][1], "disease-cluster bootstrap SD of AUROC", samp[i], SVS + ":sample_sd_cluster"),
             ("A", STR[s][1], "ratio (label rounded)", ratio, "derived")]
a1.set_xticks(x); a1.set_xticklabels(["%s\n(n = %s)" % (STR[s][1], format(N[s], ",")) for s in STR], fontsize=8)
a1.set_ylabel("SD of AUROC"); a1.set_ylim(0, max(samp) * 1.45)
a1.legend(frameon=False, fontsize=7.5, loc="upper left")
a1.set_title("Test-set sampling, not the seed,\ndrives the uncertainty", loc="left", fontsize=9)
a1.text(-0.12, 1.06, "A", transform=a1.transAxes, fontsize=11, fontweight="bold")

# ---- B
ys, labs = [], []; y = 0
for s in STR:
    a2.text(0.41, y - 0.3, "%s (n = %s)" % (STR[s][1], format(N[s], ",")), fontsize=7.5, fontweight="bold", va="center")
    for m in ("AUROC", "AUPRC"):
        y -= 1; r = hc[(hc.stratum == s) & (hc.metric == m)].iloc[0]
        ns = 0.5 if m == "AUROC" else NP[s] / N[s]
        a2.plot([ns, ns], [y - 0.32, y + 0.32], color=GREY, lw=1.4)
        a2.plot([r.lo, r.hi], [y, y], color=BLUE, lw=2.4, solid_capstyle="butt"); a2.plot(r.value, y, "o", color=BLUE, ms=6)
        a2.text(r.hi + 0.012, y, "%.3f (%.3f\u2013%.3f)" % (r.value, r.lo, r.hi), va="center", fontsize=7.5)
        ys.append(y); labs.append(m)
        rows += [("B", STR[s][1], m, r.value, HC), ("B", STR[s][1], m + " 95% CI low", r.lo, HC),
                 ("B", STR[s][1], m + " 95% CI high", r.hi, HC), ("B", STR[s][1], m + " no-skill line", ns,
                 "0.5" if m == "AUROC" else HC + ":n_pos/n")]
    y -= 0.75
a2.text(0.505, ys[0] + 0.34, "no skill", fontsize=6.5, color=GREY, va="bottom")
a2.set_yticks(ys); a2.set_yticklabels(labs, fontsize=8); a2.tick_params(axis="y", length=0)
a2.set_xlim(0.4, 1.0); a2.set_ylim(min(ys) - 0.5, -0.1); a2.set_xlabel("Score")
a2.spines["left"].set_visible(True)
a2.set_title("DisAtlas on IxIDN-Matched\n(95% CI, disease-cluster bootstrap)", loc="left", fontsize=9)
a2.text(-0.14, 1.06, "B", transform=a2.transAxes, fontsize=11, fontweight="bold")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, bbox_inches="tight", facecolor="white"); fig.savefig(OUT[:-4] + ".pdf", bbox_inches="tight", facecolor="white")
rows = [r[:4] + (str(r[4]).replace(str(REPO) + "/", ""),) if isinstance(r, tuple) else r[:4] + [str(r[4]).replace(str(REPO) + "/", "")] for r in rows]
pd.DataFrame(rows, columns=["panel", "stratum", "quantity", "value", "source"]).to_csv(OUT[:-4] + "_values.tsv", sep="\t", index=False)
print("wrote", OUT)
