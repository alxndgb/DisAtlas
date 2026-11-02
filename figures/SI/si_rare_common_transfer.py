#!/usr/bin/env python
"""SI figure: DisAtlas performance by disease rarity (test-split strata) and by anchor disease, current model.

Re-render of FIGURES/SI_rare_common_transfer.png (
frozen previous-release inputs). Layout, panels, markers, colours and annotations are kept; only the inputs change.
  a  test-split IxIDN-Matched pairs by rarity stratum: DisAtlas AUROC with 95% disease-cluster bootstrap CI
     (B = 2000; ixidn_ccplus_strata.tsv, benchmark IxIDN-Matched, scope 'Test split pairs'); HPO Resnik-BMA and
     GO Resnik gene-BMA AUROC on the same pairs (r4mc3_stratified_auroc_v3.2.0.tsv). Strata with < 20 pairs in a
     class get no AUROC ("n too small").
  b  anchored analysis, all pairs: mean per-anchor AUROC over qualifying anchors (>= 1 positive and >= 1 negative
     partner of the stated class), from the per-anchor tables r4mc3_*_per_anchor_v3.2.0.tsv; 95% CI = anchor
     bootstrap recomputed here with B = 2000 (seed 4707; the on-disk SI_anchored_four_groups table used B = 1000).
     Groups with < 50 negative pairs are drawn hollow.
  c  (omitted with --two-panel) all IxIDN-Matched pairs by stratum plus the supplementary common-common benchmark
     with added matched negatives: DisAtlas AUROC (95% CI, B = 2000) against the combined matched-shortcut stack;
     the difference and its 95% CI are printed above each stratum (ccplus_strata.tsv). The IxIDN-Matched all-pairs
     intervals (DisAtlas and the paired difference) are on the Figure 5c disease-cluster draws (938 clusters);
     test-split and common-common benchmark rows keep their per-stratum disease-cluster intervals.
Outputs: next to this script: SI_rare_common_transfer.{png,pdf} (the three-panel Supplementary Figure) and
    SI_rare_common_transfer_values.tsv; --two-panel writes panels a-b only as SI_rare_common_transfer_2panel.*
    python si_rare_common_transfer.py [--two-panel]
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from disatlas.disatlas_style import apply_style, panel_label, ANNOT_SIZE  # house style of the old figure

L = DERIVED / "si_rare_common"
SRC = dict(ccplus=f"{L}/ccplus_strata.tsv", r4mc3=f"{L}/r4mc3.tsv", four=f"{L}/four.tsv", per_anchor=f"{L}/{{g}}.tsv")
OUTDIR = str(HERE)
WITH_C = "--two-panel" not in sys.argv
NAME = "SI_rare_common_transfer" if WITH_C else "SI_rare_common_transfer_2panel"

C = pd.read_csv(SRC["ccplus"], sep="\t")
T = pd.read_csv(SRC["r4mc3"], sep="\t")
A4 = pd.read_csv(SRC["four"], sep="\t").set_index("group")

DA_COL = "#d9822b"; SC_COL = "0.45"
BL = [("HPO Resnik-BMA", "HPO similarity", "#7fa6cf", "s"), ("GO Resnik gene-BMA", "GO similarity", "#4c9a6a", "^")]
ST = [("rare-rare", "Rare–\nrare"), ("rare-common", "Rare–\ncommon"), ("common-common", "Common–\ncommon"), ("all", "All\npairs")]
DIR = [("rare_anchor_rank_common_partners", "Rare\nanchor,\ncommon\npartners"),
       ("common_anchor_rank_rare_partners", "Common\nanchor,\nrare\npartners"),
       ("rare_anchor_rank_rare_partners", "Rare\nanchor,\nrare\npartners"),
       ("common_anchor_rank_common_partners", "Common\nanchor,\ncommon\npartners")]
MIN_N, MIN_ANCH, SPARSE_NEG, B_ANCH = 20, 5, 50, 2000


def cc_row(bench, scope, st):
    r = C[(C.benchmark == bench) & (C.scope == scope) & (C.stratum == st)]
    assert len(r) == 1, (bench, scope, st); return r.iloc[0]


def t_row(scope, st):
    r = T[(T.scope == scope) & (T.stratum == st)]
    assert len(r) == 1, (scope, st); return r.iloc[0]


# ---- consistency checks: point estimates and counts agree between the two stratified tables ----
for sc_c, sc_t in [("Test split pairs", "Held-out test-test pairs"), ("All pairs", "All 4364 pairs")]:
    for st, _ in ST:
        a, b = cc_row("IxIDN-Matched", sc_c, st), t_row(sc_t, st)
        assert (a.n_pos, a.n_neg) == (b.n_pos, b.n_neg), (sc_c, st)
        assert (pd.isna(a.auroc) and pd.isna(b.auroc)) or abs(a.auroc - b.auroc) < 1e-6, (sc_c, st, a.auroc, b.auroc)


# ---- anchored groups: mean over anchors, anchor bootstrap B = 2000 ----
def anchored(g):
    d = pd.read_csv(SRC["per_anchor"].format(g=g), sep="\t")
    a4 = A4.loc[g]
    assert len(d) == a4.n_anchors and abs(d.auroc.mean() - a4.mean_auroc) < 6e-5, (g, len(d), d.auroc.mean())
    rng = np.random.default_rng(4707); v = d.auroc.values
    bs = v[rng.integers(0, len(v), (B_ANCH, len(v)))].mean(1)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return dict(mean=v.mean(), lo=lo, hi=hi, hpo=d.auroc_hpo.mean(), go=d.auroc_go.mean(),
                n_anchors=len(d), n_pairs=int(a4.n_pairs), n_neg=int(a4.n_neg), lo_b1000=a4.ci_lo, hi_b1000=a4.ci_hi)


apply_style()
ncol = 3 if WITH_C else 2
fig, ax = plt.subplots(1, ncol, figsize=(7.0 * ncol, 6.4), gridspec_kw=dict(wspace=0.22))
MS, CAP, LW, S = 10, 4, 2.0, 60
PLOTTED = []


def da_point(a, x, y, lo, hi, hollow=False, col=DA_COL, fmt="o"):
    a.errorbar(x, y, yerr=[[y - lo], [hi - y]], fmt=fmt, color=col, ms=MS, capsize=CAP, lw=LW, zorder=3,
               mfc="white" if hollow else col, mec=col, mew=LW)


def bl_point(a, x, y, col, m, hollow=False):
    a.scatter(x, y, marker=m, s=S, zorder=2, facecolors="white" if hollow else col, edgecolors=col,
              linewidths=1.6 if hollow else 0)


# panel a
a = ax[0]
for i, (st, lab) in enumerate(ST):
    r, t = cc_row("IxIDN-Matched", "Test split pairs", st), t_row("Held-out test-test pairs", st)
    nm = lab.replace("\n", "")
    if min(r.n_pos, r.n_neg) < MIN_N:
        a.text(i + 0.1, 0.62, "n too small", ha="center", va="center", fontsize=ANNOT_SIZE, color="0.35", style="italic")
        PLOTTED.append(("a", nm, "", "", "", "", "", int(r.n_pos), int(r.n_neg), "test split; < 20 negatives, not plotted"))
    else:
        da_point(a, i, r.auroc, r.ci_lo, r.ci_hi)
        for k, (c_, _, col, m) in enumerate(BL):
            bl_point(a, i + 0.14 + 0.08 * k, t[f"{c_} | auroc"], col, m)
        PLOTTED.append(("a", nm, r.auroc, r.ci_lo, r.ci_hi, t["HPO Resnik-BMA | auroc"], t["GO Resnik gene-BMA | auroc"],
                        int(r.n_pos), int(r.n_neg), "DisAtlas+CI: ccplus strata (IxIDN-Matched, Test split pairs, B=2000); HPO/GO: r4mc3 table"))
    a.text(i + 0.1, 0.355, f"{int(r.n_pos)}/{int(r.n_neg)}", ha="center", va="bottom", fontsize=ANNOT_SIZE, color="0.35")
a.set_xticks(range(len(ST))); a.set_xticklabels([l for _, l in ST]); a.set_xlim(-0.45, len(ST) - 0.35)
a.set_title("Test-split pairs", loc="left"); a.set_ylabel("AUROC (DisAtlas: 95% CI)")

# panel b
b = ax[1]; HOLLOW_B = []
for i, (g, lab) in enumerate(DIR):
    d = anchored(g); nm = lab.replace("\n", " ")
    if d["n_anchors"] < MIN_ANCH:
        b.text(i, 0.62, "too few\nnegatives", ha="center", va="center", fontsize=ANNOT_SIZE, color="0.35", style="italic")
        PLOTTED.append(("b", nm, "", "", "", "", "", d["n_anchors"], d["n_pairs"], "not plotted"))
    else:
        hol = d["n_neg"] < SPARSE_NEG
        if hol: HOLLOW_B.append(d["n_neg"])
        da_point(b, i, d["mean"], d["lo"], d["hi"], hollow=hol)
        bl_point(b, i + 0.14, d["hpo"], BL[0][2], "s", hollow=hol); bl_point(b, i + 0.22, d["go"], BL[1][2], "^", hollow=hol)
        PLOTTED.append(("b", nm, d["mean"], d["lo"], d["hi"], d["hpo"], d["go"], d["n_anchors"], d["n_pairs"],
                        f"per-anchor table, anchor bootstrap B=2000 seed 4707 (on-disk B=1000 CI {d['lo_b1000']:.3f}-{d['hi_b1000']:.3f})"))
    b.text(i + 0.07, 0.355, f"{d['n_anchors']} anchors\n{d['n_pairs']} pairs", ha="center", va="bottom",
           fontsize=ANNOT_SIZE, color="0.35")
b.set_xticks(range(len(DIR))); b.set_xticklabels([l for _, l in DIR]); b.set_xlim(-0.45, len(DIR) - 0.35)
b.set_ylabel("Mean per-anchor AUROC"); b.set_title("Anchored, all pairs", loc="left")

# panel c (optional)
if WITH_C:
    cx = ax[2]
    SC = [("IxIDN-Matched", "rare-rare", "Rare–\nrare"), ("IxIDN-Matched", "rare-common", "Rare–\ncommon"),
          ("IxIDN-Matched", "common-common", "Common–\ncommon"),
          ("IxIDN-Matched-CC+", "common-common", "Common–\ncommon,\nadded\nnegatives")]
    for i, (bench, st, lab) in enumerate(SC):
        r = cc_row(bench, "All pairs", st)
        sig = r.delta_vs_stack_lo > 0
        hol = r.n_neg < SPARSE_NEG
        da_point(cx, i, r.auroc, r.ci_lo, r.ci_hi, hollow=hol)
        cx.scatter(i + 0.18, r.shortcut_stack_matched, marker="D", s=S * 0.8, zorder=2, color=SC_COL,
                   facecolors="white" if hol else SC_COL, edgecolors=SC_COL, linewidths=1.6 if hol else 0)
        sgn = lambda v: f"{v:+.3f}".replace("-", "−")
        cx.text(i + 0.07, 0.383, f"Δ {sgn(r.delta_vs_stack)}\n({sgn(r.delta_vs_stack_lo)} to\n{sgn(r.delta_vs_stack_hi)})",
                ha="center", va="bottom", fontsize=ANNOT_SIZE, color="0.15" if sig else "0.45",
                fontweight="bold" if sig else "normal")
        cx.text(i + 0.07, 0.355, f"{int(r.n_pos)}/{int(r.n_neg)}", ha="center", va="bottom", fontsize=ANNOT_SIZE, color="0.35")
        PLOTTED.append(("c", lab.replace("\n", " ").replace("– ", "–"), r.auroc, r.ci_lo, r.ci_hi, r.shortcut_stack_matched,
                        f"delta {r.delta_vs_stack:.3f} ({r.delta_vs_stack_lo:.3f} to {r.delta_vs_stack_hi:.3f})",
                        int(r.n_pos), int(r.n_neg), f"ccplus strata ({bench}, All pairs, B=2000)"))
    cx.set_xticks(range(len(SC))); cx.set_xticklabels([l for _, _, l in SC]); cx.set_xlim(-0.45, len(SC) - 0.35)
    cx.set_ylabel("AUROC (DisAtlas: 95% CI)"); cx.set_title("All pairs, against matched shortcuts", loc="left")

for a_, letter in zip(ax, "abc"):
    a_.axhline(0.5, color="0.8", ls="--", lw=1.0, zorder=1)
    a_.set_ylim(0.33, 1.0)
    a_.set_yticks(np.arange(0.4, 1.01, 0.1))
    panel_label(a_, letter)
if WITH_C:
    ax[2].text(-0.42, 0.475, "Δ DisAtlas − stack (95% CI); bold = CI excludes 0", fontsize=ANNOT_SIZE, color="0.35", va="bottom")

h = [Line2D([], [], color=DA_COL, marker="o", ms=MS, ls="", label="DisAtlas")] + \
    [Line2D([], [], color=col, marker=m, ms=8, ls="", label=lab) for _, lab, col, m in BL]
if WITH_C:
    h.append(Line2D([], [], color=SC_COL, marker="D", ms=7, ls="", label="matched-shortcut stack"))
assert len(set(HOLLOW_B)) <= 1
if HOLLOW_B:
    h.append(Line2D([], [], color="0.35", marker="o", mfc="white", mec="0.35", mew=1.6, ms=MS, ls="",
                    label=f"hollow = few negatives (n = {HOLLOW_B[0]})" if not WITH_C else "hollow = fewer than 50 negatives"))
fig.legend(handles=h, frameon=False, ncol=len(h), loc="upper center", bbox_to_anchor=(0.5, -0.04))

os.makedirs(OUTDIR, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(f"{OUTDIR}/{NAME}.{ext}", dpi=300, bbox_inches="tight")
cols = ["panel", "stratum", "disatlas_auroc", "ci_lo", "ci_hi", "comparator1", "comparator2", "n_pos_or_anchors",
        "n_neg_or_pairs", "source"]
V = pd.DataFrame(PLOTTED, columns=cols)
for c_ in ["disatlas_auroc", "ci_lo", "ci_hi", "comparator1"]:
    V[c_] = pd.to_numeric(V[c_], errors="coerce").round(3)
V["comparator1"] = V["comparator1"].astype(object)
V.loc[V.panel != "c", "comparator2"] = pd.to_numeric(V.loc[V.panel != "c", "comparator2"], errors="coerce").round(3)
V.to_csv(f"{OUTDIR}/{NAME}_values.tsv", sep="\t", index=False)
print(f"{OUTDIR}/{NAME}.png")
