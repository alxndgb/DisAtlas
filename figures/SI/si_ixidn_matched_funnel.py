"""Supplementary Figure SI_ixidn_matched_funnel: how IxIDN-Matched is built.

  a  pairs left after each construction step, IxIDN- (left) and IxIDN+ (right); bar width proportional to sqrt(pairs)
  b  shortest-path hop distance in the trial co-investigation graph for all disease pairs and for the matched pairs

Inputs (repo-relative):
  data/derived/ixidn_funnel/funnel_counts.tsv      one row per stage and class
  data/derived/ixidn_funnel/hop_distribution.tsv   hop histogram for all pairs and each stage
Outputs next to this script: si_ixidn_matched_funnel.{png,pdf}, si_ixidn_matched_funnel_values.json
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from disatlas.figstyle_si import apply_figure_style

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DV = ROOT / "data" / "derived" / "ixidn_funnel"
apply_figure_style()

NEG, POS = "#c0392b", "#2e6db4"
NEG_L, POS_L = "#ebb3ad", "#a9c4e8"

fc = pd.read_csv(DV / "funnel_counts.tsv", sep="\t")
hd = pd.read_csv(DV / "hop_distribution.tsv", sep="\t")


def n(stage, cls):
    v = fc[(fc.stage.astype(str) == str(stage)) & (fc["class"] == cls)].pairs
    assert len(v) == 1, (stage, cls)
    return int(v.iloc[0])


n_dis = int(fc[fc.stage.astype(str) == "0"].what.str.extract(r"(\d+) diseases")[0].iloc[0])
m_dis = int(fc[(fc.stage.astype(str) == "6") & (fc["class"] == "IxIDN+")].what.str.extract(r"(\d+) diseases")[0].iloc[0])
neg = [(n(0, "all"), f"All disease pairs\n({n_dis:,} trial-graph diseases)"),
       (n(1, "IxIDN-"), "≥3 hops apart"),
       (n(2, "IxIDN-"), "Co-investigation veto"),
       (n(3, "IxIDN-"), "Both diseases in the\n985-disease universe"),
       (n(6, "IxIDN-"), "1:1 exact match to a\nQC-passed positive")]
pos = [(n(1, "IxIDN+"), "Co-investigated\n(shared trial drug, 1 hop)"),
       (n(3, "IxIDN+"), "Both diseases in the\n985-disease universe"),
       (n(4, "IxIDN+"), "Has a negative in its\n7-key shortcut cell"),
       (n(5, "IxIDN+"), "LLM mechanism QC:\nSTRONG or PLAUSIBLE"),
       (n(6, "IxIDN+"), "1:1 exact match")]

fig = plt.figure(figsize=(9.2, 8.6))
axA = fig.add_axes([0.0, 0.50, 1.0, 0.47])
axB = fig.add_axes([0.13, 0.06, 0.84, 0.36])

# ---- panel a: twin funnels, width ∝ sqrt(pairs) ----
axA.set_xlim(0, 10); axA.set_ylim(-1.9, 6.5); axA.axis("off")
W0 = 1.9 / np.sqrt(neg[0][0])
cxN, cxP, h = 3.6, 5.95, 0.85


def funnel(cx, rows, light, dark, side):
    for i, (v, lab) in enumerate(rows):
        y1 = 5.0 - i * 1.0; y0 = y1 - h
        wt = W0 * np.sqrt(v) / 2
        wb = W0 * np.sqrt(rows[i + 1][0]) / 2 if i + 1 < len(rows) else wt
        last = i == len(rows) - 1
        col = dark if last else light
        axA.add_patch(plt.Polygon([(cx - wt, y1), (cx + wt, y1), (cx + wb, y0), (cx - wb, y0)], closed=True, color=col, lw=0))
        tcol = dark if (last or side == "pos") else "black"
        if side == "neg" and not last:
            axA.text(cx, (y0 + y1) / 2, f"{v:,}", ha="center", va="center", fontsize=11, weight="bold", color="black")
        else:
            off = max(wt, 0.12) + 0.08
            axA.text(cx - off if side == "neg" else cx + off, (y0 + y1) / 2, f"{v:,}", ha="right" if side == "neg" else "left",
                     va="center", fontsize=11, weight="bold", color=tcol)
        if side == "neg":
            axA.text(2.75, (y0 + y1) / 2, lab, ha="right", va="center", fontsize=9.5, color="#333333")
        else:
            axA.text(6.85, (y0 + y1) / 2, lab, ha="left", va="center", fontsize=9.5, color="#333333")


funnel(cxN, neg, NEG_L, NEG, "neg")
funnel(cxP, pos, POS_L, POS, "pos")
axA.text(cxN, 5.45, "IxIDN⁻ (negatives)", ha="center", fontsize=11, weight="bold", color=NEG)
axA.text(cxP + 0.0, 5.45, "IxIDN⁺ (positives)", ha="center", fontsize=11, weight="bold", color=POS)
yb = 5.0 - 4 * 1.0 - h
for x0, x1 in [(cxN, 4.15), (cxP, 5.45)]:
    axA.annotate("", xy=(x1, yb - 0.45), xytext=(x0, yb - 0.08), arrowprops=dict(arrowstyle="->", color="#777777", lw=1))
tot = n(6, "IxIDN+") + n(6, "IxIDN-")
axA.text(4.8, yb - 0.95, f"IxIDN-Matched: {tot:,} pairs ({n(6, 'IxIDN+'):,} IxIDN⁺ + {n(6, 'IxIDN-'):,} IxIDN⁻), {m_dis:,} diseases\n"
         "out-of-distribution transfer test", ha="center", va="center", fontsize=10,
         bbox=dict(boxstyle="square,pad=0.5", fc="#f3f3f3", ec="#999999"))
axA.text(0.05, 6.45, "a", fontsize=13, weight="bold", va="top")
axA.text(1.6, 6.45, "How IxIDN-Matched is built: pairs left after each step (width ∝ √pairs)", fontsize=11, va="top")

# ---- panel b: hop distance histogram ----
order = ["disconnected", "1", "2", "3", "4", "5", "6"]
hd["hops"] = hd.hops.astype(str)
hd = hd.set_index("hops").loc[order]
x = np.arange(len(order))
light = ["#c5c3d2", POS_L, "#9b9b9b"] + [NEG_L] * 4
axB.bar(x, hd.all_pairs, width=0.8, color=light)
axB.bar(x, hd.matched_pos, width=0.45, color=POS)
axB.bar(x, hd.matched_neg, width=0.45, color=NEG)
for i, k in enumerate(order):
    a = int(hd.all_pairs[k]); axB.text(i, a * 1.15, f"{a:,}", ha="center", va="bottom", fontsize=8.5)
    for col, c in [("matched_pos", POS), ("matched_neg", NEG)]:
        v = int(hd[col][k])
        if v >= 10:
            axB.text(i, v * 0.7, f"{v:,}", ha="center", va="top", fontsize=8.5, color="white", weight="bold")
        elif v > 0:
            axB.text(i + 0.3, v * 1.2, f"{v:,}", ha="left", va="bottom", fontsize=8.5, color=c, weight="bold")
axB.tick_params(labelsize=10)
axB.set_yscale("log"); axB.set_ylim(1, 4e7)
axB.set_yticks([1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7]); axB.set_yticklabels(["1", "10", "100", "1k", "10k", "100k", "1M", "10M"])
axB.set_xticks(x); axB.set_xticklabels(["no path"] + order[1:])
axB.axvline(2.5, ls="--", color="#555555", lw=1)
axB.text(2.57, 2e7, "IxIDN⁻ sampled from ≥3 hops", color=NEG, fontsize=9, va="center")
axB.set_xlabel("Shortest path between the two diseases in the trial co-investigation graph (hops)", fontsize=11)
axB.set_ylabel("Disease pairs (log scale)", fontsize=11)
for s in ("top", "right"): axB.spines[s].set_visible(False)
from matplotlib.patches import Patch
axB.legend(handles=[Patch(color="#9b9b9b", label="All disease pairs (light shades)"), Patch(color=NEG, label="IxIDN-Matched IxIDN⁻"),
                    Patch(color=POS, label="IxIDN-Matched IxIDN⁺")], loc="upper right", frameon=False, fontsize=9)
fig.text(0.0, 0.455, "b", fontsize=13, weight="bold")
fig.text(0.13, 0.455, "Hop distance: matched positives are 1 hop, matched negatives ≥3 hops", fontsize=11)

fig.savefig(HERE / "si_ixidn_matched_funnel.png", dpi=300, facecolor="white", bbox_inches="tight")
fig.savefig(HERE / "si_ixidn_matched_funnel.pdf", facecolor="white", bbox_inches="tight")
vals = {"neg": [v for v, _ in neg], "pos": [v for v, _ in pos], "matched_pairs": tot, "matched_diseases": m_dis,
        "trial_graph_diseases": n_dis, "hops_all_pairs": {k: int(hd.all_pairs[k]) for k in order},
        "hops_matched_neg": {k: int(hd.matched_neg[k]) for k in order}}
json.dump(vals, open(HERE / "si_ixidn_matched_funnel_values.json", "w"), indent=1)
print(json.dumps(vals))
