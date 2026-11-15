"""Supplementary Figure SI_mechdist_embedding_pair_score: from the MechDist graph to the DisAtlas embedding and pair score.

  a real example from the gene-anchored MechDist training graph (benchmark pairs stripped): atypical hemolytic-uremic
    syndrome (aHUS, MONDO:0016244) and paroxysmal nocturnal hemoglobinuria (PNH, MONDO:0100244), an IxIDN+ pair with no
    direct edge, connected through shared neighbours (12 named; count and hub degree from s6_inputs.json)
  b embedding cosine of all 4,364 IxIDN-Matched pairs (2,182 per class) with the tune-split threshold
  c-d pair-score schematic and the embedding cosine of six known control pairs. These two panels are a static drawing
    (data/derived/si_mechdist_pair_score/schematic_cd_static.png); the six cosines it shows are listed in
    schematic_cd_control_cosines.tsv and written to the values JSON.
Panels a-b are drawn here and stacked over panels c-d at the same 2,481-pixel width.
Inputs: data/derived/si_mechdist_pair_score/{s6_inputs.json, ixidn_matched_cosine_all.tsv, schematic_cd_static.png,
schematic_cd_control_cosines.tsv}.

    python figures/SI/si_mechdist_embedding_pair_score.py
Writes SI_mechdist_embedding_pair_score.png and SI_mechdist_embedding_pair_score_values.json next to this script.
"""
import os, sys, io
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'src'))
import matplotlib; matplotlib.use('Agg')
from PIL import Image
from disatlas.figstyle_helpers import *
D = os.path.join(REPO, 'data', 'derived', 'si_mechdist_pair_score')
import json, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib.patches as mpatches
apply_figure_style()
CT = pd.read_csv(os.path.join(D, "schematic_cd_control_cosines.tsv"), sep="\t")
O = json.load(open(os.path.join(D, "s6_inputs.json")))
T = O["training_graph"]; G = dict(direct_AB=T["direct_AB"], hub_shared=T["hub_in_shared"], n_shared=T["n_shared"], hub_deg=T["hub_deg"])  # training graph (benchmark pairs stripped)
assert not G["direct_AB"] and G["hub_shared"]
NB = {s: None for s in O["figure_shared"]}; assert len(NB) == 12
GENETIC = ["MONDO:0009335", "MONDO:0013040", "MONDO:0013041", "MONDO:0013042", "MONDO:0013043", "MONDO:0035290"]
SHOWN = [("GENETIC", "Genetic aHUS subtypes (6)\nCFH susceptibility · MCP/CD46 · CFI\nCFB · C3 · complement-gene"),
         ("MONDO:0019737", "Thrombotic microangiopathy"), ("MONDO:0018896", "Thrombotic thrombocytopenic\npurpura"),
         ("MONDO:0008585", "HELLP syndrome"), ("MONDO:0005081", "Preeclampsia"), ("MONDO:0007407", "Cryoglobulinemic vasculitis")]
HUB = "MONDO:0005070"   # 'neoplasm' umbrella term: 12th shared neighbour, drawn as a footnote not a node
assert set(GENETIC) | {k for k, _ in SHOWN[1:]} | {HUB} == set(NB)

fig = plt.figure(figsize=(8.6, 4.0)); gs = fig.add_gridspec(1, 2, width_ratios=[1.7, 1], wspace=0.22)
axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
axA.set_xlim(-3.05, 4.05); axA.set_ylim(-2.2, 1.75); axA.axis("off")
focal_c, _ = focal_palette(["pair", "other"], focal="pair", focal_color="#c0392b"); edge_c = META_GREY
pA, pB = (-1.85, 0.0), (0.95, 0.0)
ys = np.linspace(1.25, -1.25, len(SHOWN)); xm = -0.45
for (k, lab), y in zip(SHOWN, ys):
    for p in (pA, pB): axA.plot([p[0], xm], [p[1], y], color=edge_c, lw=0.8, alpha=0.55, zorder=1)
    if k == "GENETIC":
        axA.add_patch(mpatches.FancyBboxPatch((xm - 0.09, y - 0.09), 0.18, 0.18, boxstyle="round,pad=0.02", fc="white", ec="#4a4a4a", lw=1.1, zorder=3))
    else:
        axA.scatter([xm], [y], s=70, color="white", edgecolor="#4a4a4a", lw=1.1, zorder=3)
    axA.text(xm, y + 0.13, lab, ha="center", va="bottom", fontsize=5.6, color="#222222", linespacing=1.1, zorder=4,
             bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.9))
# held-out direct edge
axA.plot([pA[0], pB[0]], [-1.75, -1.75], color=focal_c, lw=1.0, ls=(0, (3, 2)), zorder=1)
for p in (pA, pB): axA.plot([p[0], p[0]], [p[1], -1.75], color=focal_c, lw=1.0, ls=(0, (3, 2)), zorder=1)
axA.text((pA[0] + pB[0]) / 2, -1.82, "aHUS–PNH: IxIDN⁺ benchmark pair, no direct edge in the training graph",
         ha="center", va="top", fontsize=5.6, color=focal_c)
for p, lab in ((pA, "Atypical hemolytic-\nuremic syndrome\n(aHUS)"), (pB, "Paroxysmal nocturnal\nhemoglobinuria\n(PNH)")):
    axA.scatter([p[0]], [p[1]], s=150, color=focal_c, edgecolor=focal_c, zorder=5)
    L = p is pA
    axA.text(p[0] + (-0.17 if L else 0.1), p[1] + (0 if L else -0.2), lab, ha="right" if L else "left", va="center" if L else "top", fontsize=6.0, fontweight="bold", color=focal_c, linespacing=1.1, zorder=6,
             bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.9))
# Node2Vec -> embedding (no model-score label in panel a)
axA.annotate("", xy=(2.55, 0.08), xytext=(1.25, 0.08), arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.3))
axA.text(1.9, 0.16, "Node2Vec", ha="center", va="bottom", fontsize=6.3, fontweight="bold")
rng = np.random.default_rng(7); bw = 0.05
for row, (y0, lab) in enumerate(((0.6, "aHUS"), (-0.6, "PNH"))):
    v = rng.normal(0, 1, 14); v = v / np.abs(v).max() * 0.28
    for i, x in enumerate(v): axA.add_patch(mpatches.Rectangle((3.1 + i * (bw + 0.01), y0 + min(0, x)), bw, max(abs(x), 0.01), fc=focal_c, ec="none", alpha=0.85))
    axA.plot([3.05, 3.1 + 14 * (bw + 0.01)], [y0, y0], color="#999999", lw=0.5)
    axA.text(3.02, y0, lab, ha="right", va="center", fontsize=5.6, color="#222222")
axA.text(3.52, -1.2, "DisAtlas embeddings\n(32-d, 5-seed ensemble;\nbars illustrative)", ha="center", va="top", fontsize=5.6, color="#222222")
axA.set_title("Real example: aHUS and PNH in the MechDist training graph", loc="left", fontsize=8)
axA.text(-2.9, -2.05, f"{G['n_shared']} shared neighbours; 11 drawn in six groups. Not drawn: 'neoplasm' umbrella term (degree {G['hub_deg']:,})\n"
         "and the remaining shared neighbours.",
         ha="left", va="top", fontsize=5.4, color="#555555")

# Panel B (gene-anchored embeddings): embedding cosine on all 4,364 IxIDN-Matched pairs, from the Figure 4a score file
FE = os.path.join(D, "ixidn_matched_cosine_all.tsv")
SC = pd.read_csv(FE, sep="\t", usecols=["pool", "y", "c"]); RR = dict(threshold_youden_tune=O["threshold_youden_tune"], test_sensitivity=O["test_sensitivity"], test_specificity=O["test_specificity"]); thr = RR["threshold_youden_tune"]
te = SC[SC.pool == "test"]
assert len(SC) == 4364 and (SC.y == 1).sum() == 2182 and (SC.y == 0).sum() == 2182 and SC.c.notna().all()
assert (te.y == 1).sum() == 491 and (te.y == 0).sum() == 603 and round(thr, 2) == 0.20
assert abs(((te.c >= thr) & (te.y == 1)).sum() / 491 - RR["test_sensitivity"]) < 1e-9 and abs(((te.c < thr) & (te.y == 0)).sum() / 603 - RR["test_specificity"]) < 1e-9
dpos = SC.loc[SC.y == 1, "c"].values; dneg = SC.loc[SC.y == 0, "c"].values
mpos, mneg = float(np.median(dpos)), float(np.median(dneg))
json.dump(dict(n_pos=len(dpos), n_neg=len(dneg), median_pos=mpos, median_neg=mneg, threshold=thr, min=float(SC.c.min()), max=float(SC.c.max()),
               test_n_pos=491, test_n_neg=603, test_sensitivity=RR["test_sensitivity"], test_specificity=RR["test_specificity"], n_shared_train=G["n_shared"], hub_deg_train=G["hub_deg"],
               control_cosines={r.tag: round(float(r.cosine_score), 3) for r in CT.itertuples()}),
          open(os.path.join(HERE, "SI_mechdist_embedding_pair_score_values.json"), "w"), indent=1)
def _fm(v): return f"{v:.2f}".replace("-", "\u2212")
bins = np.linspace(-0.5, 1.0, 31)
pos_c, neg_c = focal_palette(["IxIDN+ (positive)", "IxIDN- (negative)"], focal="IxIDN+ (positive)", focal_color="#c0392b")
axB.hist(dneg, bins=bins, color=neg_c, alpha=0.6, label=f"IxIDN\u207b (negative), median {_fm(mneg)}", density=True)
axB.hist(dpos, bins=bins, color=pos_c, alpha=0.65, label=f"IxIDN\u207a (positive), median {_fm(mpos)}", density=True)
ytop = axB.get_ylim()[1] * 1.32; axB.set_ylim(0, ytop)
axB.vlines([mneg], 0, ytop * 0.66, color=neg_c, lw=1.1, ls=":"); axB.vlines([mpos], 0, ytop * 0.66, color=pos_c, lw=1.1, ls=":")
axB.vlines([thr], 0, ytop * 0.66, color="#333333", lw=0.9, ls="--")
axB.text(thr, ytop * 0.68, f"threshold {thr:.2f}", ha="center", va="bottom", fontsize=5.8, color="#333333")
axB.set_xlabel("Cosine similarity of the two diseases'\nDisAtlas embeddings (32-d, mean over 5 seeds)"); axB.set_ylabel("Density")
axB.set_title("IxIDN-Matched pairs separate by embedding cosine", loc="left", fontsize=8)
axB.legend(frameon=False, fontsize=6.3, loc="upper left", handlelength=1.2, handletextpad=0.5, borderaxespad=0.1)
axB.text(0.02, 0.83, f"n = {len(dpos):,} pairs each; dotted lines: medians;\ndashed line: tune-split threshold", transform=axB.transAxes, ha="left", va="top", fontsize=5.8, color="#555555")
axB.margins(x=0.02); set_frame(axB)
panel_letter(axA, "a"); panel_letter(axB, "b")
buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=300, bbox_inches="tight"); plt.close(fig); buf.seek(0)
top = Image.open(buf)
bg = Image.new("RGB", top.size, "white"); bg.paste(top, mask=top.split()[-1] if top.mode == "RGBA" else None); top = bg
bottom = Image.open(os.path.join(D, "schematic_cd_static.png")).convert("RGB"); W = bottom.size[0]; assert W == 2481
if top.size[0] != W:
    assert top.size[0] < W; pad = Image.new("RGB", (W, top.size[1]), "white"); pad.paste(top, (0, 0)); top = pad
out = Image.new("RGB", (W, top.size[1] + bottom.size[1]), "white"); out.paste(top, (0, 0)); out.paste(bottom, (0, top.size[1]))
out.save(os.path.join(HERE, "SI_mechdist_embedding_pair_score.png"), dpi=(300, 300))
print("SI_mechdist_embedding_pair_score", out.size)
