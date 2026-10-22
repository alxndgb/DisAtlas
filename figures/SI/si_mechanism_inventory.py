"""Supplementary Figure SI_mechanism_sensitivity: mechanism-inventory sensitivity to GO-CAM and Reactome.

Appendix section: Mechanism inventory sensitivity (#sec:SI_mechanism_sensitivity_gocam,
#sec:SI_mechanism_sensitivity_reactome, #sec:SI_mechanism_sensitivity_landscape), Figure #fig:SI_mechanism_sensitivity.
  a  coverage of the four adopted evidence channels (GO, HPA, Pfam, HPO) and of the two pathway sources tested and
     not adopted (GO-CAM, Reactome): % of the 4,364 IxIDN-Matched pairs scored and % of the 13,021 network diseases.
  b  the two tested sources: Reactome as a fifth channel (coverage, change in AUROC with 95% disease-cluster CI over
     five retrained seeds, diseases whose edges change) and GO-CAM (coverage; dropped before training).
  c  overlap of each disease's 50 nearest neighbours: with vs without Reactome (same seed), two seeds of the same
     four-channel model, and the hypergeometric chance level.
The three panels are drawn separately and stacked (a on top, b and c side by side), as in the manuscript image.
Inputs: figures/SI/source_data/inventory_coverage.tsv, inventory_reactome_auroc.tsv, inventory_knn50_overlap.tsv,
        inventory_reactome_summary.json
Output: figures/SI/mechanism_inventory_sensitivity.{png,pdf}
Run:    python figures/SI/si_mechanism_inventory.py      (needs Pillow)
"""
import io, json, os, sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _style import src  # noqa: E402

CV = pd.read_csv(src("inventory_coverage.tsv"), sep="\t").set_index("key")
AU = pd.read_csv(src("inventory_reactome_auroc.tsv"), sep="\t")
KNN = pd.read_csv(src("inventory_knn50_overlap.tsv"), sep="\t")
META = json.load(open(src("inventory_reactome_summary.json")))
DPI = 250


def pc(k, num, den):
    return 100 * CV.loc[k, num] / CV.loc[k, den]


NNET = int(CV.loc["frozen", "n_net"])
RX = AU[AU.scope.str.startswith("IxIDN-Matched") & (AU.arm == "+ Reactome")].iloc[0]
chance = META["chance_hypergeom"]
same_seed = KNN.loc[KNN.kind == "with vs without Reactome, same seed", "overlap"].to_numpy()
floor = KNN.loc[KNN.kind == "two seeds of the four-channel model", "overlap"].to_numpy()
assert len(same_seed) == 5 and len(floor) == 10

# ---- value check: appendix #sec:SI_mechanism_sensitivity_* ----
assert [round(pc("gocam_causal", *c), 1) for c in (("bench_genes", "n_bench_genes"), ("pairs_scored", "n_pairs"),
        ("net_diseases", "n_net"))] == [28.2, 42.0, 28.8]
assert (round(pc("reactome_channel", "pairs_scored", "n_pairs"), 1), round(pc("reactome_channel", "net_diseases", "n_net"), 1)) == (89.0, 69.1)
assert META["n_diseases_touched"] == 556 and round(100 * 556 / NNET, 1) == 4.3
assert (round(RX.delta_vs_control, 3), round(RX.d_ci_lo, 3), round(RX.d_ci_hi, 3)) == (-0.003, -0.010, 0.004)
assert round(same_seed.mean(), 3) == 0.684 and (round(floor.min(), 3), round(floor.max(), 3)) == (0.675, 0.681)
assert round(chance, 3) == 0.004
print("appendix values: OK")

B, M, S = 9, 8, 7
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": B})
OLD = "#e6e6e6"; OLDE = "#666666"; NEW = "#fdf1c7"; NEWE = "#c98a00"; TEST = "#f2f2f2"; TESTE = "#8c8c8c"


def to_img(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=DPI, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    plt.close(fig); buf.seek(0)
    return Image.open(buf).convert("RGB")


def autocrop(im, pad=10):
    arr = np.asarray(im); nz = np.any(arr < 250, axis=2)
    r, c = np.where(nz.any(axis=1))[0], np.where(nz.any(axis=0))[0]
    return im.crop((max(c.min() - pad, 0), max(r.min() - pad, 0), min(c.max() + pad, im.size[0] - 1),
                    min(r.max() + pad, im.size[1] - 1)))


# ---- panel a: coverage ----
rows = [("go_only", "GO"), ("hpa_only", "HPA"), ("dom_only", "Pfam"), ("hpo_only", "HPO"),
        ("reactome_channel", "Reactome"), ("gocam_causal", "GO-CAM")]
pairs_pct = [pc(k, "pairs_scored", "n_pairs") for k, _ in rows]
dis_pct = [pc(k, "net_diseases", "n_net") for k, _ in rows]
face = [OLD, NEW, NEW, OLD, TEST, TEST]; edge = [OLDE, NEWE, NEWE, OLDE, TESTE, TESTE]; ls = ["-"] * 4 + ["--"] * 2
y = list(range(len(rows)))[::-1]; sep_y = y[3] - 0.5; h = 0.6
fig = plt.figure(figsize=(11.0, 3.0))
lab = fig.add_axes([0, 0, 1, 1]); lab.set_xlim(0, 110); lab.set_ylim(0, 1); lab.axis("off")
lab.text(0.009, 0.97, "a", fontsize=13, weight="bold", va="top", transform=lab.transAxes)
lab.text(0.035, 0.97, "Coverage of the four adopted evidence channels (GO, HPA, Pfam, HPO) and of the two\n"
         "pathway sources tested and not adopted (GO-CAM, Reactome)", fontsize=B + 1, va="top", transform=lab.transAxes)


def draw_panel(ax, vals, xlabel):
    for yi, fc, ec, l, v in zip(y, face, edge, ls, vals):
        ax.barh(yi, v, height=h, color=fc, edgecolor=ec, linestyle=l, linewidth=1.1)
        ax.text(v + 1.5, yi, f"{v:.1f}%", va="center", fontsize=S)
    ax.set_yticks(y); ax.set_ylim(min(y) - 0.7, max(y) + 0.7)
    ax.set_xlim(0, 130); ax.set_xticks([0, 20, 40, 60, 80, 100]); ax.set_xlabel(xlabel, fontsize=S)
    ax.spines[["top", "right"]].set_visible(False)
    ax.axhline(sep_y, color="#999999", lw=0.8, ls=":")
    ax.text(129, y[1], "adopted", ha="right", va="center", fontsize=S, color="#444444", style="italic")
    ax.text(129, y[5], "tested, not\nadopted", ha="right", va="center", fontsize=S, color="#444444", style="italic")


ax = fig.add_axes([0.14, 0.14, 0.345, 0.62]); draw_panel(ax, pairs_pct, "% of 4,364 IxIDN-Matched pairs scored")
ax.set_yticklabels([l for _, l in rows], fontsize=M)
ax2 = fig.add_axes([0.58, 0.14, 0.345, 0.62]); draw_panel(ax2, dis_pct, "% of 13,021 network diseases placed")
ax2.set_yticklabels([])
img_a = autocrop(to_img(fig))

# ---- panel b: the two tested pathway sources ----
fig = plt.figure(figsize=(11.0, 8.0))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 110); ax.set_ylim(0, 80); ax.axis("off")


def box(x, yb, w, hb, text, fc, ec, fs=M, ls="-", lw=1.1):
    ax.add_patch(FancyBboxPatch((x, yb), w, hb, boxstyle="round,pad=0.25,rounding_size=0.8", fc=fc, ec=ec, lw=lw, ls=ls))
    lines = text.split("\n")
    ax.text(x + w / 2, yb + hb - 1.15, lines[0], ha="center", va="top", fontsize=fs, weight="bold")
    ax.text(x + w / 2, yb + hb - 3.2, "\n".join(lines[1:]), ha="center", va="top", fontsize=S, linespacing=1.25)


cvt = (f"{pc('reactome_channel', 'pairs_scored', 'n_pairs'):.1f}% of pairs,\n"
       f"{pc('reactome_channel', 'net_diseases', 'n_net'):.1f}% of diseases")
rxn_text = ("Reactome pathway membership\nadded as a fifth channel,\nshared-gene disease pairs only\n"
            f"channel defined for {cvt}\n\n"
            + f"ΔAUROC {RX.delta_vs_control:+.3f}\n(95% CI {RX.d_ci_lo:+.3f} to {RX.d_ci_hi:+.3f}; 5 seeds,\nIxIDN-Matched)\n\n".replace("-0", "−0")
            + f"changes the edges of {META['n_diseases_touched']:,}\nof {NNET:,} diseases ({100 * META['n_diseases_touched'] / NNET:.1f}%)")
ax.text(1, 27.0, "b", fontsize=13, weight="bold", va="top")
ax.text(4, 27.0, f"{pc('reactome_channel', 'net_diseases', 'n_net'):.0f}% and {pc('gocam_causal', 'net_diseases', 'n_net'):.0f}% "
        "disease coverage: pathway channels\n(Reactome, GO-CAM) tested, not adopted", fontsize=B + 1, va="top")
box(1.5, 3.0, 25.5, 20.5, rxn_text, TEST, TESTE, ls="--")
box(29.0, 3.0, 23.5, 20.5, "GO-CAM causal activity graph\n\n"
    f"covers {pc('gocam_causal', 'bench_genes', 'n_bench_genes'):.1f}% of benchmark genes,\n"
    f"{pc('gocam_causal', 'pairs_scored', 'n_pairs'):.1f}% of benchmark pairs\nand "
    f"{pc('gocam_causal', 'net_diseases', 'n_net'):.1f}% of {NNET:,} diseases\n\ndropped on coverage,\nbefore training a model",
    TEST, TESTE, ls="--")
img_b = autocrop(to_img(fig))

# ---- panel c: 50-nearest-neighbour overlap ----
groups = [("chance\n(hypergeometric null)", np.array([chance]), "#999999", "#555555"),
          ("two seeds of the same\nfour-channel model\n(retraining floor)", floor, "#9ecae1", "#2b6ba3"),
          ("with vs without Reactome\n(same seed)", same_seed, NEWE, "#8a5a00")]
fig = plt.figure(figsize=(5.4, 3.0))
lab = fig.add_axes([0, 0, 1, 1]); lab.set_xlim(0, 1); lab.set_ylim(0, 1); lab.axis("off")
lab.text(0.015, 0.97, "c", fontsize=13, weight="bold", va="top", transform=lab.transAxes)
lab.text(0.09, 0.97, "Mechanism-inventory choice barely moves the\nembedding's neighbourhood structure",
         fontsize=B + 1, va="top", transform=lab.transAxes)
axc = fig.add_axes([0.40, 0.14, 0.56, 0.60]); yc = [2, 1, 0]; rng = np.random.default_rng(0)
for yi, (gl, vals, fc, ec) in zip(yc, groups):
    jit = rng.uniform(-0.12, 0.12, size=len(vals)) if len(vals) > 1 else np.array([0.0])
    axc.scatter(vals, np.full(len(vals), yi) + jit, s=34, facecolor=fc, edgecolor=ec, linewidth=0.9, zorder=3)
    if len(vals) > 1:
        axc.plot([vals.min(), vals.max()], [yi, yi], color=ec, lw=1.4, zorder=2); txt = f"{vals.min():.3f}-{vals.max():.3f}"
    else:
        txt = f"{vals[0]:.4f}"
    axc.text(1.0, yi, txt, va="center", fontsize=S, transform=axc.get_yaxis_transform())
axc.set_yticks(yc); axc.set_yticklabels([g[0] for g in groups], fontsize=M)
axc.set_xlim(0, 0.78); axc.set_xticks([0, 0.2, 0.4, 0.6])
axc.set_xlabel(f"fraction of 50 nearest neighbours shared\n({META['n_common_nodes']:,} common diseases, 5 seeds)", fontsize=S)
axc.set_ylim(-0.6, 2.6); axc.spines[["top", "right"]].set_visible(False)
img_c = autocrop(to_img(fig))

# ---- compose: a on top; b and c side by side ----
th = max(img_b.size[1], img_c.size[1])
b_r = img_b.resize((round(img_b.size[0] * th / img_b.size[1]), th))
c_r = img_c.resize((round(img_c.size[0] * th / img_c.size[1]), th))
gap = 40; rw = b_r.size[0] + gap + c_r.size[0]
row = Image.new("RGB", (rw, th), (255, 255, 255)); row.paste(b_r, (0, 0)); row.paste(c_r, (b_r.size[0] + gap, 0))
tw = max(img_a.size[0], rw)
a_r = img_a.resize((tw, round(img_a.size[1] * tw / img_a.size[0])))
row_r = row.resize((tw, round(th * tw / rw)))
out = Image.new("RGB", (tw, a_r.size[1] + row_r.size[1] + 20), (255, 255, 255))
out.paste(a_r, (0, 0)); out.paste(row_r, (0, a_r.size[1] + 20))
for ext in ("png", "pdf"):
    p = os.path.join(HERE, f"mechanism_inventory_sensitivity.{ext}")
    out.save(p, dpi=(DPI, DPI)) if ext == "png" else out.save(p, resolution=DPI)
    print(p)
