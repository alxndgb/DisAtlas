"""Appendix figure #fig:SI_mechdist_embedding_schematic (section: From the MechDist network to DisAtlas
embeddings).

(a) Real example from the MechDist training graph: atypical hemolytic-uremic syndrome (aHUS, MONDO:0016244)
    and paroxysmal nocturnal hemoglobinuria (PNH, MONDO:0100244). Their direct MechDist edge is an IxIDN+
    benchmark pair and is removed before training; in the training graph they stay two hops apart through
    12 shared neighbours.
(b) Euclidean distance between the two diseases' DisAtlas embeddings on the IxIDN-Matched benchmark.

Inputs (all shipped in git):
    data/03_network/mechdist_network_edges.parquet.gzip        full MechDist network
    figures/SI/source_data/mechdist_training_strip_edges.tsv.gz benchmark pairs (IxIDN+ and IxIDN-) present
                                                                 in the network and removed before training
    figures/SI/source_data/mechdist_schematic_pair_distance.tsv per-pair embedding distance (mean over the 5
                                                                 seeds of the per-seed Euclidean distance)
    figures/SI/source_data/mechdist_schematic_example.json      model score of the aHUS-PNH pair and its
                                                                 position among random disease pairs
The training graph is rebuilt and checked against the training edge-set SHA-256 recorded for the frozen model.

    python figures/SI/si_mechdist_schematic.py               # figure from the shipped tables
    python figures/SI/si_mechdist_schematic.py --from-model  # also recompute distances and the example
                                                             # score from the frozen model package in
                                                             # data/zenodo/models/ and check them
"""
import argparse
import gzip
import hashlib
import json
import os
import sys

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import REPO_ROOT, apply_compact_style, compact_panel_letter as panel_letter, muted_focal_pair, src, save

ap = argparse.ArgumentParser()
ap.add_argument("--from-model", action="store_true")
args = ap.parse_args()

A, B = "MONDO:0016244", "MONDO:0100244"
# edge-set hash over the shipped (blinded) disease IDs; same edge set as the frozen model's training graph
# (FROZEN_MECHDIST_MODEL.json: expected_train_edgeset_sha256_blinded_ids)
TRAIN_EDGES, TRAIN_SHA = 325837, "b7ae48e4bfae735956838b76982f0577ff3f29efcbb75ec8833da312d5234a83"

# ---------------------------------------------------------------- training graph (panel a)
E = pd.read_parquet(os.path.join(REPO_ROOT, "data/03_network/mechdist_network_edges.parquet.gzip"))
strip = pd.read_csv(src("mechdist_training_strip_edges.tsv.gz"), sep="\t")
G = nx.from_pandas_edgelist(E, "disease_a", "disease_b", ["channel_class"])
assert G.has_edge(A, B), "direct aHUS-PNH edge missing from the full network"
G.remove_edges_from(zip(strip.disease_a, strip.disease_b))
h = hashlib.sha256("\n".join(sorted("\t".join(sorted((u, v))) for u, v in G.edges())).encode()).hexdigest()
assert G.number_of_edges() == TRAIN_EDGES and h == TRAIN_SHA, "training graph differs from the frozen model's"
assert not G.has_edge(A, B)
shared = sorted(nx.common_neighbors(G, A, B))
NB = {n: {"train_degree": G.degree(n), "channel_to_a": G[A][n]["channel_class"],
          "channel_to_b": G[B][n]["channel_class"]} for n in shared}
assert len(NB) == 12 and nx.shortest_path_length(G, A, B) == 2           # appendix: 12 shared neighbours, 2 hops
assert (G.degree(A), G.degree(B)) == (37, 27)

GENETIC = ["MONDO:0009335", "MONDO:0013040", "MONDO:0013041", "MONDO:0013042", "MONDO:0013043", "MONDO:0035290"]
SHOWN = [("GENETIC", "Genetic aHUS subtypes (6)\nCFH susceptibility · MCP/CD46 · CFI\nCFB · C3 · complement-gene"),
         ("MONDO:0019737", "Thrombotic microangiopathy"), ("MONDO:0018896", "Thrombotic thrombocytopenic\npurpura"),
         ("MONDO:0008585", "HELLP syndrome"), ("MONDO:0005081", "Preeclampsia"), ("MONDO:0007407", "Cryoglobulinemic vasculitis")]
HUB = "MONDO:0005070"   # 'neoplasm' umbrella term: 12th shared neighbour, drawn as a footnote, not a node
assert set(GENETIC) | {k for k, _ in SHOWN[1:]} | {HUB} == set(NB)

# ---------------------------------------------------------------- distances and example score
pdist = pd.read_csv(src("mechdist_schematic_pair_distance.tsv"), sep="\t")
EX = json.load(open(src("mechdist_schematic_example.json")))

if args.from_model:
    pkg = os.path.join(REPO_ROOT, "data/zenodo/models")
    sys.path.insert(0, pkg)
    from frozen_model_loader import FrozenMechDistModel
    m = FrozenMechDistModel(pkg)
    pairs = pd.read_csv(os.path.join(REPO_ROOT, "data/08S_eval/ixidn4364_scores.tsv"), sep="\t")
    pairs = pairs.set_index("blind_id").loc[pdist.blind_id]
    X = []
    for s in (4707, 11, 202, 3407, 8191):                                  # the 5 ensemble seeds
        z = np.load(os.path.join(pkg, "embeddings", f"node2vec_s{s}.npz"))
        ix = {d: i for i, d in enumerate(z["ids"])}
        X.append(np.linalg.norm(z["X"][pairs.a.map(ix).values] - z["X"][pairs.b.map(ix).values], axis=1))
    d = np.mean(X, axis=0)
    assert np.allclose(d, pdist.euclidean_distance.values, atol=1e-6) and (pairs.y.values == pdist.y.values).all()
    ids = list(m.ids); n = len(ids); pos = {d: i for i, d in enumerate(ids)}
    p = float(m.score(A, B))
    rowA, rowB = m.score_pairs([A] * n, ids), m.score_pairs([B] * n, ids)
    i = np.random.default_rng(0).integers(0, n, 200000); j = np.random.default_rng(1).integers(0, n, 200000)
    frac = float((m.score_pairs([ids[x] for x in i], [ids[x] for x in j]) < p).mean())
    got = {"p": p, "rank_of_b_among_a_partners": int((rowA > rowA[pos[B]]).sum()) + 1,
           "rank_of_a_among_b_partners": int((rowB > rowB[pos[A]]).sum()) + 1, "frac_random_pairs_below_p": frac}
    for k, v in got.items():
        assert np.isclose(v, EX[k], atol=1e-6), (k, v, EX[k])
    print("model recompute matches the shipped tables")

dpos = pdist.loc[pdist.y == 1, "euclidean_distance"].values
dneg = pdist.loc[pdist.y == 0, "euclidean_distance"].values
assert len(dpos) == len(dneg) == 2182                                        # appendix: n = 2,182 pairs each
assert (round(np.median(dpos), 2), round(np.median(dneg), 2)) == (1.14, 1.40)  # appendix: median 1.14 vs 1.40
assert f"{EX['p']:.2f}" == "0.89" and f"{100 * (1 - EX['frac_random_pairs_below_p']):.1f}" == "0.6"
print("OK: appendix values reproduced (12 shared neighbours, p = 0.89, top 0.6%, medians 1.14 / 1.40)")


# ---------------------------------------------------------------- figure style
apply_compact_style()
META_GREY = "#888888"

# ---------------------------------------------------------------- panel a
fig = plt.figure(figsize=(8.6, 4.0)); gs = fig.add_gridspec(1, 2, width_ratios=[1.7, 1], wspace=0.22)
axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
axA.set_xlim(-3.05, 4.05); axA.set_ylim(-2.2, 1.75); axA.axis("off")
focal_c, _ = muted_focal_pair("#c0392b"); edge_c = META_GREY
pA, pB = (-1.85, 0.0), (0.95, 0.0)
ys = np.linspace(1.25, -1.25, len(SHOWN)); xm = -0.45
for (k, lab), y in zip(SHOWN, ys):
    for p_ in (pA, pB): axA.plot([p_[0], xm], [p_[1], y], color=edge_c, lw=0.8, alpha=0.55, zorder=1)
    if k == "GENETIC":
        axA.add_patch(mpatches.FancyBboxPatch((xm - 0.09, y - 0.09), 0.18, 0.18, boxstyle="round,pad=0.02",
                                              fc="white", ec="#4a4a4a", lw=1.1, zorder=3))
    else:
        axA.scatter([xm], [y], s=70, color="white", edgecolor="#4a4a4a", lw=1.1, zorder=3)
    axA.text(xm, y + 0.13, lab, ha="center", va="bottom", fontsize=5.6, color="#222222", linespacing=1.1, zorder=4,
             bbox=dict(fc="white", ec="none", pad=0.4, alpha=0.9))
# held-out direct edge
axA.plot([pA[0], pB[0]], [-1.75, -1.75], color=focal_c, lw=1.0, ls=(0, (3, 2)), zorder=1)
for p_ in (pA, pB): axA.plot([p_[0], p_[0]], [p_[1], -1.75], color=focal_c, lw=1.0, ls=(0, (3, 2)), zorder=1)
axA.text((pA[0] + pB[0]) / 2, -1.82, "direct aHUS–PNH edge: IxIDN⁺ benchmark pair, removed before training",
         ha="center", va="top", fontsize=5.6, color=focal_c)
for p_, lab in ((pA, "Atypical hemolytic-\nuremic syndrome\n(aHUS)"), (pB, "Paroxysmal nocturnal\nhemoglobinuria\n(PNH)")):
    axA.scatter([p_[0]], [p_[1]], s=150, color=focal_c, edgecolor=focal_c, zorder=5)
    L = p_ is pA
    axA.text(p_[0] + (-0.17 if L else 0.1), p_[1] + (0 if L else -0.2), lab, ha="right" if L else "left",
             va="center" if L else "top", fontsize=6.0, fontweight="bold", color=focal_c, linespacing=1.1, zorder=6,
             bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.9))
# Node2Vec -> embedding -> score
axA.annotate("", xy=(2.55, 0.08), xytext=(1.25, 0.08), arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.3))
axA.text(1.9, 0.16, "Node2Vec", ha="center", va="bottom", fontsize=6.3, fontweight="bold")
rng = np.random.default_rng(7); bw = 0.05                    # illustrative bars only
for row, (y0, lab) in enumerate(((0.6, "aHUS"), (-0.6, "PNH"))):
    v = rng.normal(0, 1, 14); v = v / np.abs(v).max() * 0.28
    for i_, x in enumerate(v):
        axA.add_patch(mpatches.Rectangle((3.1 + i_ * (bw + 0.01), y0 + min(0, x)), bw, max(abs(x), 0.01),
                                         fc=focal_c, ec="none", alpha=0.85))
    axA.plot([3.05, 3.1 + 14 * (bw + 0.01)], [y0, y0], color="#999999", lw=0.5)
    axA.text(3.02, y0, lab, ha="right", va="center", fontsize=5.6, color="#222222")
axA.text(3.52, -1.2, "DisAtlas embeddings\n(32-d, 5-seed ensemble;\nbars illustrative)", ha="center", va="top",
         fontsize=5.6, color="#222222")
axA.text(3.52, 1.12, f"Model score p = {EX['p']:.2f}\n(top {100 * (1 - EX['frac_random_pairs_below_p']):.1f}% "
         "of random disease pairs)", ha="center", va="bottom", fontsize=5.8, color=focal_c, fontweight="bold")
axA.set_title("Real example: aHUS and PNH in the MechDist training graph", loc="left", fontsize=8)
axA.text(-2.9, -2.05, "Edges: shared HPO phenotype and/or HPA/Pfam protein evidence (rank-mean fused).\n"
         f"12th shared neighbour not drawn: 'neoplasm' umbrella term (degree {NB[HUB]['train_degree']:,}).",
         ha="left", va="top", fontsize=5.4, color="#555555")

# ---------------------------------------------------------------- panel b
bins = np.linspace(0, 1.75, 36)
pos_c, neg_c = muted_focal_pair("#c0392b")
axB.hist(dneg, bins=bins, color=neg_c, alpha=0.6, label="IxIDN\u207b (negative)", density=True)
axB.hist(dpos, bins=bins, color=pos_c, alpha=0.65, label="IxIDN\u207a (positive)", density=True)
axB.axvline(np.median(dneg), color=neg_c, lw=1.1, ls="--"); axB.axvline(np.median(dpos), color=pos_c, lw=1.1, ls="--")
axB.set_xlabel("Euclidean distance between the two diseases'\nDisAtlas embeddings (32-d, ensembled over 5 seeds)")
axB.set_ylabel("Density")
axB.set_title("IxIDN-Matched pairs separate by embedding distance", loc="left", fontsize=8)
axB.legend(frameon=False, fontsize=6.3, loc="upper left", handlelength=1.2, handletextpad=0.5, borderaxespad=0.1)
axB.text(0.02, 0.80, f"n = {len(dpos):,} pairs each", transform=axB.transAxes, ha="left", va="top", fontsize=5.8,
         color="#555555")
axB.margins(x=0.02)
for side in ("bottom", "left"): axB.spines[side].set_linewidth(0.6)
axB.tick_params(direction="out", length=3, width=0.6)
panel_letter(axA, "a"); panel_letter(axB, "b")
print(save(fig, "mechdist_schematic", bbox_inches="tight"))
