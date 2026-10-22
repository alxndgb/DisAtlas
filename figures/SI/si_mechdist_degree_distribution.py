"""Supplementary Figure SI_mechdist_degree_dist: degree distribution of the MechDist disease-disease network.

Appendix section: IxIDN-Matched benchmark construction (#sec:SI_ixidn_matched_build), Figure #fig:SI_mechdist_degree_dist.
  Log-log scatter of node degree against the number of diseases with that degree, annotated with the network
  statistics quoted in the caption (N = 13,021 diseases, L = 332,082 edges, mean degree 51.0, average clustering
  coefficient 0.47).
The statistics are recomputed from the shipped network edge and node lists and checked against the caption.
Uses its own light plot style (kept to match the manuscript image), not the SI house style.
Inputs: data/03_network/mechdist_network_edges.parquet.gzip, data/03_network/mechdist_network_nodes.parquet.gzip
Output: figures/SI/mechdist_degree_distribution.{png,pdf}
Run:    python figures/SI/si_mechdist_degree_distribution.py      (needs networkx; about 30 s for square clustering)
"""
import os, sys
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _style import save  # noqa: E402

NET = os.path.join(HERE, "..", "..", "data", "03_network")
PAL = dict(blue="#2a78d6", ink="#0b0b0b", ink2="#52514e", muted="#898781", grid="#e1e0d9", axis="#c3c2b7",
           surface="#fcfcfb")

# ---- network and statistics ----
E = pd.read_parquet(os.path.join(NET, "mechdist_network_edges.parquet.gzip"), columns=["disease_a", "disease_b"])
V = pd.read_parquet(os.path.join(NET, "mechdist_network_nodes.parquet.gzip"), columns=["disease_id"])
G = nx.Graph()
G.add_nodes_from(V.disease_id)
G.add_edges_from(zip(E.disease_a, E.disease_b))
degrees = np.array([d for _, d in G.degree()])
N, L = G.number_of_nodes(), G.number_of_edges()
mean_degree = float(degrees.mean())
avg_clustering = nx.average_clustering(G)
avg_square_clustering = float(np.mean(list(nx.square_clustering(G).values())))

# ---- value check: caption of #fig:SI_mechdist_degree_dist ----
assert (N, L) == (13021, 332082), (N, L)
assert round(mean_degree, 1) == 51.0 and round(avg_clustering, 2) == 0.47, (mean_degree, avg_clustering)
assert (int(degrees.min()), int(degrees.max())) == (1, 3123)
print("caption values: OK")

# ---- figure ----
plt.rcParams.update({
    "figure.facecolor": PAL["surface"], "axes.facecolor": PAL["surface"],
    "axes.edgecolor": PAL["axis"], "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": PAL["grid"], "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "text.color": PAL["ink"], "axes.labelcolor": PAL["ink2"], "axes.titlecolor": PAL["ink"],
    "xtick.color": PAL["muted"], "ytick.color": PAL["muted"], "font.size": 11, "figure.dpi": 110,
})
fig, ax = plt.subplots(figsize=(5.4, 4.6))
vals, counts = np.unique(degrees, return_counts=True)
ax.scatter(vals, counts, s=12, color=PAL["blue"], zorder=3)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("degree"); ax.set_ylabel("frequency"); ax.set_title("MechDist Network Degree Distribution")
ax.text(0.04, 0.22,
        f"N: {N:,}  L: {L:,}  <k>: {mean_degree:.3f}\n"
        f"min/max-deg: {degrees.min()}/{degrees.max()}\n"
        f"Avg clustering coefficient: {avg_clustering:.2f}\n"
        f"Avg square clustering coefficient: {avg_square_clustering:.2f}",
        transform=ax.transAxes, ha="left", va="top", fontsize=9.5, color=PAL["ink2"])
fig.tight_layout()
print(save(fig, "mechdist_degree_distribution"))
