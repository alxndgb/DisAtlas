"""Appendix figure #fig:SI_cluster_validation_layout (display layout fidelity, Spearman rho = 0.50).

Hexagonal-binned UPNA-DDAL distance (1 - p) against the distance in the two-dimensional display layout of
Figure 3 for random disease pairs. Pairs are drawn exactly as in the layout-quality check: 300,000 index pairs
from numpy default_rng(0) over the 13,019 diseases in display-map order (self-pairs dropped); map distance is
rescaled by the least-squares factor s; 20,000 of the pairs are shown.

Inputs
    figures/main/data/fig3_display_map.tsv.gz              display-layout coordinates (disease_id, cluster, x, y)
    source_data/layout_fidelity_random_pairs.tsv.gz        the drawn pairs (a, b) and their UPNA-DDAL distance d = 1 - p
Options
    --from-model   recompute d with the frozen model package in data/zenodo/models for every pair whose two
                   diseases carry the same identifier in the package and in the display map, and check it

    python figures/SI/si_layout_fidelity.py [--from-model]
"""
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import REPO_ROOT, src, save

dm = pd.read_csv(os.path.join(REPO_ROOT, "figures/main/data/fig3_display_map.tsv.gz"), sep="\t")
T = pd.read_csv(src("layout_fidelity_random_pairs.tsv.gz"), sep="\t", float_precision="round_trip")
ids = dm.disease_id.to_numpy(); n = len(ids)
rng = np.random.default_rng(0)
P = rng.integers(0, n, (300_000, 2)); P = P[P[:, 0] != P[:, 1]]
assert (T.a.to_numpy() == ids[P[:, 0]]).all() and (T.b.to_numpy() == ids[P[:, 1]]).all(), "pair draw differs"

if "--from-model" in sys.argv:
    pkg = os.path.join(REPO_ROOT, "data/zenodo/models"); sys.path.insert(0, pkg)
    from frozen_model_loader import FrozenMechDistModel
    m = FrozenMechDistModel(pkg); known = set(m.ids)
    ok = (T.a.isin(known) & T.b.isin(known)).to_numpy()
    p = np.asarray(m.score_pairs(list(T.a[ok]), list(T.b[ok])), dtype=np.float64)
    err = np.abs((1 - p) - T.d.to_numpy()[ok])
    assert err.max() < 1e-6, err.max()
    print(f"from model: {ok.sum():,} of {len(T):,} pairs recomputed, max |diff| = {err.max():.1e}")

Y = dm[["x", "y"]].to_numpy()
dt = T.d.to_numpy(np.float64)
d2 = np.linalg.norm(Y[P[:, 0]] - Y[P[:, 1]], axis=1).astype(np.float64)
s = (dt @ d2) / (d2 @ d2)
rho = float(np.corrcoef(pd.Series(dt).rank(), pd.Series(d2).rank())[0, 1])
assert f"{rho:.2f}" == "0.50", rho
rng.choice(n, 3000, replace=False); rng.choice(n, 4000, replace=False)   # kNN / trustworthiness draws of the check
sh = rng.choice(len(dt), 20_000, replace=False)
S = np.c_[dt[sh], s * d2[sh]].astype(np.float32)
print(f"OK: Spearman rho = {rho:.4f} over {len(dt):,} pairs; scale s = {s:.4f}")

mpl.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 300})
FOC = "#1f3b73"
fig = plt.figure(figsize=(11.0, 5.2))          # geometry of the original multi-panel layout
ax = fig.add_axes([0.06, 0.10, 0.25, 0.335])
ax.hexbin(S[:, 0], S[:, 1], gridsize=40, bins="log", cmap="Greys", mincnt=1, linewidths=0, extent=(0, 1.02, 0, 1.6))
ax.plot([0, 1.02], [0, 1.02], color=FOC, lw=0.8); ax.set_xlim(-0.02, 1.04)
ax.set_xlabel("UPNA-DDAL distance (1 − p)"); ax.set_ylabel("Map distance (rescaled)")
ax.set_title("Map vs model distance", loc="left")
ax.text(0.03, 0.86, f"Spearman ρ = {rho:.2f}", transform=ax.transAxes, fontsize=7)
# The appendix image is a crop of a larger multi-panel canvas; this window reproduces that crop
# (including its sub-pixel placement) on the same 11.0 x 5.2 in canvas.
from matplotlib.transforms import Bbox
CROP = Bbox([[0.019443, 0.028257], [3.502777, 2.578257]])
print(save(fig, "SI_cluster_map_vs_model_distance", bbox_inches=CROP, pad_inches=0))
