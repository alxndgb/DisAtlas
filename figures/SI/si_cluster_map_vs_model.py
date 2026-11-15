"""SI figure: distance on the two-dimensional disease map vs model distance (#fig:SI_cluster_map_vs_model).
Pairs: 2,000,000 random index pairs (numpy default_rng(0), pairs with i = j dropped) over the 13,021 diseases; model
distance = 1 - cosine of the five-seed ensemble (data/zenodo/model_outputs/cosine_matrix), map distance = Euclidean
distance on the Isomap layout (cosine, 15 neighbours; isomap_x/y of data/zenodo/model_outputs/clusters/
disease_clusters_k55.parquet), both stored as float32 as in the analysis run. Spearman rho on all pairs.
Usage: python figures/SI/si_cluster_map_vs_model.py
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import sys, os, json
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import pandas as pd
from disatlas import figstyle_si as fs
import disatlas.deterministic  # noqa: F401
MO = REPO / "data" / "zenodo" / "model_outputs"; VJ, OUT = DERIVED / "si_cluster" / "sifu_cluster_13021_values.json", HERE; fs.apply_figure_style()
ids = json.load(open(MO / "cosine_matrix" / "ids_13021.json")); N13 = len(ids); assert N13 == 13021
C = np.load(MO / "cosine_matrix" / "cosine_ens5_13021.npy", mmap_mode="r")
XY = pd.read_parquet(MO / "clusters" / "disease_clusters_k55.parquet").set_index("disease_id").loc[ids, ["isomap_x", "isomap_y"]].to_numpy(float)
rg = np.random.default_rng(0); i_, j_ = rg.integers(0, N13, 2_000_000), rg.integers(0, N13, 2_000_000); k_ = i_ != j_; i_, j_ = i_[k_], j_[k_]
D = 1.0 - np.asarray(C, dtype=np.float32); np.fill_diagonal(D, 0.0)
dmod = D[i_, j_].astype(np.float32).astype(float); del D
dmap = np.linalg.norm(XY[i_] - XY[j_], axis=1).astype(np.float32).astype(float)
V = json.load(open(VJ))["layout"]; rho = spearmanr(dmod, dmap)[0]
assert len(dmod) == V["n_pairs"] and abs(rho - V["spearman"]) < 1e-4, (len(dmod), rho, V)
xmax = float(np.ceil(np.percentile(dmod, 99.9) * 10) / 10)  # cosine distance can exceed 1
fig, ax = plt.subplots(figsize=(3.6, 3.2))
hb = ax.hexbin(dmod, dmap / np.percentile(dmap, 99.9), gridsize=60, bins="log", cmap="Blues", mincnt=1, linewidths=0, extent=(0, xmax, 0, 1))
q = np.linspace(0, 1, 21); bx = np.quantile(dmod, q); mids, meds = [], []
for lo, hi in zip(bx[:-1], bx[1:]):
    m = (dmod >= lo) & (dmod <= hi); mids.append(np.median(dmod[m])); meds.append(np.median(dmap[m]) / np.percentile(dmap, 99.9))
ax.plot(mids, meds, "-", color="#D55E00", lw=1.3, label="median map distance\n(per 5% of pairs)")
ax.set_xlabel("Model distance (cosine)"); ax.set_ylabel("Map distance (scaled)"); ax.set_xlim(0, xmax); ax.set_ylim(0, 1)
ax.set_title(f"Spearman ρ = {rho:.2f}, {len(dmod):,} random pairs")
ax.legend(loc="upper left", fontsize=6)
cb = fig.colorbar(hb, ax=ax, fraction=.046, pad=.03); cb.set_label("Pairs per bin")
N = dict(png_bytes=fs.save(fig, OUT, "SI_cluster_map_vs_model_distance"), overlaps=fs.overlap_check(fig), rho=rho, n_pairs=int(len(dmod)),
         xmax=xmax, share_dmod_above_xmax=float(np.mean(dmod > xmax)))
json.dump(N, open(f"{OUT}/SI_cluster_map_vs_model_distance_values.json", "w"), indent=1, default=float); print(N)
