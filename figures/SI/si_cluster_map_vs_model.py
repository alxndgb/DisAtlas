"""SI figure: distance on the two-dimensional disease map vs model distance (#fig:SI_cluster_map_vs_model).
Pairs: the random disease pairs saved by sifu_cluster_13021.py (13,021-disease index; Isomap, cosine, 15 neighbours,
as in the main-text map); model distance = 1 - p (5-seed mean decoder probability). Spearman rho on all pairs.
Usage: si_cluster_map_vs_model.py <layout_pairs_13021.npz> <values.json> <outdir>
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import sys, os, json
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from disatlas import figstyle_si as fs
NPZ, VJ, OUT = ZENODO_SIBIO / "layout_pairs_13021.npz", DERIVED / "si_cluster" / "sifu_cluster_13021_values.json", HERE; fs.apply_figure_style()
z = np.load(NPZ); dmod, dmap = z["dmod"].astype(float), z["dmap"].astype(float)
V = json.load(open(VJ))["layout"]; rho = spearmanr(dmod, dmap)[0]
assert len(dmod) == V["n_pairs"] and abs(rho - V["spearman"]) < 1e-4, (len(dmod), rho, V)
fig, ax = plt.subplots(figsize=(3.6, 3.2))
hb = ax.hexbin(dmod, dmap / np.percentile(dmap, 99.9), gridsize=60, bins="log", cmap="Blues", mincnt=1, linewidths=0, extent=(0, 1, 0, 1))
q = np.linspace(0, 1, 21); bx = np.quantile(dmod, q); mids, meds = [], []
for lo, hi in zip(bx[:-1], bx[1:]):
    m = (dmod >= lo) & (dmod <= hi); mids.append(np.median(dmod[m])); meds.append(np.median(dmap[m]) / np.percentile(dmap, 99.9))
ax.plot(mids, meds, "-", color="#D55E00", lw=1.3, label="median map distance\n(per 5% of pairs)")
ax.set_xlabel("Model distance (1 − p)"); ax.set_ylabel("Map distance (scaled)"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_title(f"Spearman ρ = {rho:.2f}, {len(dmod):,} random pairs")
ax.legend(loc="upper left", fontsize=6)
cb = fig.colorbar(hb, ax=ax, fraction=.046, pad=.03); cb.set_label("Pairs per bin")
N = dict(png_bytes=fs.save(fig, OUT, "SI_cluster_map_vs_model_distance"), overlaps=fs.overlap_check(fig), rho=rho, n_pairs=int(len(dmod)))
json.dump(N, open(f"{OUT}/SI_cluster_map_vs_model_distance_values.json", "w"), indent=1, default=float); print(N)
