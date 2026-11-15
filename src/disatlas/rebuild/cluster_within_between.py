"""Rebuild data/derived/si_cluster/within_between_K55_13021.tsv.gz and the "internal" block of
sifu_cluster_13021_values.json (Figure S16 b; K = 55 silhouette on the cosine distance D = 1 - cosine).

    python -m disatlas.rebuild.cluster_within_between OUTDIR

Inputs: data/zenodo/model_outputs/cosine_matrix/{cosine_ens5_13021.npy, ids_13021.json} and
data/zenodo/model_outputs/clusters/disease_clusters_k55.parquet (cluster of each disease).
within = mean distance to the other members of the own cluster; nearest_other = smallest mean distance to another
cluster; silhouette = (nearest_other - within) / max(within, nearest_other). Same arithmetic (float32 distance matrix,
float64 sums) as the authors' generator. Test: tests/test_rebuild_cluster_within_between.py.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
MO = ROOT / "data/zenodo/model_outputs"


def compute():
    ids = json.load(open(MO / "cosine_matrix/ids_13021.json")); N = len(ids); assert N == 13021
    K = pd.read_parquet(MO / "clusters/disease_clusters_k55.parquet").set_index("disease_id")
    lab = K.loc[ids, "cluster"].to_numpy(int); assert len(np.unique(lab)) == 55
    D = 1.0 - np.asarray(np.load(MO / "cosine_matrix/cosine_ens5_13021.npy", mmap_mode="r"), dtype=np.float32)
    np.fill_diagonal(D, 0.0)
    uniq = np.unique(lab); sums = np.zeros((N, len(uniq)))
    for j, c in enumerate(uniq):
        sums[:, j] = D[:, lab == c].sum(1)
    del D
    size = np.array([(lab == c).sum() for c in uniq]); own = np.searchsorted(uniq, lab)
    a = sums[np.arange(N), own] / np.maximum(size[own] - 1, 1)
    mo = sums / size[None, :]; mo[np.arange(N), own] = np.inf; b = mo.min(1)
    sil = np.where(size[own] > 1, (b - a) / np.maximum(a, b), 0.0)
    cs = pd.Series(sil).groupby(lab).mean()
    internal = dict(silhouette_all=float(sil.mean()), median_within=float(np.median(a)),
                    median_nearest_other=float(np.median(b)), share_diseases_negative_silhouette=float((sil < 0).mean()),
                    n_clusters_negative_mean_silhouette=int((cs < 0).sum()), min_cluster_size=int(size.min()),
                    max_cluster_size=int(size.max()))
    T = pd.DataFrame(dict(mondo_id=ids, cluster=lab, within=np.round(a, 5), nearest_other=np.round(b, 5),
                          silhouette=np.round(sil, 5)))
    return T, internal


if __name__ == "__main__":
    o = Path(sys.argv[1]); o.mkdir(parents=True, exist_ok=True)
    T, internal = compute()
    T.to_csv(o / "within_between_K55_13021.tsv.gz", sep="\t", index=False)
    json.dump(dict(internal=internal), open(o / "sifu_cluster_internal.json", "w"), indent=1)
