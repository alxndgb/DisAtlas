"""Rebuild data/derived/si_cluster/k_selection_boot.tsv.gz (Figure S16 a: K selection by bootstrap silhouette).

    python -m disatlas.rebuild.k_selection_boot OUTDIR [REP ...]     (default: all 100 replicates; hours on 8 CPUs)
    sbatch scripts/slurm/k_selection_boot_array.sh OUTDIR; ... OUTDIR --merge      (150-task array, about 10 min each)

Needs the optional cluster extra:  .venv/bin/pip install -r requirements-cluster.txt   (kmedoids, FasterPAM)
Design (the authors' select_k_bootstrap, cosine kind): D = max(1 - cosine, 0), diagonal 0, float64, from
data/zenodo/model_outputs/cosine_matrix/cosine_ens5_13021.npy. Replicate r: 80 % random subsample of the 13,021
diseases without replacement (numpy default_rng(1000 + r), sorted), k-medoids FasterPAM at K = 10, 15, ..., 80
(random_state = r), exact silhouette on the subsample (precomputed distance), smallest cluster size, adjusted Rand
index against the full-data FasterPAM labels (random_state = 0) restricted to the same diseases.
Test: tests/test_rebuild_k_selection_boot.py (replicate 0 at every K).
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
MO = ROOT / "data/zenodo/model_outputs/cosine_matrix"
K_GRID = list(range(10, 85, 5)); N_BOOT = 100; FRAC = 0.8
NCPU = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))


def distance():
    ids = json.load(open(MO / "ids_13021.json"))
    D = 1.0 - np.load(MO / "cosine_ens5_13021.npy").astype(np.float64)
    D = np.maximum(D, 0.0); np.fill_diagonal(D, 0.0)
    return ids, D


def _km(D, k, seed):
    import kmedoids
    return np.asarray(kmedoids.fasterpam(D, k, random_state=seed, n_cpu=NCPU).labels, dtype=np.int32)


def compute(reps=range(N_BOOT), k_grid=K_GRID):
    from sklearn.metrics import silhouette_score, adjusted_rand_score
    ids, D = distance(); n = len(ids)
    full0 = {k: _km(D, k, 0) for k in k_grid}
    rows = []
    for r in reps:
        sub = np.sort(np.random.default_rng(1000 + r).choice(n, int(FRAC * n), replace=False))
        Ds = np.ascontiguousarray(D[np.ix_(sub, sub)])
        for k in k_grid:
            lab = _km(Ds, k, r)
            rows.append(dict(kind="boot", rep=r, K=k, silhouette=float(silhouette_score(Ds, lab, metric="precomputed")),
                             min_size=int(np.bincount(lab).min()), sub_seed=1000 + r,
                             ari_vs_full=float(adjusted_rand_score(full0[k][sub], lab))))
    return pd.DataFrame(rows)


def shard(task):
    """SLURM array task 0..149 -> (K, block of 10 replicates); writes one file per task (see scripts/slurm/)."""
    return K_GRID[task // 10], list(range((task % 10) * 10, (task % 10) * 10 + 10))


def merge(outdir):
    return pd.concat([pd.read_csv(f, sep="\t") for f in sorted(Path(outdir).glob("K*_r*.tsv"))]).sort_values(["rep", "K"]).reset_index(drop=True)


if __name__ == "__main__":
    o = Path(sys.argv[1]); o.mkdir(parents=True, exist_ok=True)
    if sys.argv[2:3] == ["--task"]:
        k, reps = shard(int(sys.argv[3])); compute(reps, [k]).to_csv(o / f"K{k}_r{reps[0]:03d}.tsv", sep="\t", index=False)
    elif sys.argv[2:3] == ["--merge"]:
        merge(o).to_csv(o / "k_selection_boot.tsv.gz", sep="\t", index=False)
    else:
        reps = [int(x) for x in sys.argv[2:]] or list(range(N_BOOT))
        compute(reps).to_csv(o / "k_selection_boot.tsv.gz", sep="\t", index=False)
