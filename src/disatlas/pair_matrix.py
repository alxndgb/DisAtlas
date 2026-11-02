"""All-pairs ensemble scores for the 13,021 embedded diseases.

    python -m disatlas.pair_matrix seed <seed> <tmpdir>     # one seed -> <tmpdir>/P_s<seed>.npy (float32, n x n)
    python -m disatlas.pair_matrix combine <tmpdir> <out.npz>

The output npz holds ``ids`` (sorted disease ids) and ``p_upper`` (float32), the strict upper triangle
(row-major, k=1) of the five-seed mean probability matrix; ``square(npz)`` restores the full symmetric
matrix with NaN diagonal.
"""
import os
import sys

import numpy as np


def seed_matrix(seed: int, out: str, model_dir=None, block: int = 512):
    from disatlas.model import DEFAULT_DIR, load_ensemble
    P = seed_P(load_ensemble(model_dir or DEFAULT_DIR), seed, block)
    os.makedirs(out, exist_ok=True)
    np.save(os.path.join(out, f"P_s{seed}.npy"), P)


def seed_P(m, seed: int, block: int = 512):
    import torch
    ids = m.ids
    X = m.emb[seed][[m.idx[seed][i] for i in ids]]
    n = len(ids)
    P = np.empty((n, n), np.float32)
    net = m.nets[seed]
    torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    with torch.no_grad():
        Xt = torch.as_tensor(X)
        for i in range(0, n, block):
            F = (Xt[i:i + block, None, :] * Xt[None, :, :]).reshape(-1, X.shape[1])
            P[i:i + block] = torch.sigmoid(net(F)).reshape(-1, n).numpy()
    return P


def combine(tmp: str, out: str, seeds=(4707, 4708, 4709, 4710, 4711)):
    from disatlas.model import DEFAULT_DIR, load_ensemble
    ids = load_ensemble(DEFAULT_DIR, verify=False).ids
    acc = None
    for s in seeds:
        P = np.load(os.path.join(tmp, f"P_s{s}.npy"), mmap_mode="r")
        acc = P.astype(np.float64) if acc is None else acc + P
    acc /= len(seeds)
    iu = np.triu_indices(len(ids), k=1)
    np.savez_compressed(out, ids=np.array(ids), p_upper=acc[iu].astype(np.float32))


def square(npz) -> tuple:
    z = np.load(npz)
    n = len(z["ids"])
    M = np.full((n, n), np.nan, np.float32)
    iu = np.triu_indices(n, k=1)
    M[iu] = z["p_upper"]
    M.T[iu] = z["p_upper"]
    return list(z["ids"]), M


if __name__ == "__main__":
    if sys.argv[1] == "seed":
        seed_matrix(int(sys.argv[2]), sys.argv[3])
    else:
        combine(sys.argv[2], sys.argv[3])
