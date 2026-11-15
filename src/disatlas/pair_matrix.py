"""All-pairs cosine matrix of the 13,021 embedded diseases (five-seed mean cosine; distance = 1 - cosine).

    python -m disatlas.pair_matrix build <out_dir> [--model DIR]   # cosine_ens5_13021.npy + ids_13021.json
    python -m disatlas.pair_matrix compare <out_dir>                # against data/zenodo/model_outputs/cosine_matrix

``build`` writes the n x n float32 matrix (rows/columns in the order of ``ids_13021.json``, the order of the model's
embedding files) in the layout of the shipped data/zenodo/model_outputs/cosine_matrix/cosine_ens5_13021.npy.
``compare`` reports the largest absolute difference on the diseases both matrices hold; the shipped matrix is
reproduced from the shipped embeddings to float32 rounding (|diff| < 1e-6). For a cold retrain (P_full graph, a
different disease set and new embeddings) the difference is informational: the retrain is judged by the parity rule of
tests/test_cold_tolerance.py, not by score identity.
"""
import argparse
import glob
import json
import os

import numpy as np

from disatlas import paths

SHIPPED = paths.ZENODO / "model_outputs" / "cosine_matrix"


def build(out: str, model_dir=None):
    from disatlas.model import load_ensemble
    m = load_ensemble(model_dir)
    os.makedirs(out, exist_ok=True)
    C = m.matrix()
    np.save(os.path.join(out, f"cosine_ens5_{len(m.ids)}.npy"), C)
    json.dump(m.ids, open(os.path.join(out, f"ids_{len(m.ids)}.json"), "w"))
    return m.ids, C


def compare(out: str, ref=SHIPPED, block: int = 1024) -> dict:
    def one(d, pat):
        f = sorted(glob.glob(os.path.join(d, pat)))
        assert len(f) == 1, (d, pat, f)
        return f[0]
    ids_o = json.load(open(one(out, "ids_*.json")))
    ids_r = json.load(open(one(ref, "ids_*.json")))
    A = np.load(one(out, "cosine_ens5_*.npy"), mmap_mode="r")
    B = np.load(one(ref, "cosine_ens5_*.npy"), mmap_mode="r")
    if ids_o != ids_r:
        p = {d: i for i, d in enumerate(ids_o)}
        common = [d for d in ids_r if d in p]
        q = {d: i for i, d in enumerate(ids_r)}
        ia, ib = np.array([p[d] for d in common]), np.array([q[d] for d in common])
    else:
        ia = ib = np.arange(len(ids_r))
    mx = 0.0
    for i in range(0, len(ia), block):
        a = np.asarray(A[ia[i:i + block]][:, ia], np.float64); b = np.asarray(B[ib[i:i + block]][:, ib], np.float64)
        mx = max(mx, float(np.abs(a - b).max()))
    return dict(n=int(len(ia)), same_order=ids_o == ids_r, max_abs_diff=mx)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["build", "compare"])
    ap.add_argument("out")
    ap.add_argument("--model", default=None, help="model directory (default: disatlas.model.final_model_dir())")
    a = ap.parse_args()
    if a.cmd == "build":
        build(a.out, a.model)
    if not (SHIPPED / "cosine_ens5_13021.npy").exists():
        print("data/zenodo/model_outputs/cosine_matrix absent: no shipped matrix to compare with"); return
    print(json.dumps(compare(a.out)))


if __name__ == "__main__":
    main()
