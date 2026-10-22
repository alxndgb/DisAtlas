"""Optional: recompute the saved TxGNN pair scores from the published checkpoint (Protocol 4).

Not needed to reproduce the paper: score_txgnn.py reads data/txgnn_pair_scores.tsv.gz. Run
fetch_txgnn.sh first, then

    python benchmarks/txgnn/run_txgnn_inference.py      # writes benchmarks/txgnn/rerun/txgnn_pair_scores.tsv.gz

Needs torch and pandas; CPU is enough. The checkpoint is used unchanged (not retrained on MechDist or on
IxIDN-Matched): released disease embeddings (node_emb.pkl, rows = TxGNN disease index) and relation weights
(model.pt, w_rels). Score = TxGNN's disease-disease DistMult, sigmoid(sum h_a * w_dd * h_b).
MONDO ids map to TxGNN disease nodes through the PrimeKG full_graph split (grouped node ids such as
'123_456' map every part). Pairs whose two diseases fall on the same TxGNN node are flagged `collapsed`;
unmapped and collapsed pairs are ranked last in scoring.
"""
import gzip
import os
import pickle
import sys

import numpy as np
import pandas as pd
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from _common import load_key  # noqa: E402

X = os.path.join(HERE, "external")
CKPT = os.path.join(X, "checkpoint")
SPLIT = os.path.join(X, "primekg", "full_graph_42")
ET = ("disease", "disease_disease", "disease")


def idmap(df, t):
    m = pd.concat([df.loc[df.x_type == t, ["x_id", "x_idx"]].set_axis(["id", "idx"], axis=1),
                   df.loc[df.y_type == t, ["y_id", "y_idx"]].set_axis(["id", "idx"], axis=1)]).drop_duplicates()
    return dict(zip(m.id.astype(str), m.idx.astype(int)))


def mondo_lookup(dis):
    out = {}
    for nid, idx in dis.items():
        for part in nid.split("_"):
            try:
                out.setdefault(str(int(float(part))), set()).add(idx)
            except ValueError:
                pass
    return {k: sorted(v) for k, v in out.items()}


def map_pairs(pairs, lookup):
    def nodes_of(curie):
        prefix, code = curie.split(":", 1)
        return lookup.get(str(int(code)), []) if prefix == "MONDO" else []
    na, nb = pairs.a.map(nodes_of), pairs.b.map(nodes_of)
    a, b = [], []
    for x, y in zip(na, nb):  # first node per side; a collapsed pair uses the shared node on both sides
        shared = sorted(set(x) & set(y))
        a.append(shared[0] if shared else (x[0] if x else -1))
        b.append(shared[0] if shared else (y[0] if y else -1))
    return np.array(a), np.array(b), np.array([bool(set(x) & set(y)) for x, y in zip(na, nb)])


def distmult(ea, w, eb):
    return 1 / (1 + np.exp(-(ea * w * eb).sum(1)))


def main():
    df = pd.concat([pd.read_csv(os.path.join(SPLIT, f"{s}.csv")) for s in ["train", "valid", "test"]])
    dis, drug = idmap(df, "disease"), idmap(df, "drug")
    drug_id = pickle.load(open(os.path.join(CKPT, "drug_id.pkl"), "rb"))
    assert all(drug.get(str(x)) == i for i, x in enumerate(drug_id)), "node indexing differs from the published model"
    sd = torch.load(os.path.join(CKPT, "model.pt"), map_location="cpu")
    etypes = sorted({tuple(r) for r in df[["x_type", "relation", "y_type"]].drop_duplicates().values})
    assert len(etypes) == sd["w_rels"].shape[0]
    W = sd["w_rels"].numpy()
    E = pickle.load(open(os.path.join(CKPT, "node_emb.pkl"), "rb"))["disease"].numpy()
    assert E.shape[0] == len(set(dis.values()))

    pairs = load_key()[["blind_id", "a", "b"]]
    ia, ib, collapsed = map_pairs(pairs, mondo_lookup(dis))
    ok = (ia >= 0) & (ib >= 0)
    ea, eb = E[np.where(ok, ia, 0)], E[np.where(ok, ib, 0)]
    ok &= (np.abs(ea).sum(1) > 0) & (np.abs(eb).sum(1) > 0)
    score = np.where(ok, distmult(ea, W[etypes.index(ET)], eb), np.nan)
    cos = np.where(ok, (ea * eb).sum(1) / (np.linalg.norm(ea, axis=1) * np.linalg.norm(eb, axis=1) + 1e-12), np.nan)
    out = pairs.assign(idx_a=ia, idx_b=ib, mapped=ok, collapsed=collapsed, score=score, cos=cos)
    os.makedirs(os.path.join(HERE, "rerun"), exist_ok=True)
    p = os.path.join(HERE, "rerun", "txgnn_pair_scores.tsv.gz")
    with gzip.GzipFile(p, "wb", mtime=0) as g:
        g.write(out.to_csv(sep="\t", index=False, lineterminator="\n").encode())
    saved = pd.read_csv(os.path.join(HERE, "data", "txgnn_pair_scores.tsv.gz"), sep="\t")
    same = np.allclose(saved.score.values, out.score.values, rtol=1e-5, equal_nan=True)
    print(p, "matches saved scores" if same else "differs from saved scores")


if __name__ == "__main__":
    main()
