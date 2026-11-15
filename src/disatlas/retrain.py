"""Optional retrain path for the DisAtlas ensemble (the recipe of the frozen run; not needed to reproduce results).

The shipped ensemble (data/zenodo/model) is the model of record. This module re-runs the training recipe that
produced it, step for step and with the same random streams:

  training graph (build/graph/pfull_graph_edges.parquet from python -m disatlas.graph.pfull_build, file row
  order) minus every edge that is an IxIDN+ or IxIDN-
  pair -> weighted second-order random walks on w_arm (exact rejection sampling, p = 1, q = 0.5, walks of 20
  nodes, 20 per start node, regenerated every batch; start nodes drawn per epoch with probability ~ deg^-0.5)
  -> skip-gram on windows of 3 nodes with one embedding table (dim 32, Adam lr 0.01, 25 epochs, 128 start nodes
  per batch; 5 negative walks per positive walk, nodes drawn ~ deg^1) -> rows L2 unit-normalised, stored float32.
  Seeds 4707-4711; the pair score is the mean over seeds of the cosine of the two disease embeddings
  (disatlas.model). No decoder is trained.

The frozen walker was compiled with numba; its random stream is the Mersenne Twister of np.random.RandomState,
which this pure-python walker reads in the same order, so walks, splits and initial weights equal the frozen
ones. Floating-point sums in torch depend on the thread count and on the BLAS build, so retrained embeddings match
the frozen ones to about 2e-7, not bit for bit. Each seed reseeds every random stream it uses (torch.manual_seed,
its own torch.Generator and RandomState) before training, so the embedding of a seed does not depend on what ran
before it in the same process.

Needs torch (requirements-retrain.txt; not part of the warm environment).

    python -m disatlas.retrain --out runs/retrain --seeds 4707 4708 4709 4710 4711
    python -m disatlas.retrain --out runs/smoke --seeds 4707 --epochs 1 --num-walks 2
"""
from __future__ import annotations

import argparse
import bisect
import json
import os

import numpy as np
import pandas as pd

from disatlas.paths import BUILD, RAW, processed

RECIPE = dict(dim=32, walk_length=20, num_walks=20, window=2, p=1.0, q=0.5, epochs=25, lr=0.01, neg=5,
              neg_alpha=1.0, start_beta=0.5, batch_size=128)


def load_graph(union_path, ixidn_dir):
    """Training graph as an adjacency dict in networkx insertion order (rows in file order, IxIDN edges removed).
    Returns (adj, excluded pairs). Duplicate rows add up their walk weight, as in the frozen build."""
    E = pd.read_parquet(union_path, columns=["u", "v", "w_arm"])
    adj = {}
    for u, v, w in zip(E.u.tolist(), E.v.tolist(), E.w_arm.tolist()):
        adj.setdefault(u, {})
        adj.setdefault(v, {})
        adj[u][v] = adj[v][u] = adj[u].get(v, 0.0) + w
    ex = set()
    for f in ("ixidn_plus_pairs.parquet", "ixidn_minus_pairs.parquet"):
        P = pd.read_parquet(os.path.join(ixidn_dir, f))
        ex |= {frozenset(p) for p in zip(P.disease_a, P.disease_b)}
    for p in ex:
        a, b = tuple(p)
        if a in adj and b in adj[a]:
            del adj[a][b], adj[b][a]
    return adj, ex


def edges_in_order(adj):
    """Edges in the order networkx Graph.edges() yields them."""
    seen, out = set(), []
    for u, nb in adj.items():
        out += [(u, v) for v in nb if v not in seen]
        seen.add(u)
    return out


class Walker:
    """Weighted node2vec walk: propose x ~ w(v, x), accept with alpha / max(1/p, 1, 1/q), alpha = 1/p if x is the
    previous node, 1 if x is a neighbour of it, 1/q otherwise. One RandomState stream for the whole run."""

    def __init__(self, adj, ids, p, q, seed):
        ix = {k: i for i, k in enumerate(ids)}
        self.nbr = [sorted(ix[v] for v in adj[k]) for k in ids]
        w = [adj[k][ids[j]] for k, nb in zip(ids, self.nbr) for j in nb]
        self.cw = np.r_[0.0, np.cumsum(w)].tolist()  # cumulative weights over the sorted CSR, as in the frozen run
        self.lo = np.r_[0, np.cumsum([len(nb) for nb in self.nbr])].tolist()
        self.nset = [set(nb) for nb in self.nbr]
        self.ip, self.iq = 1.0 / p, 1.0 / q
        self.amax = max(self.ip, 1.0, self.iq)
        self.rs, self.buf, self.k = np.random.RandomState(seed), [], 0

    def _r(self):
        if self.k == len(self.buf):
            self.buf, self.k = self.rs.random_sample(1 << 20).tolist(), 0
        self.k += 1
        return self.buf[self.k - 1]

    def _pick(self, v):
        lo, hi, cw = self.lo[v], self.lo[v] + len(self.nbr[v]), self.cw
        u = cw[lo] + self._r() * (cw[hi] - cw[lo])
        return self.nbr[v][min(bisect.bisect_right(cw, u, lo + 1, hi + 1) - lo - 1, hi - lo - 1)]

    def walks(self, starts, steps):
        out = []
        for v in starts:
            w = [v, self._pick(v)]
            for _ in range(steps - 1):
                t, v = w[-2], w[-1]
                while True:
                    x = self._pick(v)
                    a = self.ip if x == t else (1.0 if x in self.nset[t] else self.iq)
                    if self._r() * self.amax < a:
                        break
                w.append(x)
            out.append(w)
        return out


def windows(rw, c):
    import torch
    return torch.cat([rw[:, j:j + c] for j in range(rw.shape[1] - c + 1)], dim=0)


def embed(adj, ids, cfg, seed):
    """Skip-gram node2vec (the loss and sampling of torch_geometric Node2Vec). Returns float32 (n, dim)."""
    import torch
    torch.manual_seed(seed)
    n, L, c, nw = len(ids), cfg["walk_length"], cfg["window"] + 1, cfg["num_walks"]
    emb = torch.nn.Embedding(n, cfg["dim"])
    emb.reset_parameters()  # Node2Vec initialises its table twice
    opt = torch.optim.Adam(emb.parameters(), lr=cfg["lr"])
    gen = torch.Generator().manual_seed(seed)
    deg = np.array([len(adj[k]) for k in ids], float)
    neg_p = torch.as_tensor(deg ** cfg["neg_alpha"] / (deg ** cfg["neg_alpha"]).sum(), dtype=torch.float)
    start_p = torch.as_tensor(deg ** -cfg["start_beta"] / (deg ** -cfg["start_beta"]).sum(), dtype=torch.float)
    walker = Walker(adj, ids, cfg["p"], cfg["q"], seed)

    def nll(rw, sign):
        h0 = emb(rw[:, 0]).view(len(rw), 1, -1)
        h = emb(rw[:, 1:].contiguous().view(-1)).view(len(rw), -1, cfg["dim"])
        s = torch.sigmoid((h0 * h).sum(dim=-1).view(-1))
        return -torch.log((s if sign > 0 else 1 - s) + 1e-15).mean()

    for ep in range(cfg["epochs"]):
        perm = torch.multinomial(start_p, n, replacement=True, generator=gen)
        tot = 0.0
        for i in range(0, n, cfg["batch_size"]):
            b = perm[i:i + cfg["batch_size"]]
            pos = windows(torch.tensor(walker.walks(b.repeat(nw).tolist(), L - 1)), c)
            bb = b.repeat(nw * cfg["neg"])
            rw = torch.multinomial(neg_p, bb.numel() * L, replacement=True, generator=gen).view(-1, L)
            neg = windows(torch.cat([bb.view(-1, 1), rw], dim=-1), c)
            loss = nll(pos, 1) + nll(neg, -1)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            tot += loss.item() * len(b)
        print(f"seed {seed} epoch {ep + 1} loss {tot / n:.5f}", flush=True)
    X = emb.weight.detach().numpy().astype(np.float64)
    return (X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)).astype(np.float32)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph", default=str(BUILD / "graph" / "pfull_graph_edges.parquet"),
                    help="training graph (u, v, w_arm); default: the P_full graph from disatlas.graph.pfull_build")
    ap.add_argument("--ixidn", default=str(RAW / "benchmarks" / "ixidn"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[4707, 4708, 4709, 4710, 4711])
    ap.add_argument("--threads", type=int, default=None, help="torch threads (default: torch's choice)")
    for k in ("epochs", "num_walks"):
        ap.add_argument("--" + k.replace("_", "-"), type=int, default=None)
    ap.add_argument("--config-only", action="store_true", help="only (re)write model_config.json for the seeds in --out")
    a = ap.parse_args()
    import torch
    if a.threads:
        torch.set_num_threads(a.threads)
    cfg = dict(RECIPE)
    for k in ("epochs", "num_walks"):
        if getattr(a, k) is not None:
            cfg[k] = getattr(a, k)
    adj, ex = load_graph(a.graph, a.ixidn)
    ids = sorted(k for k in adj if adj[k])
    if a.config_only:
        return write_config(a.out, ids, cfg)
    os.makedirs(os.path.join(a.out, "embeddings"), exist_ok=True)
    for s in a.seeds:
        X = embed(adj, ids, cfg, s)
        np.savez(os.path.join(a.out, "embeddings", f"node2vec_s{s}.npz"), ids=np.array(ids), X=X)
        info = dict(seed=s, nodes=len(ids), threads=torch.get_num_threads(), recipe=cfg)
        json.dump(info, open(os.path.join(a.out, f"retrain_s{s}.json"), "w"), indent=1)
        print(json.dumps({k: v for k, v in info.items() if k != "recipe"}), flush=True)
    write_config(a.out, ids, cfg)


def write_config(out, ids, cfg):
    """model_config.json + nodes.txt so disatlas.model.load_ensemble(out) reads the retrained seeds."""
    from disatlas.model import _sha256
    seeds = sorted(int(f[10:-4]) for f in os.listdir(os.path.join(out, "embeddings")) if f.startswith("node2vec_s") and f.endswith(".npz"))
    files = [f"embeddings/node2vec_s{s}.npz" for s in seeds]
    open(os.path.join(out, "nodes.txt"), "w").write("\n".join(ids) + "\n")
    json.dump(dict(model="DisAtlas ensemble retrained with disatlas.retrain (node2vec embeddings; pair score = 5-seed mean cosine)", seeds=seeds, recipe=cfg,
                   files={f: _sha256(os.path.join(out, f)) for f in files}),
              open(os.path.join(out, "model_config.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
