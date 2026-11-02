"""Optional retrain path for the DisAtlas ensemble (documented recipe; not needed to reproduce results).

The shipped ensemble (data/zenodo/model) is the model of record. This module re-runs the published recipe
on the shipped training graph so the training procedure can be inspected and re-executed:

  training graph (data/zenodo/graph/training_graph_edges.parquet; benchmark IxIDN+/IxIDN- pairs already
  removed) -> weighted second-order random walks on w_arm (20 walks per node, length 20, p = 1, q = 0.5)
  -> skip-gram with negative sampling (dim 32, window 2, 5 negatives, 25 epochs, lr 0.01, batch 128)
  -> isolated nodes re-attached at the origin, rows L2 unit-normalised
  -> Hadamard MLP decoder (hidden 128, dropout 0.3, lr 1e-3, weight decay 1e-4; positives = training
     edges, 5 uniform random negatives per positive; disease-disjoint 10 % validation split; early stopping)
  for seeds 4707-4711. Ensemble probability = mean over seeds of sigmoid(decoder(x_a * x_b)).

The frozen run used a compiled walker; this pure numpy/torch implementation follows the same recipe but
does not reproduce the frozen embeddings bit for bit (different random streams), so retrained scores agree
with the shipped ones statistically, not exactly.

    python -m disatlas.retrain --out runs/retrain --seeds 4707 4708 4709 4710 4711
    python -m disatlas.retrain --out runs/smoke --seeds 4707 --max-nodes 600 --epochs 1 --num-walks 2
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd

from disatlas.paths import processed

RECIPE = dict(dim=32, walk_length=20, num_walks=20, window=2, p=1.0, q=0.5, epochs=25, lr=0.01, neg=5,
              batch_size=128, dec_hidden=128, dec_dropout=0.3, dec_lr=1e-3, dec_wd=1e-4, dec_neg=5,
              dec_epochs=200, dec_patience=10, val_frac=0.1)


def load_graph(path, max_nodes=None, seed=0):
    E = pd.read_parquet(path)[["u", "v", "w_arm"]]
    if max_nodes:  # smoke test: induced subgraph on a BFS ball around a random node
        rng = np.random.default_rng(seed)
        adj = {}
        for a, b in zip(E.u, E.v):
            adj.setdefault(a, []).append(b)
            adj.setdefault(b, []).append(a)
        start = sorted(adj)[rng.integers(len(adj))]
        keep, frontier = {start}, [start]
        while frontier and len(keep) < max_nodes:
            nxt = []
            for x in frontier:
                for y in adj[x]:
                    if y not in keep and len(keep) < max_nodes:
                        keep.add(y)
                        nxt.append(y)
            frontier = nxt
        E = E[E.u.isin(keep) & E.v.isin(keep)]
    ids = sorted(set(E.u) | set(E.v))
    ix = {k: i for i, k in enumerate(ids)}
    u, v = E.u.map(ix).values, E.v.map(ix).values
    w = E.w_arm.values.astype(float)
    n = len(ids)
    src = np.concatenate([u, v])
    dst = np.concatenate([v, u])
    ww = np.concatenate([w, w])
    o = np.lexsort((dst, src))
    src, dst, ww = src[o], dst[o], ww[o]
    indptr = np.zeros(n + 1, np.int64)
    np.add.at(indptr, src + 1, 1)
    indptr = np.cumsum(indptr)
    return ids, indptr, dst, ww, np.stack([u, v], 1)


def walks(indptr, nbr, w, num_walks, length, p, q, rng):
    """Weighted node2vec walks (second order, unnormalised bias 1/p return, 1 common neighbour, 1/q out)."""
    n = len(indptr) - 1
    out = []
    nset = [set(nbr[indptr[i]:indptr[i + 1]].tolist()) for i in range(n)]
    for _ in range(num_walks):
        for s in rng.permutation(n):
            if indptr[s] == indptr[s + 1]:
                continue
            walk = [s]
            while len(walk) < length:
                cur = walk[-1]
                lo, hi = indptr[cur], indptr[cur + 1]
                nb, pw = nbr[lo:hi], w[lo:hi]
                if len(walk) == 1:
                    pr = pw
                else:
                    prev = walk[-2]
                    bias = np.where(nb == prev, 1.0 / p,
                                    np.where(np.fromiter((x in nset[prev] for x in nb), bool, len(nb)), 1.0, 1.0 / q))
                    pr = pw * bias
                walk.append(int(nb[rng.choice(len(nb), p=pr / pr.sum())]))
            out.append(walk)
    return out


def skipgram(walk_list, n, cfg, seed):
    import torch
    torch.manual_seed(seed)
    W = np.array([x for x in walk_list if len(x) == cfg["walk_length"]], np.int64)
    pairs = []
    for k in range(1, cfg["window"] + 1):
        pairs.append(np.stack([W[:, :-k].ravel(), W[:, k:].ravel()], 1))
    P = torch.as_tensor(np.concatenate(pairs))
    emb = torch.nn.Embedding(n, cfg["dim"], sparse=True)
    ctx = torch.nn.Embedding(n, cfg["dim"], sparse=True)
    torch.nn.init.normal_(emb.weight, std=0.1)
    torch.nn.init.zeros_(ctx.weight)
    opt = torch.optim.SparseAdam(list(emb.parameters()) + list(ctx.parameters()), lr=cfg["lr"])
    g = torch.Generator().manual_seed(seed)
    for _ in range(cfg["epochs"]):
        for b in torch.randperm(len(P), generator=g).split(cfg["batch_size"] * 64):
            a, c = P[b, 0], P[b, 1]
            neg = torch.randint(0, n, (len(b), cfg["neg"]), generator=g)
            ea = emb(a)
            pos = torch.nn.functional.logsigmoid((ea * ctx(c)).sum(-1))
            ng = torch.nn.functional.logsigmoid(-(ea[:, None, :] * ctx(neg)).sum(-1)).sum(-1)
            loss = -(pos + ng).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
    X = emb.weight.detach().numpy().astype(np.float32)
    return X


def finalize(X, has_edge):
    X = X.copy()
    X[~has_edge] = 0.0
    nrm = np.linalg.norm(X, axis=1, keepdims=True)
    X = np.where(nrm > 0, X / np.clip(nrm, 1e-12, None), 0.0).astype(np.float32)
    return X


def train_decoder(X, edges, cfg, seed):
    import torch
    from torch import nn
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    n = len(X)
    val_nodes = rng.random(n) < cfg["val_frac"]
    tr = edges[~(val_nodes[edges[:, 0]] | val_nodes[edges[:, 1]])]
    va = edges[val_nodes[edges[:, 0]] & val_nodes[edges[:, 1]]]
    if len(va) < 10:
        va = tr[rng.choice(len(tr), max(10, len(tr) // 10), replace=False)]
    net = nn.Sequential(nn.Linear(X.shape[1], cfg["dec_hidden"]), nn.ReLU(), nn.Dropout(cfg["dec_dropout"]),
                        nn.Linear(cfg["dec_hidden"], 1))
    opt = torch.optim.Adam(net.parameters(), lr=cfg["dec_lr"], weight_decay=cfg["dec_wd"])
    Xt = torch.as_tensor(X)
    lossf = nn.BCEWithLogitsLoss()

    def batch(E, r):
        neg = r.integers(0, n, (len(E) * cfg["dec_neg"], 2))
        A = np.concatenate([E[:, 0], neg[:, 0]])
        B = np.concatenate([E[:, 1], neg[:, 1]])
        y = np.concatenate([np.ones(len(E)), np.zeros(len(neg))]).astype(np.float32)
        return Xt[A] * Xt[B], torch.as_tensor(y)

    vF, vy = batch(va, np.random.default_rng(seed + 1))
    best, bad, state = np.inf, 0, None
    for _ in range(cfg["dec_epochs"]):
        net.train()
        for b in np.array_split(rng.permutation(len(tr)), max(1, len(tr) // 4096)):
            F, y = batch(tr[b], rng)
            opt.zero_grad()
            lossf(net(F).squeeze(1), y).backward()
            opt.step()
        net.eval()
        with torch.no_grad():
            vl = float(lossf(net(vF).squeeze(1), vy))
        if vl < best - 1e-4:
            best, bad, state = vl, 0, {k: t.clone() for k, t in net.state_dict().items()}
        else:
            bad += 1
            if bad >= cfg["dec_patience"]:
                break
    net.load_state_dict(state)
    return net.eval(), best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph", default=str(processed("graph", "training_graph_edges.parquet")))
    ap.add_argument("--out", required=True)
    ap.add_argument("--seeds", type=int, nargs="+", default=[4707, 4708, 4709, 4710, 4711])
    ap.add_argument("--max-nodes", type=int, default=None, help="smoke test on a BFS subgraph")
    for k in ("epochs", "num_walks", "dec_epochs"):
        ap.add_argument("--" + k.replace("_", "-"), type=int, default=None)
    ap.add_argument("--config-only", action="store_true", help="only (re)write model_config.json for the seeds in --out")
    a = ap.parse_args()
    cfg = dict(RECIPE)
    for k in ("epochs", "num_walks", "dec_epochs"):
        if getattr(a, k) is not None:
            cfg[k] = getattr(a, k)
    if a.config_only:
        return write_config(a.out, load_graph(a.graph, a.max_nodes, 0)[0], cfg)
    from safetensors.torch import save_file
    for s in a.seeds:
        ids, indptr, nbr, w, edges = load_graph(a.graph, a.max_nodes, s)
        rng = np.random.default_rng(s)
        W = walks(indptr, nbr, w, cfg["num_walks"], cfg["walk_length"], cfg["p"], cfg["q"], rng)
        X = finalize(skipgram(W, len(ids), cfg, s), np.diff(indptr) > 0)
        net, vl = train_decoder(X, edges, cfg, s)
        D = os.path.join(a.out, "embeddings")
        os.makedirs(D, exist_ok=True)
        os.makedirs(os.path.join(a.out, "decoders"), exist_ok=True)
        np.savez(os.path.join(D, f"node2vec_s{s}.npz"), ids=np.array(ids), X=X)
        dec_cfg = dict(edge_op="hadamard", hidden=cfg["dec_hidden"], in_dim=cfg["dim"], layers=1,
                       dropout=cfg["dec_dropout"], lr=cfg["dec_lr"], wd=cfg["dec_wd"], seed=s)
        sd = {k: t.contiguous() for k, t in net.state_dict().items()}
        save_file(sd, os.path.join(a.out, "decoders", f"decoder_s{s}.safetensors"),
                  metadata={"config": json.dumps(dec_cfg)})
        info = dict(seed=s, nodes=len(ids), edges=int(len(edges)), walks=len(W), val_loss=vl, recipe=cfg)
        json.dump(info, open(os.path.join(a.out, f"retrain_s{s}.json"), "w"), indent=1)
        print(json.dumps({k: info[k] for k in ("seed", "nodes", "edges", "walks", "val_loss")}))
        write_config(a.out, ids, cfg)


def write_config(out, ids, cfg):
    """model_config.json + nodes.txt so disatlas.model.load_ensemble(out) reads the retrained seeds."""
    from disatlas.model import _sha256
    seeds = sorted(int(f[9:-12]) for f in os.listdir(os.path.join(out, "decoders")) if f.endswith(".safetensors"))
    files = [f"embeddings/node2vec_s{s}.npz" for s in seeds] + [f"decoders/decoder_s{s}.safetensors" for s in seeds]
    open(os.path.join(out, "nodes.txt"), "w").write("\n".join(ids) + "\n")
    json.dump(dict(model="DisAtlas ensemble retrained with disatlas.retrain", seeds=seeds, recipe=cfg,
                   files={f: _sha256(os.path.join(out, f)) for f in files}),
              open(os.path.join(out, "model_config.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
