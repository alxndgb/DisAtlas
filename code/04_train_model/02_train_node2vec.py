"""Train the per-seed Node2Vec embedding of G_train (PyTorch Geometric backend).

Reproduces embed_node2vec.node2vec_embeddings from the source repo: dim=32, walk_length=20,
window_size=2 (context_size=3), num_walks=20, p=1, q=0.5, epochs=25, lr=0.01, batch_size=128,
num_negative_samples=5, neg_alpha=1.0 (degree-proportional negative sampling), start_beta=0.5
(low-degree-biased walk starts). These non-default negative-sampling/walk-start knobs are part of
the frozen-model contract (FROZEN_MECHDIST_MODEL.json) -- changing them changes the embedding.

Run once per seed in {4707, 11, 202, 3407, 8191}. CPU, single-threaded, for bitwise reproducibility
(the tuning runs used 8 threads and are therefore only statistically, not bitwise, reproducible).
"""
import argparse
import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sort_order import original_order  # noqa: E402

import numpy as np

HP = dict(dim=32, walk_length=20, num_walks=20, window_size=2, p=1, q=0.5,
          epochs=25, lr=0.01, num_negative_samples=5, neg_alpha=1.0, start_beta=0.5,
          batch_size=128, min_count=1, weighted=False, extend=False)


def _edge_index(G_str):
    import torch
    node_ids = original_order([n for n in G_str.nodes() if G_str.degree(n) > 0])
    idx = {n: i for i, n in enumerate(node_ids)}
    src, dst = [], []
    for u, v in G_str.edges():
        if u == v:
            continue
        src += [idx[u], idx[v]]
        dst += [idx[v], idx[u]]
    return torch.tensor([src, dst], dtype=torch.long), node_ids


def node2vec_embeddings(G_str, seed, dim=32, walk_length=20, num_walks=20, window_size=2,
                        p=1, q=0.5, epochs=25, lr=0.01, num_negative_samples=5,
                        neg_alpha=1.0, start_beta=0.5, batch_size=128, device="cpu"):
    os.environ.setdefault("PYTHONHASHSEED", "0")
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    import random
    import torch
    from torch_geometric.nn import Node2Vec

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    edge_index, node_ids = _edge_index(G_str)
    context_size = min(window_size + 1, walk_length)
    model = Node2Vec(edge_index, embedding_dim=dim, walk_length=walk_length, context_size=context_size,
                     walks_per_node=num_walks, p=p, q=q, num_negative_samples=num_negative_samples,
                     num_nodes=len(node_ids), sparse=False).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    gen = torch.Generator().manual_seed(seed)
    deg = np.array([G_str.degree(v) for v in node_ids], float)
    neg_p = torch.as_tensor(deg ** neg_alpha / (deg ** neg_alpha).sum(), dtype=torch.float)
    start_p = torch.as_tensor(deg ** -start_beta / (deg ** -start_beta).sum(), dtype=torch.float)
    n_nodes = len(node_ids)

    def _neg_sample(batch):
        if neg_alpha == 0:
            return model.neg_sample(batch)
        bb = batch.repeat(num_walks * num_negative_samples)
        rw = torch.multinomial(neg_p, bb.numel() * walk_length, replacement=True,
                               generator=gen).view(-1, walk_length)
        rw = torch.cat([bb.view(-1, 1), rw], dim=-1)
        return torch.cat([rw[:, j:j + context_size] for j in range(1 + walk_length + 1 - context_size)], dim=0)

    model.train()
    for _ in range(epochs):
        perm = (torch.randperm(n_nodes, generator=gen) if start_beta == 0 else
                torch.multinomial(start_p, n_nodes, replacement=True, generator=gen))
        for i in range(0, len(perm), batch_size):
            batch = perm[i:i + batch_size]
            loss = model.loss(model.pos_sample(batch).to(device), _neg_sample(batch).to(device))
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()

    mat = model.embedding.weight.detach().cpu().numpy().astype(np.float32)
    emb_dict = {nid: mat[i] for i, nid in enumerate(node_ids)}
    return emb_dict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", default="../../data/zenodo/network/G_train.pickle")
    ap.add_argument("--seed", type=int, required=True, choices=[4707, 11, 202, 3407, 8191])
    ap.add_argument("--output-dir", default="../../data/zenodo/models")
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    import torch
    torch.set_num_threads(args.threads)

    with open(args.graph, "rb") as fh:
        G_train = pickle.load(fh)

    emb_dict = node2vec_embeddings(G_train, seed=args.seed, **HP)
    # L2 unit-norm post-processing (modeling_common.finalize_embedding), float32
    ids = original_order(list(emb_dict))
    X = np.stack([emb_dict[i].astype(np.float64) for i in ids])
    X = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    X = X.astype(np.float32)

    out_path = os.path.join(args.output_dir, f"node2vec_s{args.seed}_rebuilt.npz")
    np.savez(out_path, ids=np.array(ids), X=X)
    print("wrote", out_path, "shape", X.shape)


if __name__ == "__main__":
    main()
