"""Train the UPNA-DDAL MLP decoder on top of a per-seed Node2Vec embedding.

Reproduces upna_ddal.{sample_training_pairs, train_val_split, make_decoder, train_decoder} from
the source repo. Positives = G_train edges; negatives = sampled non-edges (5x, uniform scheme).
Early stopping on validation AUROC (disease-disjoint split, val_frac=0.15), best weights restored.

Note on pos_weight: the `disclosure` field of FROZEN_MECHDIST_MODEL.json states pos_weight was changed
from 1 to 10 relative to the tuning-time selection record, while its `decoder.pos_weight` field is 1.0
and the per-seed build records (build/decoder_s<seed>.json) do not record it. This script defaults to
the recorded value (1.0); pass --pos-weight 10 for the disclosed setting.
"""
import argparse
import json
import os
import sys

import numpy as np
import torch
import torch.nn as nn

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sort_order import original_order  # noqa: E402

HP = dict(edge_op="hadamard", hidden=128, layers=1, dropout=0.3, lr=0.001, wd=0.0001,
          neg_ratio=5.0, batch_size=1024, patience=8, epochs=100, val_frac=0.15)


class Decoder:
    def __init__(self, in_dim, hidden=128, layers=1, dropout=0.3, lr=0.001, wd=0.0001, seed=4707, device="cpu"):
        torch.manual_seed(seed)
        mods = []
        d = in_dim
        for _ in range(layers):
            mods += [nn.Linear(d, hidden), nn.ReLU(), nn.Dropout(dropout)]
            d = hidden
        mods.append(nn.Linear(d, 1))
        self.net = nn.Sequential(*mods).to(device)
        self.opt = torch.optim.Adam(self.net.parameters(), lr=lr, weight_decay=wd)
        self.device = device


def edge_features(emb, pairs, op="hadamard"):
    a = np.stack([emb[u] for u, v in pairs]).astype(np.float32)
    b = np.stack([emb[v] for u, v in pairs]).astype(np.float32)
    if op == "hadamard":
        return a * b
    raise ValueError(op)


def sample_training_pairs(G, nodes, neg_ratio, seed):
    rng = np.random.default_rng(seed)
    pos = [tuple(e) for e in G.edges()]
    node_arr = np.array(nodes)
    neg, seen = [], {frozenset(e) for e in pos}
    n_neg_target = int(len(pos) * neg_ratio)
    while len(neg) < n_neg_target:
        batch = rng.choice(node_arr, size=(n_neg_target - len(neg), 2))
        for u, v in batch:
            if u == v:
                continue
            fs = frozenset((u, v))
            if fs in seen:
                continue
            seen.add(fs)
            neg.append((u, v))
    return pos, neg[:n_neg_target]


def disease_disjoint_split(pos_pairs, neg_pairs, val_frac, seed):
    rng = np.random.default_rng(seed)
    nodes = original_order(list({n for p in list(pos_pairs) + list(neg_pairs) for n in p}))
    perm = list(nodes)
    rng.shuffle(perm)
    cut = int(round(len(perm) * (1 - val_frac)))
    tr_nodes, va_nodes = set(perm[:cut]), set(perm[cut:])

    def _restrict(pairs, label):
        tr, va = [], []
        for u, v in pairs:
            if u in tr_nodes and v in tr_nodes:
                tr.append((u, v, label))
            elif u in va_nodes and v in va_nodes:
                va.append((u, v, label))
        return tr, va

    pos_tr, pos_va = _restrict(pos_pairs, 1)
    neg_tr, neg_va = _restrict(neg_pairs, 0)
    train = pos_tr + neg_tr
    val = pos_va + neg_va
    rng.shuffle(train)
    return train, val


def auroc(y_true, y_score):
    order = np.argsort(-y_score)
    y = np.asarray(y_true)[order]
    n_pos, n_neg = y.sum(), len(y) - y.sum()
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    tp = np.cumsum(y)
    fp = np.cumsum(1 - y)
    tpr = tp / n_pos
    fpr = fp / n_neg
    return float(np.trapz(tpr, fpr))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--embedding", required=True, help=".npz from 02_train_node2vec.py (keys 'ids','X')")
    ap.add_argument("--graph", default="../../data/zenodo/network/G_train.pickle")
    ap.add_argument("--seed", type=int, required=True, choices=[4707, 11, 202, 3407, 8191])
    ap.add_argument("--pos-weight", type=float, default=1.0,
                    help="see module docstring CAVEAT -- recorded value is 1.0, disclosed value is 10")
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    import pickle
    with open(args.graph, "rb") as fh:
        G = pickle.load(fh)
    npz = np.load(args.embedding, allow_pickle=True)
    emb = {k: v for k, v in zip(npz["ids"], npz["X"])}

    nodes = [n for n in G.nodes() if n in emb]
    pos_pairs, neg_pairs = sample_training_pairs(G, nodes, HP["neg_ratio"], args.seed)
    train, val = disease_disjoint_split(pos_pairs, neg_pairs, HP["val_frac"], args.seed)
    print(f"n_train={len(train)} n_val={len(val)}")

    dim = next(iter(emb.values())).shape[0]
    dec = Decoder(dim, hidden=HP["hidden"], layers=HP["layers"], dropout=HP["dropout"],
                 lr=HP["lr"], wd=HP["wd"], seed=args.seed)

    def _batch_xy(rows):
        pairs = [(u, v) for u, v, _ in rows]
        y = np.array([lab for _, _, lab in rows], dtype=np.float32)
        X = edge_features(emb, pairs, HP["edge_op"])
        return torch.as_tensor(X), torch.as_tensor(y)

    Xval, yval = _batch_xy(val)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(args.pos_weight))
    best_auroc, best_state, wait = -1.0, None, 0
    rng = np.random.default_rng(args.seed)
    for epoch in range(HP["epochs"]):
        order = rng.permutation(len(train))
        for i in range(0, len(order), HP["batch_size"]):
            idx = order[i:i + HP["batch_size"]]
            rows = [train[j] for j in idx]
            Xb, yb = _batch_xy(rows)
            logit = dec.net(Xb).squeeze(-1)
            loss = loss_fn(logit, yb)
            dec.opt.zero_grad(set_to_none=True)
            loss.backward()
            dec.opt.step()
        with torch.no_grad():
            val_logit = dec.net(Xval).squeeze(-1).numpy()
        val_auc = auroc(yval.numpy(), val_logit)
        improved = val_auc > best_auroc
        if improved:
            best_auroc, best_state, wait = val_auc, {k: v.clone() for k, v in dec.net.state_dict().items()}, 0
        else:
            wait += 1
        print(f"epoch {epoch} val_auroc {val_auc:.4f} best {best_auroc:.4f}")
        if wait >= HP["patience"]:
            break
    dec.net.load_state_dict(best_state)

    out = args.output or f"../../data/zenodo/models/decoder_s{args.seed}_rebuilt.pt"
    torch.save({"state_dict": dec.net.state_dict(),
               "config": {**HP, "pos_weight": args.pos_weight, "seed": args.seed, "in_dim": dim}}, out)
    print("wrote", out, "best_val_auroc", best_auroc)


if __name__ == "__main__":
    main()
