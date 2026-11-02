"""Frozen DisAtlas v3.2 ensemble: load and score disease pairs.

The model bundle (data/zenodo/model/disatlas_v3.2_ensemble) holds, for each of the 5 seeds
(4707-4711), a 32-d node2vec embedding (``embeddings/node2vec_s<seed>.npz``: ``ids``, ``X`` float32,
unit norm) and a Hadamard MLP decoder (``decoders/decoder_s<seed>.safetensors``). The ensemble
probability of a pair (a, b) is the mean over seeds of sigmoid(decoder_s(X_s[a] * X_s[b])); it is
symmetric and NaN when an id is not embedded.

    from disatlas.model import load_ensemble
    m = load_ensemble()                       # default: data/zenodo/model/disatlas_v3.2_ensemble
    m.score("MONDO:0005147", "MONDO:0005148")
    m.score_pairs(A, B, per_seed=True)        # (5, n)
    m.top_partners("MONDO:0005147", k=10)     # highest-scoring diseases for one disease
"""
from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import pandas as pd

from disatlas.paths import model_dir as _model_dir

DEFAULT_DIR = str(_model_dir())  # data/zenodo/model/... (warm) or build/model/... (cold)


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _build_net(cfg: dict):
    from torch import nn
    h, layers, p = int(cfg["hidden"]), int(cfg.get("layers", 1)), float(cfg.get("dropout", 0.0))
    if layers == 1 and not p:
        return nn.Sequential(nn.Linear(cfg["in_dim"], h), nn.ReLU(), nn.Linear(h, 1))
    mods = [nn.Linear(cfg["in_dim"], h), nn.ReLU()]
    for _ in range(layers - 1):
        mods += [nn.Dropout(p), nn.Linear(h, h), nn.ReLU()]
    mods += [nn.Dropout(p), nn.Linear(h, 1)]
    return nn.Sequential(*mods)


class Ensemble:
    def __init__(self, model_dir: str = DEFAULT_DIR, verify: bool = True):
        from safetensors import safe_open
        from safetensors.torch import load_file
        self.dir = model_dir
        self.config = json.load(open(os.path.join(model_dir, "model_config.json")))
        if verify:
            bad = [r for r, s in self.config["files"].items() if _sha256(os.path.join(model_dir, r)) != s]
            if bad:
                raise ValueError(f"checksum mismatch: {bad}")
        self.seeds = [int(s) for s in self.config["seeds"]]
        self.emb, self.idx, self.nets = {}, {}, {}
        for s in self.seeds:
            z = np.load(os.path.join(model_dir, "embeddings", f"node2vec_s{s}.npz"))
            ids = [str(i) for i in z["ids"]]
            self.emb[s] = z["X"].astype(np.float32)
            self.idx[s] = {v: i for i, v in enumerate(ids)}
            p = os.path.join(model_dir, "decoders", f"decoder_s{s}.safetensors")
            with safe_open(p, framework="pt") as f:
                cfg = json.loads(f.metadata()["config"])
            assert cfg["edge_op"] == "hadamard"
            net = _build_net(cfg)
            net.load_state_dict(load_file(p, device="cpu"))
            self.nets[s] = net.eval()

    @property
    def ids(self) -> list[str]:
        return sorted(set.intersection(*[set(d) for d in self.idx.values()]))

    def _logits(self, s, A, B, batch_size=65536):
        import torch
        idx, X = self.idx[s], self.emb[s]
        ia = np.array([idx.get(str(a), -1) for a in A])
        ib = np.array([idx.get(str(b), -1) for b in B])
        ok = (ia >= 0) & (ib >= 0)
        out = np.full(len(A), np.nan)
        if ok.any():
            F = X[ia[ok]] * X[ib[ok]]
            with torch.no_grad():
                z = [self.nets[s](torch.as_tensor(F[i:i + batch_size])).squeeze(1).numpy()
                     for i in range(0, len(F), batch_size)]
            out[ok] = np.concatenate(z)
        return out

    def logits(self, A, B) -> np.ndarray:
        """Per-seed logits, shape (n_seeds, n)."""
        A, B = list(A), list(B)
        return np.stack([self._logits(s, A, B) for s in self.seeds])

    def score_pairs(self, A, B, per_seed: bool = False) -> np.ndarray:
        p = 1.0 / (1.0 + np.exp(-self.logits(A, B)))
        return p if per_seed else p.mean(0)

    def score(self, a: str, b: str) -> float:
        return float(self.score_pairs([a], [b])[0])

    def top_partners(self, disease: str, k: int = 10, exclude=()) -> pd.DataFrame:
        others = [d for d in self.ids if d != disease and d not in set(exclude)]
        p = self.score_pairs([disease] * len(others), others)
        o = np.argsort(-p, kind="stable")[:k]
        return pd.DataFrame({"disease": [others[i] for i in o], "score": p[o]})


def load_ensemble(model_dir: str = DEFAULT_DIR, verify: bool = True) -> Ensemble:
    return Ensemble(model_dir, verify=verify)
