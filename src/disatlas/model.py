"""DisAtlas 5-seed ensemble: load the disease embeddings and score disease pairs by cosine.

A DisAtlas pair score is the cosine of the two 32-d disease embeddings, averaged over the five seeds
(4707-4711); the clustering distance is one minus this score. There is no decoder.

The default model is the final (gene-anchored) model of the manuscript:
  warm  data/zenodo/model_outputs/embeddings/gene_anchored/embedding_s<seed>.npz
  cold  build/model/disatlas_v3.2_ensemble/embeddings/node2vec_s<seed>.npz (python -m disatlas.retrain)
Any directory with either layout can be passed (``ids``, ``X`` float32 per seed); the earlier released v3.2
embeddings are data/zenodo/model/disatlas_v3.2_ensemble (same layout as a retrain).

    from disatlas.model import load_ensemble
    m = load_ensemble()
    m.score("MONDO:0016244", "MONDO:0100244")    # aHUS - PNH
    m.score_pairs(A, B, per_seed=True)           # (5, n) per-seed cosines
    m.top_partners("MONDO:0005147", k=10)        # highest-cosine diseases for one disease
    m.matrix()                                   # n x n float32 mean cosine (ids order = m.ids)
"""
from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import pandas as pd

from disatlas import paths

SEEDS = (4707, 4708, 4709, 4710, 4711)


def final_model_dir():
    """Final (gene-anchored) model: shipped embeddings (warm) or the cold retrain in build/model."""
    if paths.MODE == "cold":
        return paths.BUILD / "model" / paths.MODEL_NAME
    return paths.ZENODO / "model_outputs" / "embeddings" / "gene_anchored"


DEFAULT_DIR = str(final_model_dir())


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _emb_file(model_dir, s):
    for f in (os.path.join(model_dir, "embeddings", f"node2vec_s{s}.npz"), os.path.join(model_dir, f"embedding_s{s}.npz")):
        if os.path.exists(f):
            return f
    raise FileNotFoundError(f"no embedding for seed {s} in {model_dir}")


def _expected_sums(model_dir):
    """{path: sha256} from model_config.json (model bundle, retrain) or model_outputs/SHA256SUMS (shipped)."""
    cfg = os.path.join(model_dir, "model_config.json")
    if os.path.exists(cfg):
        return {os.path.join(model_dir, k): v for k, v in json.load(open(cfg)).get("files", {}).items()}
    mo = os.path.dirname(os.path.dirname(os.path.abspath(model_dir)))
    sums = os.path.join(mo, "SHA256SUMS")
    if not os.path.exists(sums):
        return {}
    out = {}
    for line in open(sums):
        h, f = line.split()
        p = os.path.join(mo, f)
        if os.path.dirname(os.path.abspath(p)) == os.path.abspath(model_dir):
            out[p] = h
    return out


class Ensemble:
    def __init__(self, model_dir: str | None = None, verify: bool = True, seeds=SEEDS):
        model_dir = str(model_dir or DEFAULT_DIR)
        self.dir = model_dir
        files = {s: _emb_file(model_dir, s) for s in seeds}
        if verify:
            exp = _expected_sums(model_dir)
            bad = [f for f in files.values() if f in exp and _sha256(f) != exp[f]]
            if bad:
                raise ValueError(f"checksum mismatch: {bad}")
        self.seeds = [int(s) for s in seeds]
        self.emb, self.idx = {}, {}
        order = None
        for s in self.seeds:
            z = np.load(files[s], allow_pickle=True)
            ids = [str(i) for i in z["ids"]]
            X = z["X"].astype(np.float64)
            self.emb[s] = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
            self.idx[s] = {v: i for i, v in enumerate(ids)}
            order = ids if order is None else order
        shared = set.intersection(*[set(d) for d in self.idx.values()])
        self._ids = [d for d in order if d in shared]

    @property
    def ids(self) -> list[str]:
        """Diseases embedded by every seed, in the order of the first seed's file."""
        return list(self._ids)

    def _cos(self, s, A, B):
        idx, X = self.idx[s], self.emb[s]
        ia = np.array([idx.get(str(a), -1) for a in A], dtype=np.int64)
        ib = np.array([idx.get(str(b), -1) for b in B], dtype=np.int64)
        ok = (ia >= 0) & (ib >= 0)
        out = np.full(len(ia), np.nan)
        out[ok] = (X[ia[ok]] * X[ib[ok]]).sum(1)
        return out

    def score_pairs(self, A, B, per_seed: bool = False) -> np.ndarray:
        """Cosine score of each pair (A[i], B[i]): per seed (n_seeds, n) or the 5-seed mean (n,); NaN if not embedded."""
        A, B = list(A), list(B)
        c = np.stack([self._cos(s, A, B) for s in self.seeds])
        return c if per_seed else c.mean(0)

    def score(self, a: str, b: str) -> float:
        return float(self.score_pairs([a], [b])[0])

    def matrix(self, ids=None) -> np.ndarray:
        """Mean cosine matrix (float32, len(ids) x len(ids)); per-seed products in float64, mean of the 5 seeds."""
        ids = self.ids if ids is None else list(ids)
        acc = None
        for s in self.seeds:
            X = self.emb[s][[self.idx[s][d] for d in ids]]
            M = X @ X.T
            acc = M if acc is None else acc + M
        return (acc / len(self.seeds)).astype(np.float32)

    def top_partners(self, disease: str, k: int = 10, exclude=()) -> pd.DataFrame:
        others = [d for d in self.ids if d != disease and d not in set(exclude)]
        p = self.score_pairs([disease] * len(others), others)
        o = np.argsort(-p, kind="stable")[:k]
        return pd.DataFrame({"disease": [others[i] for i in o], "score": p[o]})


def load_ensemble(model_dir: str | None = None, verify: bool = True) -> Ensemble:
    return Ensemble(model_dir, verify=verify)
