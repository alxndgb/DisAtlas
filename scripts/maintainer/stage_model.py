"""Maintainer-only: stage the frozen v3.2 ensemble (per-seed node2vec embeddings) into data/zenodo/model with a
path-free config, then verify that the loader's per-seed cosines reproduce the released v3.2 IxIDN-Matched cosine table
(data/zenodo/model_outputs/scores/ixidn_matched_cosine_released_v32.parquet). DisAtlas scores a pair by the cosine of
its embeddings; no decoder is staged.

usage: python scripts/maintainer/stage_model.py <frozen package dir> <data dir>
"""
import glob
import hashlib
import json
import os
import shutil
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))

pkg, data = sys.argv[1], sys.argv[2]
out = os.path.join(data, "zenodo", "model", "disatlas_v3.2_ensemble")
for sub in ("embeddings",):
    os.makedirs(os.path.join(out, sub), exist_ok=True)
    for f in sorted(glob.glob(os.path.join(pkg, sub, "*"))):
        shutil.copyfile(f, os.path.join(out, sub, os.path.basename(f)))
shutil.copyfile(os.path.join(pkg, "nodes.txt"), os.path.join(out, "nodes.txt"))

fz = json.load(open(os.path.join(pkg, "FROZEN_MECHDIST_MODEL.json")))
g = fz["graph"]
cfg = dict(product="DisAtlas", model="DisAtlas v3.2 ensemble (node2vec embeddings, cosine score)",
           version="3.2.0", seeds=fz["seeds"], ensemble="mean of the 5 per-seed cosines of the unit-norm embeddings", node2vec=fz["node2vec"],
           walker="weighted second-order random walks on the training graph, walk weight w_arm",
           training_graph=dict(file=None, shipped=False,
                               note="earlier 498,155-edge v3.2 training graph, not redistributed; the edge-set hash identifies it",
                               n_nodes=g["n_nodes"], n_edges=g["n_edges"], walk_weight="w_arm",
                               train_edgeset_sha256=g["train_edgeset_sha256"]))
_sig = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tests", "pfull_graph_signature.json")))
cfg["final_model_training_graph"] = dict(
    file="build/graph/pfull_graph_edges.parquet", build="python -m disatlas.graph.pfull_build --out build/graph",
    note="training graph of the final (gene-anchored) model: P_full, rebuilt from shipped inputs and checked against tests/pfull_graph_signature.json",
    n_nodes=_sig["nodes"], n_edges=_sig["edges"], walk_weight="w_arm", edgeset_sha256=_sig["edgeset_sha256"], w_arm_sha256=_sig["w_arm_sha256"])
cfg = {k: cfg[k] for k in list(cfg)[:list(cfg).index("training_graph") + 1] + ["final_model_training_graph"] + [k for k in list(cfg)[list(cfg).index("training_graph") + 1:] if k != "final_model_training_graph"]}
cfg["node2vec"].pop("num_threads", None)
cfg["node2vec"].pop("device", None)
cfg["node2vec"]["post"] = "embeddings re-attached for isolated nodes, then L2 unit-normalised; stored float32"


def sha(p):
    h = hashlib.sha256()
    h.update(open(p, "rb").read())
    return h.hexdigest()


cfg["files"] = {os.path.relpath(p, out): sha(p) for p in sorted(glob.glob(os.path.join(out, "*", "*")))
                + [os.path.join(out, "nodes.txt")]}
json.dump(cfg, open(os.path.join(out, "model_config.json"), "w"), indent=1)

from disatlas.model import load_ensemble  # noqa: E402

m = load_ensemble(out)
ref = pd.read_parquet(os.path.join(data, "zenodo", "model_outputs", "scores", "ixidn_matched_cosine_released_v32.parquet"))
C = m.score_pairs(ref.a, ref.b, per_seed=True)
for i, s in enumerate(m.seeds):
    dl = np.nanmax(np.abs(C[i] - ref[f"c_s{s}"].values))
    print("seed", s, "rows", len(ref), "max |cosine diff|", dl)
    assert dl < 1e-6
assert np.nanmax(np.abs(C.mean(0) - ref.c.values)) < 1e-6
print("loader reproduces the released per-seed and 5-seed mean cosines")
