"""Maintainer-only: stage the frozen v3.2 ensemble into data/zenodo/model with a path-free config,
then verify the loader reproduces the stored per-seed
logits of the frozen evaluation protocol.

usage: python scripts/maintainer/stage_model.py <frozen package dir> <data dir>
"""
import glob
import gzip
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
for sub in ("embeddings", "decoders"):
    os.makedirs(os.path.join(out, sub), exist_ok=True)
    for f in sorted(glob.glob(os.path.join(pkg, sub, "*"))):
        shutil.copyfile(f, os.path.join(out, sub, os.path.basename(f)))
shutil.copyfile(os.path.join(pkg, "nodes.txt"), os.path.join(out, "nodes.txt"))

fz = json.load(open(os.path.join(pkg, "FROZEN_MECHDIST_MODEL.json")))
g = fz["graph"]
cfg = dict(product="DisAtlas", model="DisAtlas v3.2 ensemble (node2vec + Hadamard MLP decoder)",
           version="3.2.0", seeds=fz["seeds"], ensemble=fz["ensemble"], node2vec=fz["node2vec"],
           walker="weighted second-order random walks on the training graph, walk weight w_arm",
           decoder={k: v for k, v in fz["decoder"].items() if k != "seed"},
           training_graph=dict(file="graph/union_graph_edges.parquet", n_nodes=g["n_nodes"],
                               n_edges=g["n_edges"], walk_weight="w_arm",
                               train_edgeset_sha256=g["train_edgeset_sha256"]))
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
for s in m.seeds:
    src = os.path.join(pkg, "eval", f"ixidn_frozen_protocol_scores_s{s}.tsv.gz")
    d = pd.read_csv(src, sep="\t")
    d = d.rename(columns={"a": "disease_a", "b": "disease_b"})
    lg = m._logits(s, d.disease_a.tolist(), d.disease_b.tolist())
    dl = np.nanmax(np.abs(lg - d.logit.values))
    print("seed", s, "rows", len(d), "max |logit diff|", dl)
    assert dl < 1e-5
print("loader reproduces stored per-seed logits for all seeds")
