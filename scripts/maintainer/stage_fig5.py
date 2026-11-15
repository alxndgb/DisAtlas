"""Maintainer-only: stage the inputs of Main Figure 5 (gene-anchored model, cosine readout) into
data/derived/fig5_clustering and data/derived/control_pairs.  Nothing is refitted: the K = 55 labels and the Isomap layout
are the stored clustering run, and every score is an entry of the 5-seed cosine matrix shipped on Zenodo
(cosine_matrix/cosine_ens5_13021.npy).  Source locations are passed as KEY=PATH and are never tracked.

usage: python scripts/maintainer/stage_fig5.py <data dir> KEY=PATH ...
keys:  xy           Isomap (cosine, 15 neighbours) layout, .npy (13,021 x 2, rows in the order of `ids`)
       labels       K = 55 labels, .npz with ids / labels
       display_map  frozen display map (mondo_id, disp_x, disp_y) used only to orient the layout (similarity transform)
       meta         disease metadata parquet (disease_id, disease_class, is_rare, label)
       pairs        IxIDN-Matched benchmark table (blind_id, a, b, y, pool, rare_common_pair; scores are ignored)
       cosine       5-seed cosine matrix (.npy) and cosine_ids (its row order, JSON list)
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.linalg import orthogonal_procrustes

data = sys.argv[1]
src = dict(a.split("=", 1) for a in sys.argv[2:])
out = os.path.join(data, "derived", "fig5_clustering")
ctl = os.path.join(data, "derived", "control_pairs")

z = np.load(src["labels"], allow_pickle=True)
ids = [str(d) for d in z["ids"]]; n = len(ids); ix = {d: i for i, d in enumerate(ids)}
assert n == 13021, n
XY = np.load(src["xy"]).astype(float); assert XY.shape == (n, 2)
cos_ids = json.load(open(src["cosine_ids"])); assert cos_ids == ids, "cosine matrix row order differs from the labels"
C = np.load(src["cosine"], mmap_mode="r"); assert C.shape == (n, n)

# similarity transform (rotation/reflection + uniform scale + shift) onto the frozen display map, for orientation only
c21 = pd.read_csv(src["display_map"], sep="\t").set_index("mondo_id")
com = [d for d in ids if d in c21.index]; A = XY[[ix[d] for d in com]]; B = c21.loc[com, ["disp_x", "disp_y"]].to_numpy(float)
Am, Bm = A.mean(0), B.mean(0); Rm, _ = orthogonal_procrustes(A - Am, B - Bm); sc = np.linalg.norm(B - Bm) / np.linalg.norm(A - Am)
XYa = (XY - Am) @ Rm * sc + Bm
dm = pd.read_parquet(src["meta"]).set_index("disease_id")
M = pd.DataFrame({"mondo_id": ids, "isomap_x": XY[:, 0], "isomap_y": XY[:, 1], "map_x": XYa[:, 0], "map_y": XYa[:, 1],
                  "cluster": np.asarray(z["labels"]).astype(int)}).set_index("mondo_id")
M = M.join(dm[["disease_class", "is_rare", "label"]], how="left").reset_index()
assert M.cluster.nunique() == 55
os.makedirs(out, exist_ok=True)
M.to_csv(os.path.join(out, "disease_map_k55.tsv.gz"), sep="\t", index=False,
         compression={"method": "gzip", "mtime": 0})
print("disease_map_k55.tsv.gz", M.shape, "procrustes scale", sc)

S = pd.read_csv(src["pairs"], sep="\t")[["blind_id", "a", "b", "y", "pool", "rare_common_pair"]]
assert len(S) == 4364 and S.a.isin(ix).all() and S.b.isin(ix).all()
S["p"] = [float(C[ix[a], ix[b]]) for a, b in zip(S.a, S.b)]
S.to_csv(os.path.join(out, "ixidn_matched_pair_scores.tsv"), sep="\t", index=False)
print("ixidn_matched_pair_scores.tsv", S.shape)

PAIRS = [("MONDO:0001734", "MONDO:0011705"), ("MONDO:0016168", "MONDO:0018088"), ("MONDO:0016244", "MONDO:0100244"),
         ("MONDO:0005011", "MONDO:0005101"), ("MONDO:0004782", "MONDO:0005148"), ("MONDO:0018208", "MONDO:0007039")]
P = pd.DataFrame([dict(a=a, b=b, p=float(C[ix[a], ix[b]])) for a, b in PAIRS])
P.to_csv(os.path.join(ctl, "fig5_panel_a_pairs_p.tsv"), sep="\t", index=False)
print("fig5_panel_a_pairs_p.tsv", P.shape)
