#!/usr/bin/env python
"""Pre-registered parity test of the P_full training graph (released specific edges + Habcd broad graph) against the
released model (cosine readout, 5-seed ensemble).

Rule, fixed before any arm was scored:
  endpoint  AUROC on the held-out IxIDN-Matched test pairs (pool == 'test', 1,094 pairs);
  delta     AUROC(P_full) - AUROC(released model) on the pairs both score;
  CI        disease-cluster bootstrap, clusters = endpoint a, B = 2,000, percentile 95 % interval, default_rng(0);
  PARITY    CI lower bound > -0.02 AND coverage >= 95 % of the test pairs;
  NOT YET   otherwise, if delta >= -0.05 and coverage >= 85 %; else FAR.
The released model is the final (gene-anchored) model of the manuscript: its seed 4707-4711 disease embeddings ship in
data/zenodo/model_outputs/embeddings/gene_anchored.

    python scripts/pfull_parity.py                                  # replay the authors' P_full run (shipped cosines)
    python scripts/pfull_parity.py --model build/model/disatlas_v3.2_ensemble   # your cold retrain on P_full
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score as AU

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
Z = os.path.join(ROOT, "data", "zenodo")
PAIRS = os.path.join(Z, "model_outputs", "pfull_reference", "pfull_ixidn_cosine_5seed.parquet")  # pair columns + per-seed P_full cosines
REC = os.path.join(Z, "model_outputs", "pfull_reference", "pfull_ixidn_cosine_5seed.parquet")
RELEASED = os.path.join(Z, "model_outputs", "embeddings", "gene_anchored")
RECORD = os.path.join(ROOT, "tests", "pfull_parity_record.json")
SEEDS = (4707, 4708, 4709, 4710, 4711)
MARGIN, COVER, NB = -0.02, 0.95, 2000


def cosines(model_dir, P):
    out = {}
    for s in SEEDS:
        f = os.path.join(model_dir, "embeddings", f"node2vec_s{s}.npz")  # a retrain (disatlas.retrain --out)
        if not os.path.exists(f):
            f = os.path.join(model_dir, f"embedding_s{s}.npz")  # shipped final-model embeddings
        z = np.load(f, allow_pickle=True)
        ix = {str(d): i for i, d in enumerate(z["ids"])}
        X = z["X"].astype(np.float64); X /= np.linalg.norm(X, axis=1, keepdims=True)
        ia = np.array([ix.get(x, -1) for x in P.a]); ib = np.array([ix.get(x, -1) for x in P.b]); ok = (ia >= 0) & (ib >= 0)
        c = np.full(len(P), np.nan); c[ok] = (X[ia[ok]] * X[ib[ok]]).sum(1)
        out[s] = c
    return out


def fauc(y, x):
    r = rankdata(x); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def boot(rng, df, x, y, nb=NB):
    g = df.groupby("a").indices; keys = list(g); L = df.label.values; X = df[x].values; Y = df[y].values; out = []
    for _ in range(nb):
        idx = np.concatenate([g[keys[i]] for i in rng.choice(len(keys), len(keys))])
        l = L[idx]
        if l.min() == l.max():
            continue
        out.append(fauc(l, X[idx]) - fauc(l, Y[idx]))
    return [float(v) for v in np.percentile(out, [2.5, 97.5])]


def parity(D):
    """D: blind_id, a, label, arm_s<seed>, ref_s<seed> for the test pairs."""
    D = D.copy()
    for k in ("arm", "ref"):
        D[k] = D[[f"{k}_s{s}" for s in SEEDS]].mean(axis=1, skipna=False)
    sh = D[D.arm.notna() & D.ref.notna()].reset_index(drop=True)
    r = dict(n=len(D), n_scored=int(D.arm.notna().sum()), coverage=float(D.arm.notna().mean()), shared=len(sh),
             auroc_arm=float(AU(sh.label, sh.arm)), auroc_ref=float(AU(sh.label, sh.ref)))
    r["dAUROC"] = r["auroc_arm"] - r["auroc_ref"]
    r["ci_lo"], r["ci_hi"] = boot(np.random.default_rng(0), sh, "arm", "ref")
    d = [AU(sh.label, sh[f"arm_s{s}"]) - AU(sh.label, sh[f"ref_s{s}"]) for s in SEEDS]
    r["seed_dAUROC_mean"], r["seed_dAUROC_sd"] = float(np.mean(d)), float(np.std(d, ddof=1))
    r["verdict"] = "PARITY" if (r["ci_lo"] > MARGIN and r["coverage"] >= COVER) else (
        "NOT YET" if (r["dAUROC"] >= -0.05 and r["coverage"] >= 0.85) else "FAR")
    return r


def table(model=None, released=RELEASED):
    P = pd.read_parquet(PAIRS, columns=["blind_id", "a", "b", "label", "pool"])
    D = P[["blind_id", "a", "label", "pool"]].copy()
    for s, c in cosines(released, P).items():
        D[f"ref_s{s}"] = c
    if model is None:
        R = pd.read_parquet(REC).set_index("blind_id")
        for s in SEEDS:
            D[f"arm_s{s}"] = D.blind_id.map(R[f"Pfull_cos_s{s}"]).values
    else:
        for s, c in cosines(model, P).items():
            D[f"arm_s{s}"] = c
    return D[D.pool == "test"].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default=None, help="model folder with embeddings/node2vec_s<seed>.npz (default: replay)")
    ap.add_argument("--out", default=None, help="write the result (json)")
    a = ap.parse_args()
    r = parity(table(a.model))
    if a.model is None:
        rec = json.load(open(RECORD))
        diff = max(abs(r[k] - rec[k]) for k in rec if isinstance(rec[k], float))
        print(f"max |replay - recorded|: {diff:.2e}")
        if diff > 1e-12 or r["verdict"] != rec["verdict"]:
            print("REPLAY MISMATCH"); sys.exit(1)
    print(f"P_full test: dAUROC {r['dAUROC']:+.4f} [{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}], "
          f"coverage {r['coverage']:.1%}, seed dAUROC {r['seed_dAUROC_mean']:+.4f} (sd {r['seed_dAUROC_sd']:.4f}) -> {r['verdict']}")
    if a.out:
        json.dump(r, open(a.out, "w"), indent=1)
    return r


if __name__ == "__main__":
    main()

