"""Compositionality on the released DisAtlas v3.2 embeddings (five seeds 4707-4711). No training.
(i)  the two worked examples (MCTD; thymoma): barycentre rank with the 1,000-draw random-component null (rng 7),
     saving the null draws;
(ii) systematic test over MONDO diseases with >= 2 embedded is_a parents (rng 20261011).
Outputs: data/derived/si_compositionality/ (the precomputed tables read by si_compositionality.py).
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import json, os, time, hashlib
from collections import defaultdict
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
from disatlas.sibio_common import model_dir; import os as _os
REL = str(model_dir())
OBO = _os.environ.get("MONDO_OBO")  # optional full ontology; default: data/derived/si_compositionality/mondo_nodes.tsv
OUT = str(DERIVED / "si_compositionality"); os.makedirs(OUT, exist_ok=True)
SEEDS = [4707, 4708, 4709, 4710, 4711]
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
nodes = open(f"{REL}/nodes.txt").read().split(); assert len(nodes) == 13021
XS = []
for s in SEEDS:
    z = np.load(f"{REL}/embeddings/node2vec_s{s}.npz", allow_pickle=True); assert list(z["ids"]) == nodes
    X = z["X"].astype(np.float64); XS.append(X / np.linalg.norm(X, axis=1, keepdims=True))
N = len(nodes); idx = {d: i for i, d in enumerate(nodes)}
# ---------------- MONDO (ontology file if MONDO_OBO is set, else the shipped per-node slice; identical results)
from disatlas.mondo_nodes import load_terms  # noqa: E402
terms, mondo_of, obo_md5 = load_terms(nodes, OBO)
name_of = lambda d: terms[mondo_of[d]]["name"] if mondo_of.get(d) and "name" in terms[mondo_of[d]] else d
# ---------------- (i) worked examples, exact replica of g2_embed_cluster_final.py bary_rank + null
def bary_rank(target, comps):
    S_ = [idx[c] for c in comps]
    sim = np.mean([X @ (X[S_].mean(0) / np.linalg.norm(X[S_].mean(0))) for X in XS], 0)
    keep = np.ones(N, bool); keep[S_] = False
    return int((sim[keep] > sim[idx[target]]).sum()) + 1, int(keep.sum()) - 1
EX = [("Mixed connective tissue disease", "MONDO:0005854", ["MONDO:0007915", "MONDO:0005100", "MONDO:0010030", "MONDO:0019434", "MONDO:0019127"]),
      ("Thymoma", "MONDO:0006456", ["MONDO:0009688", "MONDO:0016367"])]
exr, nulls = [], []
for lab, tgt, comps in EX:
    r, nother = bary_rank(tgt, comps)
    rg = np.random.default_rng(7); pool = [d for d in nodes if d != tgt]
    null = np.array([bary_rank(tgt, list(rg.choice(pool, len(comps), replace=False)))[0] for _ in range(1000)])
    p = (1 + (null <= r).sum()) / 1001
    exr.append(dict(example=lab, target_id=tgt, target=name_of(tgt), components="; ".join(f"{name_of(c)} ({c})" for c in comps),
                    n_components=len(comps), rank=r, out_of=nother, pct=r / nother, null_median_rank=float(np.median(null)), p=p))
    nulls += [dict(example=lab, draw=i, null_rank=int(v), null_pct=int(v) / nother) for i, v in enumerate(null)]
    log("example", exr[-1])
pd.DataFrame(exr).to_csv(f"{OUT}/examples.tsv", sep="\t", index=False)
pd.DataFrame(nulls).to_csv(f"{OUT}/examples_null_draws.tsv.gz", sep="\t", index=False,
                          compression={"method": "gzip", "mtime": 0})  # no timestamp: byte-reproducible
# ---------------- (ii) systematic test (RNG seed fixed as in the original analysis)
RNG = np.random.default_rng(20261011)
def cos_all(v): return np.mean([X @ (vv / np.linalg.norm(vv)) for X, vv in zip(XS, v)], axis=0)
def bary(S): return [X[S].mean(0) for X in XS]
sys_rows = []
for d in nodes:
    m = mondo_of[d]
    if not m: continue
    par = [p for p in terms[m]["is_a"] if p in idx and p != d]
    if len(par) < 2: continue
    sys_rows.append((idx[d], [idx[p] for p in par]))
log("systematic targets", len(sys_rows))
pool_all = np.arange(N); res = []
for ti, Si in sys_rows:
    n_c = N - len(Si)
    s = cos_all(bary(Si)); excl = np.zeros(N, bool); excl[Si] = True
    r = 1 + np.sum(s[~excl] > s[ti])
    cpar = cos_all([X[ti] for X in XS])[Si]
    Rn = RNG.choice(np.setdiff1d(pool_all, [ti]), len(Si), replace=False)
    sn = cos_all(bary(list(Rn))); exn = np.zeros(N, bool); exn[Rn] = True
    rn = 1 + np.sum(sn[~exn] > sn[ti])
    keep = Si[int(np.argmax(cpar))]
    Rh = [keep] + list(RNG.choice(np.setdiff1d(pool_all, [ti] + Si), len(Si) - 1, replace=False))
    sh = cos_all(bary(Rh)); exh = np.zeros(N, bool); exh[Rh] = True
    rh = 1 + np.sum(sh[~exh] > sh[ti])
    sp = cos_all([X[keep] for X in XS])
    rs = 1 + np.sum(sp[~excl] > sp[ti])
    res.append(dict(target=nodes[ti], target_name=name_of(nodes[ti]), n_parents=len(Si), parents="; ".join(name_of(nodes[x]) for x in Si),
                    n_candidates=n_c, rank_bary=int(r), pct_bary=r / n_c, rank_best_single_parent=int(rs), pct_best_single_parent=rs / n_c,
                    rank_null_one_real_rest_random=int(rh), pct_null_one_real_rest_random=rh / n_c,
                    rank_null_random_parents=int(rn), pct_null_random_parents=rn / n_c))
R = pd.DataFrame(res); R.round(6).to_csv(f"{OUT}/systematic_mondo_multiparent.tsv", sep="\t", index=False)
def wilson(k, n, z=1.959964):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d; return c - h, c + h
rows = []
for arm, col in [("Barycentre of real parents", "bary"), ("Best single real parent", "best_single_parent"),
                 ("One real parent + random others", "null_one_real_rest_random"), ("Random parents", "null_random_parents")]:
    k = int((R[f"pct_{col}"] <= 0.01).sum()); lo, hi = wilson(k, len(R))
    rows.append(dict(query=arm, n=len(R), median_pct=float(R[f"pct_{col}"].median()), top1pct_k=k, top1pct=k / len(R), top1pct_lo=lo, top1pct_hi=hi,
                     top10=float((R[f"rank_{col}"] <= 10).mean())))
S = pd.DataFrame(rows); S.to_csv(f"{OUT}/systematic_summary.tsv", sep="\t", index=False)
summ = dict(n_targets=len(R), parents_median=float(R.n_parents.median()),
            share_bary_better_than_best_single=float((R.rank_bary < R.rank_best_single_parent).mean()),
            share_bary_ties_best_single=float((R.rank_bary == R.rank_best_single_parent).mean()),
            share_bary_worse_than_best_single=float((R.rank_bary > R.rank_best_single_parent).mean()),
            wilcoxon_bary_vs_one_real_rest_random_p=float(wilcoxon(R.pct_bary, R.pct_null_one_real_rest_random).pvalue),
            wilcoxon_bary_vs_random_p=float(wilcoxon(R.pct_bary, R.pct_null_random_parents).pvalue),
            wilcoxon_bary_vs_best_single_p=float(wilcoxon(R.pct_bary, R.pct_best_single_parent).pvalue),
            inputs={"nodes.txt": md5(f"{REL}/nodes.txt"), "mondo.obo": md5(OBO) if OBO else obo_md5,
                    **{f"node2vec_s{s}.npz": md5(f"{REL}/embeddings/node2vec_s{s}.npz") for s in SEEDS}},
            seeds=SEEDS, rng_examples=7, rng_systematic=20261011, n_null_examples=1000)
json.dump(summ, open(f"{OUT}/compositionality_summary.json", "w"), indent=1)
log("summary", {k: v for k, v in summ.items() if k != "inputs"}); log(S.to_string()); log("done")
