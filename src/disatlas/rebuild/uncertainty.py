"""Rebuild data/derived/uncertainty/{seed_vs_sample_sd.tsv, headline_cosine_ci.tsv} (Supplementary Figure
SI_uncertainty a/b and the seed-SD against bootstrap-SD sentence of the SI uncertainty section).

    python -m disatlas.rebuild.uncertainty OUTDIR

Input: data/derived/fig6_comparison/per_run_scores.tsv (DisAtlas = five-seed mean cosine; DisAtlas|s47xx per seed).
Scopes: all 4,364 IxIDN-Matched pairs and the 1,094 test-split pairs (pool == "test").
Seed SD: SD (ddof 1) of AUROC over the five per-seed cosines.
Seed-vs-sample stream (authors' uncertainty_seed_vs_sample.py): diseases indexed in order of appearance in
concat(a, b); one default_rng(4707) stream; per replicate a disease-cluster draw (pair weight w_a * w_b), then a seed
index, then a pair draw. sample_sd_cluster = SD (ddof 1) of the ensemble AUROC over the B = 2,000 cluster draws;
AUPRC interval = 2.5/97.5 percentiles of the ensemble AUPRC over the same draws.
AUROC interval (authors' headline_auroc_g1 table, rounded to 6 decimals): diseases in order of appearance in the
row-wise (a, b) sequence, default_rng(4707), B = 2,000 cluster draws, draws without both classes skipped.
Test: tests/test_rebuild_uncertainty.py.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = Path(__file__).resolve().parents[3]
SEEDS = [f"DisAtlas|s{s}" for s in range(4707, 4712)]
B = 2000


def scopes():
    S = pd.read_csv(ROOT / "data/derived/fig6_comparison/per_run_scores.tsv", sep="\t"); assert len(S) == 4364
    return [("all", S.reset_index(drop=True)), ("test", S[S.pool == "test"].reset_index(drop=True))]


def seed_vs_sample(D):
    y = D.y.values.astype(int); n = len(D); s = D.DisAtlas.values
    dis = pd.Index(pd.unique(pd.concat([D.a, D.b]))); ia, ib = dis.get_indexer(D.a), dis.get_indexer(D.b)
    rng = np.random.default_rng(4707); au, ap = [], []
    for _ in range(B):
        c = np.bincount(rng.integers(0, len(dis), len(dis)), minlength=len(dis)); w = (c[ia] * c[ib]).astype(float); m = w > 0
        au.append(roc_auc_score(y[m], s[m], sample_weight=w[m])); ap.append(average_precision_score(y[m], s[m], sample_weight=w[m]))
        rng.integers(0, 5); rng.integers(0, n, n)          # seed pick and pair draw of the combined / pair bootstraps
    seed_sd = float(np.std([roc_auc_score(y, D[c].values) for c in SEEDS], ddof=1))
    return seed_sd, float(np.std(au, ddof=1)), np.percentile(ap, [2.5, 97.5])


def auroc_ci(D):
    y = D.y.values.astype(int); s = D.DisAtlas.values
    nodes = pd.unique(D[["a", "b"]].values.ravel()); ix = {d: i for i, d in enumerate(nodes)}
    ai, bi = D.a.map(ix).values, D.b.map(ix).values; rng = np.random.default_rng(4707); out = []
    for _ in range(B):
        w = np.bincount(rng.integers(0, len(nodes), len(nodes)), minlength=len(nodes)); ww = (w[ai] * w[bi]).astype(float); m = ww > 0
        if len(np.unique(y[m])) < 2: continue
        out.append(roc_auc_score(y[m], s[m], sample_weight=ww[m]))
    return np.percentile(out, [2.5, 97.5])


def compute():
    sv, hc = [], []
    for sc, D in scopes():
        st = {"all": "all_4364", "test": "test_test_heldout"}[sc]; y = D.y.values.astype(int)
        seed_sd, cl_sd, (plo, phi) = seed_vs_sample(D); lo, hi = auroc_ci(D)
        sv.append(dict(stratum=st, metric="auroc", seed_sd=seed_sd, sample_sd_cluster=cl_sd))
        hc.append(["AUROC p", sc, len(D), int(y.sum()), round(roc_auc_score(y, D.DisAtlas), 6), round(lo, 6), round(hi, 6), B])
        hc.append(["AUPRC p", sc, len(D), int(y.sum()), average_precision_score(y, D.DisAtlas), plo, phi, B])
    return pd.DataFrame(sv), pd.DataFrame(hc, columns=["metric", "stratum", "n", "n_pos", "value", "lo", "hi", "B"])


if __name__ == "__main__":
    o = Path(sys.argv[1]); o.mkdir(parents=True, exist_ok=True); sv, hc = compute()
    sv.to_csv(o / "seed_vs_sample_sd.tsv", sep="\t", index=False); hc.to_csv(o / "headline_cosine_ci.tsv", sep="\t", index=False)
    print("wrote", o)
