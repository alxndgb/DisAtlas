"""Shared loaders and metrics for the comparator benchmarks (Results, cross-system comparison; Figure 6).

Reads only repo-relative stored outputs. No network access, no model calls.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
BENCH = os.path.join(REPO_ROOT, "benchmarks")
KEY = os.path.join(REPO_ROOT, "data", "08S_eval", "ixidn4364_scores.tsv")


def load_key():
    """IxIDN-Matched answer key (4,364 pairs, 50% positive) with the full-DisAtlas scores and pair covariates."""
    k = pd.read_csv(KEY, sep="\t")
    assert len(k) == 4364 and k.blind_id.is_unique and int(k.y.sum()) == 2182
    return k


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def ranked_last(v):
    """Unscored pairs are ranked last (Figure 6 convention)."""
    v = np.asarray(v, float)
    return np.where(np.isnan(v), np.nanmin(v) - 1e-9, v)


def auc(y, s):
    y = np.asarray(y); r = rankdata(s); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def exp_prec_at_k(y, s, k):
    """Expected precision of the top-k pairs; ties at the cut-off shared fractionally."""
    o = np.sort(s)[::-1]; t = o[k - 1]; gt = s > t; eq = s == t
    need = k - gt.sum()
    return (y[gt].sum() + (y[eq].mean() * need if eq.any() else 0)) / k


def ap(y, s):
    d = pd.DataFrame({"s": s, "y": y}).groupby("s", sort=True)["y"].agg(["sum", "count"]).iloc[::-1]
    tp = d["sum"].cumsum().values; n = d["count"].cumsum().values
    return float(np.sum((tp / n) * (d["sum"].values / y.sum())))


METRICS = ["AUROC", "AUPRC", "F1@k=#positives", "Precision@top10%", "NPV@bottom10%"]


def ranking_metrics(y, s):
    """Figure 6b metrics. F1 at k = number of positives equals R-precision."""
    kp = int(y.sum()); kk = int(round(0.1 * len(y)))
    return {"AUROC": auc(y, s), "AUPRC": ap(y, s), "F1@k=#positives": exp_prec_at_k(y, s, kp),
            "Precision@top10%": exp_prec_at_k(y, s, kk), "NPV@bottom10%": exp_prec_at_k(1 - y, -s, kk)}


def pair_bootstrap_metrics(y, scores, n_boot=500, seed=4707):
    """Pair bootstrap (500 draws, same draws for every system); returns long table method/metric/value/lo/hi."""
    rng = np.random.default_rng(seed)
    bi = [rng.integers(0, len(y), len(y)) for _ in range(n_boot)]
    rows = []
    for name, s in scores.items():
        m = ranking_metrics(y, s); bs = pd.DataFrame([ranking_metrics(y[i], s[i]) for i in bi])
        for k, v in m.items():
            rows.append(dict(method=name, metric=k, value=v, lo=bs[k].quantile(.025), hi=bs[k].quantile(.975)))
    return pd.DataFrame(rows)


def topk_weights(v, k):
    x = np.where(np.isnan(v), -np.inf, v); w = np.zeros(len(x)); c = np.sort(x)[::-1][k - 1]
    above = x > c; tie = x == c; w[above] = 1; w[tie] = (k - above.sum()) / tie.sum()
    return w


def topk_jaccard(raw_scores, frac=0.1):
    """Overlap of top-10% pair sets (fractional ties), Figure 6d."""
    names = list(raw_scores); n = len(next(iter(raw_scores.values()))); k = int(round(frac * n))
    W = {a: topk_weights(raw_scores[a], k) for a in names}
    J = [[np.minimum(W[a], W[b]).sum() / np.maximum(W[a], W[b]).sum() for b in names] for a in names]
    return pd.DataFrame(J, index=names, columns=names)


def run_agreement(df, cols, method):
    """Pairwise Spearman rank correlation between repeat runs (or training seeds), Figure 6c."""
    out = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            a, b = df[cols[i]].values, df[cols[j]].values; ok = ~np.isnan(a) & ~np.isnan(b)
            out.append(dict(method=method, run_a=cols[i], run_b=cols[j], spearman=spearmanr(a[ok], b[ok])[0]))
    return pd.DataFrame(out)


def disease_cluster_bootstrap(df, cols, n_boot=2000, seed=4707):
    """Disease-clustered bootstrap with one Poisson(1) weight per disease (pair weight = w_a * w_b).

    Returns point AUROCs and the bootstrap AUROC draws (same resamples for every column, so paired)."""
    rng = np.random.default_rng(seed)
    dis = np.array(sorted(set(df.a) | set(df.b))); di = {d: i for i, d in enumerate(dis)}
    ia, ib = df.a.map(di).values, df.b.map(di).values; y = df.y.values
    V = {c: df[c].values for c in cols}
    pt = {c: auc(y, V[c]) for c in cols}; bs = {c: [] for c in cols}
    for _ in range(n_boot):
        w = rng.poisson(1.0, len(dis)); idx = np.repeat(np.arange(len(df)), w[ia] * w[ib])
        if len(set(y[idx])) < 2:
            continue
        for c in cols:
            bs[c].append(auc(y[idx], V[c][idx]))
    return pt, {c: np.array(v) for c, v in bs.items()}


def write_tsv(df, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, sep="\t", index=False, float_format="%.6f", lineterminator="\n")
