"""Ranking metrics used by the main figures (numpy only).

auroc(y, s, w) equals sklearn.metrics.roc_auc_score(y, s, sample_weight=w): area under the
trapezoidal ROC curve, i.e. the weighted Mann-Whitney statistic with ties credited 1/2.
"""
import numpy as np


def auroc(y, s, w=None):
    y = np.asarray(y).astype(bool); s = np.asarray(s, dtype=float)
    w = np.ones(len(y)) if w is None else np.asarray(w, dtype=float)
    order = np.argsort(s, kind="mergesort"); s, y, w = s[order], y[order], w[order]
    _, start = np.unique(s, return_index=True)
    wp = np.add.reduceat(np.where(y, w, 0.0), start)       # positive weight per distinct score
    wn = np.add.reduceat(np.where(y, 0.0, w), start)       # negative weight per distinct score
    below = np.concatenate([[0.0], np.cumsum(wn)[:-1]])    # negative weight strictly below
    P, N = wp.sum(), wn.sum()
    return float((wp * (below + 0.5 * wn)).sum() / (P * N))


def node_bootstrap_auroc(a, b, y, s, B=1000, seed=4707):
    """Disease-level bootstrap: resample diseases with replacement; a pair's weight is the product of
    its two diseases' multiplicities; draws in which one class vanishes are skipped."""
    a = np.asarray(a); b = np.asarray(b); y = np.asarray(y); s = np.asarray(s, dtype=float)
    # node order as in pandas.unique (first appearance in a0, b0, a1, b1, ...)
    inter = np.empty(2 * len(a), dtype=object); inter[0::2] = a; inter[1::2] = b
    _, first = np.unique(inter, return_index=True)
    first_nodes = inter[np.sort(first)]
    idx = {n: i for i, n in enumerate(first_nodes)}
    ai = np.array([idx[x] for x in a]); bi = np.array([idx[x] for x in b])
    r = np.random.default_rng(seed); n = len(first_nodes); out = []
    for _ in range(B):
        wt = np.bincount(r.integers(0, n, n), minlength=n); ww = wt[ai] * wt[bi]; m = ww > 0
        if len(np.unique(y[m])) < 2:
            continue
        out.append(auroc(y[m], s[m], ww[m]))
    return np.asarray(out)
