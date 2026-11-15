"""Rebuild data/derived/degree_controls/degree_controls_table.tsv (Supplementary Table S8 (SI_degree_controls): degree-matched and
within-disease AUROC of every system on the IxIDN-Matched benchmark, with the node-disjoint DisAtlas arm as reference).

    python -m disatlas.rebuild.degree_controls OUTDIR [B]

Inputs (all released): data/derived/fig6_comparison/per_run_scores.tsv (per-pair scores of every system, benchmark
split `pool`), data/derived/fig6_comparison/node_disjoint_per_seed.tsv (held-out arm, five seeds), and the training
graph data/raw/training_graph/training_graph_edges.parquet (the graph the embedding was trained on, benchmark and FDA
pairs stripped) for the degree of each disease.
Method: system score = mean over the runs of the system (DisAtlas: five-seed cosine; held-out arm: mean of its five
seeds, pairs with a missing seed unscored); Phylo Biomni and TxGNN rank their unscored pairs last, the held-out arm drops
its unscorable pairs. ONE degree-matched subset per case set (all pairs / test split), chosen without scores: quantile-
binned (10 bins) match of log1p(deg a) + log1p(deg b) between positives and negatives, numpy default_rng(4707), sorted
pair order. Within-disease AUROC = mean over diseases of the per-disease AUROC. Intervals: one shared set of B = 2000
disease-cluster draws (938 diseases, pair weight w_a * w_b, default_rng(4707)); differences use the same draws.
Test: tests/test_rebuild_degree_controls.py.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[3]
FIG6 = ROOT / "data/derived/fig6_comparison"
GRAPH = ROOT / "data/raw/training_graph/training_graph_edges.parquet"
REF = ROOT / "data/derived/degree_controls/degree_controls_table.tsv"
SEED = 4707
# system key (columns prefix in per_run_scores.tsv) -> name in the table; the order is the table's order
SYSTEMS = [("DisAtlas", "DisAtlas"), ("DisAtlas node-disjoint", "DisAtlas node-disjoint"),
           ("Claude-Chat Opus 5.5", "Claude Chat"), ("Claude-Sci Opus 5.5", "Claude Sci"), ("Phylo", "Phylo Biomni"),
           ("Gemini Pro 3.1", "Gemini Pro 3.1"), ("TxGNN published ckpt", "TxGNN (published checkpoint)")]
ND = "DisAtlas node-disjoint"


def training_degree():
    g = pd.read_parquet(GRAPH)
    return pd.concat([g.u, g.v]).value_counts().to_dict()


def system_scores():
    R = pd.read_csv(FIG6 / "per_run_scores.tsv", sep="\t"); N = pd.read_csv(FIG6 / "node_disjoint_per_seed.tsv", sep="\t")
    assert (R.blind_id == N.blind_id).all() and (R.y == N.y).all()
    S = {}
    for key, _ in SYSTEMS:
        if key == "DisAtlas": S[key] = R["DisAtlas"].to_numpy(float)
        elif key == ND: S[key] = N[[c for c in N.columns if c.startswith(ND + "|")]].mean(axis=1, skipna=False).to_numpy(float)
        else: S[key] = R[[c for c in R.columns if c.startswith(key + "|")]].mean(axis=1, skipna=False).to_numpy(float)
    return R[["blind_id", "a", "b", "y", "pool"]], S


def draws(P, B):
    nodes = pd.unique(P[["a", "b"]].values.ravel()); ix = {n: i for i, n in enumerate(nodes)}
    ai = P.a.map(ix).values; bi = P.b.map(ix).values
    rg = np.random.default_rng(SEED); W = np.empty((B, len(P)), np.int32); WN = np.empty((B, len(nodes)), np.int16)
    for k in range(B):
        w = np.bincount(rg.integers(0, len(nodes), len(nodes)), minlength=len(nodes)); WN[k] = w; W[k] = w[ai] * w[bi]
    return W, WN, {n: k for k, n in enumerate(nodes)}


def wauc_boot(y, s, W):
    y = np.asarray(y).astype(bool); s = np.asarray(s, float)
    u, g = np.unique(s, return_inverse=True); G = len(u); out = np.full(W.shape[0], np.nan)
    for k in range(W.shape[0]):
        w = W[k].astype(float)
        wp = np.bincount(g[y], w[y], G); wn = np.bincount(g[~y], w[~y], G)
        P_, N_ = wp.sum(), wn.sum()
        if P_ == 0 or N_ == 0: continue
        cn = np.cumsum(wn) - wn
        out[k] = (wp * (cn + 0.5 * wn)).sum() / (P_ * N_)
    return out


def ci(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))


def per_disease_auc(D, s):
    s = np.asarray(s, float)
    L = pd.DataFrame({"d": np.r_[D.a.values, D.b.values], "y": np.r_[D.y.values, D.y.values], "s": np.r_[s, s]})
    out = {}
    for d, g in L.groupby("d", sort=True):
        if g.y.nunique() == 2: out[d] = roc_auc_score(g.y, g.s)
    return pd.Series(out)


def degree_matched_subsample(pos, neg, deg, n_bins=10, seed=SEED):
    """Quantile-binned degree-product matching of positive and negative pairs (sets of frozensets)."""
    rng = np.random.default_rng(seed)
    dp = lambda p: np.log1p(deg.get(tuple(p)[0], 0)) + np.log1p(deg.get(tuple(p)[1], 0))
    key = lambda p: tuple(sorted(map(str, tuple(p))))
    pos, neg = sorted(pos, key=key), sorted(neg, key=key)
    pv, nv = np.array([dp(p) for p in pos]), np.array([dp(p) for p in neg])
    edges = np.quantile(np.r_[pv, nv], np.linspace(0, 1, n_bins + 1)); edges[0], edges[-1] = -np.inf, np.inf
    pi, ni = np.digitize(pv, edges[1:-1]), np.digitize(nv, edges[1:-1])
    kp, kn = [], []
    for b in range(n_bins):
        ip, inn = np.where(pi == b)[0], np.where(ni == b)[0]; k = min(len(ip), len(inn))
        if k == 0: continue
        kp += [pos[i] for i in rng.choice(ip, size=k, replace=False)]
        kn += [neg[i] for i in rng.choice(inn, size=k, replace=False)]
    return set(kp), set(kn)


def dm_select(D, deg):
    keys = [frozenset((x, z)) for x, z in zip(D.a.values, D.b.values)]; first = {}
    for i, k in enumerate(keys): first.setdefault(k, i)
    pos = {k for i, k in enumerate(keys) if D.y.values[i] == 1 and first[k] == i}
    neg = {k for i, k in enumerate(keys) if D.y.values[i] == 0 and first[k] == i}
    conf = pos & neg; pos -= conf; neg -= conf
    pm, nm = degree_matched_subsample(pos, neg, deg)
    return np.array(sorted(first[k] for k in (pm | nm)))


def compute(B=2000):
    P, S = system_scores(); deg = training_degree()
    gl = np.log1p(P.a.map(deg).fillna(0).values) + np.log1p(P.b.map(deg).fillna(0).values)
    W, WN, nix = draws(P, B)

    def sys_score(k):
        s = S[k]
        return s if k == ND else np.where(np.isfinite(s), s, np.nanmin(s) - 1e-9)

    def within(D, s, ok):
        pdv = per_disease_auc(D[ok].reset_index(drop=True), s[ok])
        di = np.array([nix[d] for d in pdv.index]); wv = WN[:, di].astype(float)
        return float(pdv.mean()), (wv * pdv.values).sum(1) / np.where(wv.sum(1) > 0, wv.sum(1), np.nan), len(pdv)

    rows = []
    for case in ("all", "test"):
        idx = np.where((P.pool == "test").values if case == "test" else np.ones(len(P), bool))[0]
        D = P.iloc[idx].reset_index(drop=True); y = D.y.values; Wd = W[:, idx]; g = gl[idx]
        selm = np.zeros(len(D), bool); selm[dm_select(D, deg)] = True
        sr = sys_score(ND)[idx]; okr = np.isfinite(sr); m = selm & okr
        dg = roc_auc_score(y[m], g[m]); dgb = wauc_boot(y[m], g[m], Wd[:, m])
        br = wauc_boot(y[m], sr[m], Wd[:, m]); wr, wbr, _ = within(D, sr, okr)
        for key, name in SYSTEMS:
            s = sys_score(key)[idx]
            a = roc_auc_score(y[m], s[m]); b = wauc_boot(y[m], s[m], Wd[:, m]); wv, wb, wn = within(D, s, okr)
            r = dict(system=name, case_set=case, n_pairs=int(okr.sum()), n_pos=int(y[okr].sum()), n_diseases=wn, dm_pairs=int(m.sum()),
                     dm_pos=int(y[m].sum()), dm_auroc=a, dm_ci_lo=ci(b)[0], dm_ci_hi=ci(b)[1], dm_degree_alone_auroc=dg,
                     dm_degree_alone_ci_lo=ci(dgb)[0], dm_degree_alone_ci_hi=ci(dgb)[1], within_auroc=wv, within_ci_lo=ci(wb)[0],
                     within_ci_hi=ci(wb)[1])
            if key != ND:
                dd = br - b; dw = wbr - wb
                r.update(dm_delta_heldout_minus_system=roc_auc_score(y[m], sr[m]) - a, dm_delta_ci_lo=ci(dd)[0], dm_delta_ci_hi=ci(dd)[1],
                         within_delta_heldout_minus_system=wr - wv, within_delta_ci_lo=ci(dw)[0], within_delta_ci_hi=ci(dw)[1])
            rows.append(r)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    B = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    compute(B).to_csv(out / "degree_controls_table.tsv", sep="\t", index=False)
    print("written", out / "degree_controls_table.tsv")
