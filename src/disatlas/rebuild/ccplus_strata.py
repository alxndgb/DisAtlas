"""Rebuild data/derived/si_rare_common/ccplus_strata.tsv (Supplementary Table SI_rare_common_strata: DisAtlas AUROC on
IxIDN-Matched and IxIDN-Matched-CC+ by rarity stratum, against the matched-shortcut stack).

    python -m disatlas.rebuild.ccplus_strata OUTDIR [B]

Inputs: data/derived/si_rare_common/ixidn_ccplus_pairs.tsv (pairs, labels, strata, shortcut features),
data/zenodo/model_outputs/cosine_matrix/{cosine_ens5_13021.npy, ids_13021.json} (DisAtlas score = five-seed mean
cosine), data/derived/fig6_comparison/per_run_scores.tsv (order of the 938 IxIDN-Matched diseases for the panel-c draws).
Shortcut stack: logistic regression on dsum, mdeg, mondo_simgic, name_tfidf_cos, genecov (median-imputed, z-scored),
5-fold out-of-fold predictions (StratifiedKFold shuffle, random_state 0), refitted within each row's pair set.
Intervals: disease-cluster bootstrap, pair weight w_a * w_b, B = 2000, draws without both classes skipped.
  IxIDN-Matched, All pairs: the panel-c draws (938 diseases in per_run_scores order, 10 shards of B/10,
  shard s seeded default_rng([4707, s])), shared by DisAtlas, the stack and their difference.
  Other rows: clusters = the diseases of the row's pair set; one stream default_rng(4707) for DisAtlas and a separate
  stream default_rng(4707) for the paired difference.
Rows with fewer than 20 pairs in one class are not scored. Test: tests/test_rebuild_ccplus_strata.py.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict, StratifiedKFold

ROOT = Path(__file__).resolve().parents[3]
MO = ROOT / "data/zenodo/model_outputs"
STACK = ["dsum", "mdeg", "mondo_simgic", "name_tfidf_cos", "genecov"]
STRATA = ["all", "rare-rare", "rare-common", "common-common"]


def load():
    X = pd.read_csv(ROOT / "data/derived/si_rare_common/ixidn_ccplus_pairs.tsv", sep="\t")
    ids = json.load(open(MO / "cosine_matrix/ids_13021.json")); ix = {d: i for i, d in enumerate(ids)}
    C = np.load(MO / "cosine_matrix/cosine_ens5_13021.npy", mmap_mode="r")
    X["p"] = [float(C[ix[a], ix[b]]) for a, b in zip(X.a, X.b)]
    X["y"] = X.label
    return X


def stack_oof(d):
    Xs = d[STACK].astype(float); Xs = Xs.fillna(Xs.median()); Xs = (Xs - Xs.mean()) / (Xs.std() + 1e-9)
    return cross_val_predict(LogisticRegression(max_iter=1000), Xs.values, d.label.values,
                             cv=StratifiedKFold(5, shuffle=True, random_state=0), method="predict_proba")[:, 1]


def _weights_stream(n, B, rng):
    for _ in range(B):
        yield np.bincount(rng.integers(0, n, n), minlength=n)


def _boot(d, nodes, weights, s1, s2=None):
    ix = {n: i for i, n in enumerate(nodes)}; ai, bi = d.a.map(ix).values, d.b.map(ix).values; y = d.y.values
    out = []
    for w in weights:
        ww = (w[ai] * w[bi]).astype(float); m = ww > 0
        if len(np.unique(y[m])) < 2: continue
        a1 = roc_auc_score(y[m], s1[m], sample_weight=ww[m])
        out.append(a1 if s2 is None else (a1, a1 - roc_auc_score(y[m], s2[m], sample_weight=ww[m])))
    return np.asarray(out)


def compute(B=2000):
    X = load(); fz = (X.source == "IxIDN-Matched").values
    order = pd.read_csv(ROOT / "data/derived/fig6_comparison/per_run_scores.tsv", sep="\t", usecols=["a", "b"])
    panel = pd.unique(order[["a", "b"]].values.ravel()); assert len(panel) == 938
    WB = []
    for s in range(10):
        rg = np.random.default_rng([4707, s]); WB += list(_weights_stream(len(panel), B // 10, rg))
    rows = []
    for bench, D0 in [("IxIDN-Matched", X[fz]), ("IxIDN-Matched-CC+", X)]:
        for scope, D1 in [("All pairs", D0), ("Test split pairs", D0[D0.pool == "test"])]:
            for st in STRATA:
                d = (D1 if st == "all" else D1[D1.stratum == st]).copy()
                npos, nneg = int(d.y.sum()), int((d.y == 0).sum())
                r = dict(benchmark=bench, scope=scope, stratum=st, n_pairs=len(d), n_pos=npos, n_neg=nneg,
                         n_added_neg=int((d.source == "CC+ added").sum()), n_diseases=len(pd.unique(d[["a", "b"]].values.ravel())))
                if min(npos, nneg) < 20:
                    r["note"] = "fewer than 20 pairs in one class: AUROC not reported"; rows.append(r); continue
                d["stk"] = stack_oof(d); p, k = d.p.values, d.stk.values
                r["auroc"] = roc_auc_score(d.y, p); r["shortcut_stack_matched"] = roc_auc_score(d.y, k)
                if bench == "IxIDN-Matched" and scope == "All pairs":
                    o = _boot(d, panel, WB, p, k); bp, bd = o[:, 0], o[:, 1]
                    r["ci_bootstrap"] = f"panel-c disease-cluster draws (938 clusters, B={B})"
                else:
                    nodes = pd.unique(d[["a", "b"]].values.ravel())
                    bp = _boot(d, nodes, _weights_stream(len(nodes), B, np.random.default_rng(4707)), p)
                    bd = _boot(d, nodes, _weights_stream(len(nodes), B, np.random.default_rng(4707)), p, k)[:, 1]
                    r["ci_bootstrap"] = f"disease-cluster bootstrap, B={B}, seed 4707"
                r["ci_lo"], r["ci_hi"] = np.percentile(bp, [2.5, 97.5])
                r["delta_vs_stack"] = r["auroc"] - r["shortcut_stack_matched"]
                r["delta_vs_stack_lo"], r["delta_vs_stack_hi"] = np.percentile(bd, [2.5, 97.5])
                r["delta_vs_stack_p_le0"] = float((bd <= 0).mean()); r["n_boot"] = float(len(bp))
                rows.append(r)
    cols = ["benchmark", "scope", "stratum", "n_pairs", "n_pos", "n_neg", "n_added_neg", "n_diseases", "auroc",
            "shortcut_stack_matched", "ci_lo", "ci_hi", "delta_vs_stack", "delta_vs_stack_lo", "delta_vs_stack_hi",
            "delta_vs_stack_p_le0", "n_boot", "ci_bootstrap", "note"]
    return pd.DataFrame(rows).reindex(columns=cols)


if __name__ == "__main__":
    o = Path(sys.argv[1]); o.mkdir(parents=True, exist_ok=True)
    compute(int(sys.argv[2]) if len(sys.argv) > 2 else 2000).to_csv(o / "ccplus_strata.tsv", sep="\t", index=False)
    print("wrote", o / "ccplus_strata.tsv")
