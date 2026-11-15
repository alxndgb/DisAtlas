"""Rebuild data/derived/si_baselines/baselines_auroc_B2000.tsv (Supplementary Table S5: DisAtlas cosine against
training-free baselines). Port of the authors' ga_contrasts_h8_h9.py (in-graph block A, own-edges contrasts) and
cos_heldout_baselines_ga.py (held-out block).

Every quantity: g1 disease-cluster bootstrap, B = 2,000, fresh default_rng(4707) per quantity; node set = diseases of the
pairs in that quantity; each draw resamples diseases with replacement and weights pair (a, b) by w_a * w_b (draws with one
class only are skipped); weighted roc_auc_score; d = DisAtlas minus score, paired over the same draws; 95% percentile
intervals. DisAtlas: c (in-graph) = data/derived/fig6_comparison/per_run_scores.tsv column DisAtlas; nd_c (diseases
held out of training) = mean of the five per-seed node-disjoint cosines in node_disjoint_per_seed.tsv.
In-graph baselines (CN, Adamic-Adar, weighted AA, degree product, TF-IDF graph SVD, weighted / binary adjacency SVD) are
computed here from data/raw/training_graph/training_graph_edges.parquet (agree with the authors' scores to 1e-13).
Annotation Jaccards and the held-out (node-disjoint) block: data/derived/si_baselines/baseline_pair_scores.tsv.gz (staged). The two MONDO simGIC rows (AUROC only, no
interval) are copied from the shipped table: they come from the ontology baseline (disatlas.ontology_baseline).

    python -m disatlas.rebuild.si_baselines OUT.tsv      (about 15 min on 4 CPUs)
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[3]
D = ROOT / "data/derived"
B = 2000; SEED = 4707
IG = ["tfidf_graph_svd_cosine__full_graph", "adjacency_svd_binary_cosine", "adjacency_svd_weighted_cosine", "adamic_adar__full_graph",
      "common_neighbours__full_graph", "weighted_adamic_adar__full_graph", "degree_product__full_graph", "allfeat_jaccard", "gene_jaccard", "hpo_jaccard"]
HO = ["tfidf_svd_place55", "aa_place55", "waa_place55", "cn_place55", "tfidf_graph_svd_cosine", "adamic_adar", "common_neighbours", "allfeat_jaccard"]
OE = ["common_neighbours", "adamic_adar", "tfidf_graph_svd_cosine"]


GRAPH = ROOT / "data/raw/training_graph/training_graph_edges.parquet"
GRAPH_COLS = ["tfidf_graph_svd_cosine__full_graph", "adjacency_svd_weighted_cosine", "adjacency_svd_binary_cosine", "adamic_adar__full_graph",
              "common_neighbours__full_graph", "weighted_adamic_adar__full_graph", "degree_product__full_graph"]


def graph_scores(a_ids, b_ids):
    """In-graph training-free scores on the released training graph (benchmark and FDA pairs stripped), as the authors'
    heuristic_check heuristics.py graph_scores(set()) (common neighbours, Adamic-Adar, weighted Adamic-Adar, degree product,
    32-d TF-IDF graph SVD cosine) and full_setting_baselines (32-d truncated SVD of the w_arm-weighted / binary adjacency)."""
    import scipy.sparse as sp
    from sklearn.decomposition import TruncatedSVD
    E = pd.read_parquet(GRAPH, columns=["u", "v", "w_arm"])
    V = sorted(set(E.u) | set(E.v)); N = len(V); ix = {d: i for i, d in enumerate(V)}
    ui, vi, w = E.u.map(ix).values, E.v.map(ix).values, E.w_arm.values.astype(float)
    A = sp.csr_matrix((np.r_[w, w], (np.r_[ui, vi], np.r_[vi, ui])), shape=(N, N)); A.sum_duplicates()
    Ab = (A > 0).astype(np.float64).tocsr()
    ia = np.array([ix[d] for d in a_ids]); ib = np.array([ix[d] for d in b_ids])
    NB = [set(Ab.indices[Ab.indptr[i]:Ab.indptr[i + 1]]) for i in range(N)]
    Wd = [dict(zip(A.indices[A.indptr[i]:A.indptr[i + 1]], A.data[A.indptr[i]:A.indptr[i + 1]])) for i in range(N)]
    deg = np.asarray(Ab.sum(1)).ravel(); stren = np.asarray(A.sum(1)).ravel(); n = len(ia)
    cn, aa, waa, dp = (np.zeros(n) for _ in range(4))
    for i in range(n):
        c = NB[ia[i]] & NB[ib[i]]; cn[i] = len(c)
        aa[i] = sum(1.0 / np.log(deg[z]) for z in c if deg[z] > 1)
        waa[i] = sum((Wd[ia[i]][z] + Wd[ib[i]][z]) / np.log1p(stren[z]) for z in c)
        dp[i] = deg[ia[i]] * deg[ib[i]]
    idf = np.log(N / (1 + np.asarray((A > 0).sum(0)).ravel()))

    def tfidf_rows(M):
        M = M.copy(); M.data = np.log1p(M.data); M = (M @ sp.diags(idf)).tocsr()
        nr = np.sqrt(np.asarray(M.multiply(M).sum(1)).ravel()); return sp.diags(1 / np.where(nr > 0, nr, 1)) @ M

    def unitcos(X):
        nr = np.linalg.norm(X, axis=1, keepdims=True); X = X / np.where(nr > 0, nr, 1); return (X[ia] * X[ib]).sum(1)

    T = tfidf_rows(A); svd = TruncatedSVD(32, algorithm="arpack", random_state=0).fit(T)
    out = {"tfidf_graph_svd_cosine__full_graph": unitcos(T @ svd.components_.T)}
    for k, M in (("adjacency_svd_weighted_cosine", A), ("adjacency_svd_binary_cosine", Ab)):
        out[k] = unitcos(TruncatedSVD(32, algorithm="arpack", random_state=0).fit_transform(M))
    out.update({"adamic_adar__full_graph": aa, "common_neighbours__full_graph": cn, "weighted_adamic_adar__full_graph": waa, "degree_product__full_graph": dp})
    return out


def pairs():
    S = pd.read_csv(D / "si_baselines/baseline_pair_scores.tsv.gz", sep="\t")
    P = pd.read_csv(D / "fig6_comparison/per_run_scores.tsv", sep="\t", usecols=["blind_id", "a", "b", "y", "pool", "DisAtlas"])
    N = pd.read_csv(D / "fig6_comparison/node_disjoint_per_seed.tsv", sep="\t")
    N = pd.DataFrame({"blind_id": N.blind_id, "nd_c": N[[c for c in N.columns if "|s" in c]].mean(axis=1)})
    S = S.merge(P.rename(columns={"DisAtlas": "c"}), on=["blind_id", "a", "b", "y"], how="left", validate="1:1").merge(N, on="blind_id", how="left", validate="1:1")
    assert len(S) == 4364 and S.c.notna().all() and S.nd_c.notna().sum() == 4357 and (S.pool == "test").sum() == 1094
    for k, v in graph_scores(S.a.astype(str).values, S.b.astype(str).values).items():
        S[k] = v
    return S


def boot(D_, cols, ref):
    D_ = D_.reset_index(drop=True); y = D_.y.values; nodes = pd.unique(D_[["a", "b"]].values.ravel()); ix = {n: i for i, n in enumerate(nodes)}
    ai, bi = D_.a.map(ix).values, D_.b.map(ix).values; rg = np.random.default_rng(SEED); bs = {c: [] for c in cols}
    for _ in range(B):
        w = np.bincount(rg.integers(0, len(nodes), len(nodes)), minlength=len(nodes)); ww = (w[ai] * w[bi]).astype(float); q = ww > 0
        if len(set(y[q])) < 2: continue
        for c in cols: bs[c].append(roc_auc_score(y[q], D_[c].values[q], sample_weight=ww[q]))
    out = []
    for c in cols:
        b = np.array(bs[c]); r = dict(score=c, auc=roc_auc_score(y, D_[c].values), lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5))
        if c != ref:
            d = np.array(bs[ref]) - b; r.update(d=roc_auc_score(y, D_[ref].values) - r["auc"], dlo=np.percentile(d, 2.5), dhi=np.percentile(d, 97.5))
        out.append(r)
    return out


def build():
    S = pairs(); test = (S.pool == "test").values; ok = S.nd_c.notna().values; rows = []
    for sub, m in (("all", np.ones(len(S), bool)), ("test", test)):
        rows += [dict(block="in_graph", subset=sub, **r) for r in boot(S[m], ["c"] + IG, "c")]
    def held(m, cols):
        X = S.loc[m, ["a", "b", "y", "nd_c"]].copy()
        for c in cols: X[c] = S.loc[m, f"ho_{c}"].values
        return X
    for sub, m in (("all", ok), ("test", ok & test)):
        rows += [dict(block="held_out", subset=sub, **r) for r in boot(held(m, HO), ["nd_c"] + HO, "nd_c")]
    rows += [dict(block="held_out_own_edges", subset="all", **r) for r in boot(held(ok, OE), ["nd_c"] + OE, "nd_c") if r["score"] != "nd_c"]
    T = pd.DataFrame(rows)[["block", "subset", "score", "auc", "lo", "hi", "d", "dlo", "dhi"]]
    REF = pd.read_csv(D / "si_baselines/baselines_auroc_B2000.tsv", sep="\t")
    return pd.concat([T, REF[REF.score == "mondo_simgic"]], ignore_index=True)


if __name__ == "__main__":
    build().to_csv(sys.argv[1], sep="\t", index=False, float_format="%.6f")
