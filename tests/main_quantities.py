"""Quantities behind the numbers in the DisAtlas main text, computed from repository data only.

Every function reads files under data/derived (or the generated tables under figures/main, tables/main)
and returns plain numbers.  tests/test_main_parity.py compares them with the manuscript tokens listed in
tests/main_text_numbers.tsv.  `all_quantities()` returns one flat dict {key: value}.
"""
import functools, glob, os, gzip
import numpy as np, pandas as pd
from scipy.stats import rankdata

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DD = os.path.join(ROOT, "data", "derived")
B_BOOT = int(os.environ.get("PARITY_B", "2000"))


def _tsv(*p, **kw):
    return pd.read_csv(os.path.join(DD, *p), sep="\t", **kw)


def auc(y, s, w=None):
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, s, sample_weight=w))


# ---------------------------------------------------------------- released five-seed cosine matrix (Zenodo)
@functools.lru_cache(None)
def cosine_matrix():
    import json
    MO = os.path.join(ROOT, "data", "zenodo", "model_outputs", "cosine_matrix")
    ids = json.load(open(os.path.join(MO, "ids_13021.json")))
    return ids, {d: i for i, d in enumerate(ids)}, np.load(os.path.join(MO, "cosine_ens5_13021.npy"), mmap_mode="r")


@functools.lru_cache(None)
def all_pair_cosines():
    ids, ix, C = cosine_matrix(); n = len(ids)
    return np.sort(np.concatenate([np.asarray(C[i, i + 1:], dtype=np.float64) for i in range(n - 1)]))


def pair_cos(a, b):
    ids, ix, C = cosine_matrix(); return float(C[ix[a], ix[b]])


def pct_all_pairs(s):
    P = all_pair_cosines(); return 100 * np.searchsorted(P, s, "left") / len(P)


def seed_mean_cos(a, b, arm):
    import glob
    v = []
    for f in sorted(glob.glob(os.path.join(ROOT, "data", "zenodo", "model_outputs", "embeddings", arm, "embedding_s*.npz"))):
        z = np.load(f); ix = {d: i for i, d in enumerate(z["ids"])}; v.append(float(z["X"][ix[a]] @ z["X"][ix[b]]))
    assert len(v) == 5, (arm, len(v)); return float(np.mean(v))


# ---------------------------------------------------------------- disease set and map (Figure 5)
@functools.lru_cache(None)
def disease_map():
    return _tsv("fig5_clustering", "disease_map_k55.tsv.gz")


def q_disease_set():
    m = disease_map(); n = len(m); rare = int(m.is_rare.sum())
    q = dict(n_diseases=n, n_rare=rare, n_common=n - rare, pct_common=100 * (n - rare) / n,
             n_clusters=int(m.cluster.nunique()))
    # clusters in which common diseases make up at least 10% of the members (every cluster holds rare diseases)
    t = m.groupby("cluster").is_rare.agg(["mean", "size"])
    n_common_c = m.groupby("cluster").is_rare.apply(lambda r: int((r == 0).sum()))
    mixed = t[10 * n_common_c >= t["size"]].index  # exact integer test (a 27/270 cluster sits at 10.0%)
    q["mixed_clusters"] = len(mixed)
    q["mixed_pct_common"] = 100 * ((m.is_rare == 0) & m.cluster.isin(mixed)).sum() / (m.is_rare == 0).sum()
    q["mixed_pct_rare"] = 100 * ((m.is_rare == 1) & m.cluster.isin(mixed)).sum() / (m.is_rare == 1).sum()
    q.update({f"cluster_size|{c}": int(s) for c, s in t["size"].items()})
    q.update({f"cluster|{i}": int(c) for i, c in zip(m.mondo_id, m.cluster)})
    return q


def q_cocluster():
    m = disease_map().set_index("mondo_id").cluster
    P = _tsv("fig5_clustering", "ixidn_matched_pair_scores.tsv"); T = P[P.pool == "test"]
    same = (T.a.map(m).values == T.b.map(m).values)
    pos, neg = same[T.y.values == 1], same[T.y.values == 0]
    return dict(test_pairs=len(T), test_pos=int(len(pos)), test_neg=int(len(neg)),
                cocluster_pos_pct=100 * pos.mean(), cocluster_neg_pct=100 * neg.mean(),
                cocluster_ratio=pos.mean() / neg.mean(),
                test_diseases=int(pd.unique(T[["a", "b"]].values.ravel()).size))


def q_strata():
    P = _tsv("fig5_clustering", "ixidn_matched_pair_scores.tsv"); q = {}
    for s, g in P.groupby("rare_common_pair"):
        q[f"stratum_auroc|{s}"] = auc(g.y, g.p); q[f"stratum_npos|{s}"] = int(g.y.sum()); q[f"stratum_nneg|{s}"] = int((g.y == 0).sum())
    q["full_auroc_pairscores"] = auc(P.y, P.p)
    return q


# ---------------------------------------------------------------- control pairs
def q_controls():
    """Positive and negative control pairs: five-seed mean cosine and its percentile among all disease pairs."""
    C = _tsv("control_pairs", "control_pairs.tsv")
    m = disease_map().set_index("mondo_id").cluster
    C["split"] = C.a.map(m).values != C.b.map(m).values
    C["cos"] = [pair_cos(a, b) for a, b in zip(C.a, C.b)]; C["pct"] = [pct_all_pairs(c) for c in C.cos]
    q = {f"ctrl_p|{t}": float(c) for t, c in zip(C.tag, C.cos)}
    q.update({f"ctrl_pct|{t}": float(v) for t, v in zip(C.tag, C.pct)})
    pos, neg = C[C.role == "pos"], C[C.role == "neg"]
    q.update(n_controls=len(C), n_pos_controls=len(pos), n_neg_controls=len(neg),
             n_pos_split=int(pos.split.sum()), n_neg_split=int(neg.split.sum()),
             min_pos_pct=float(pos.pct.min()), min_pos_p=float(pos.cos.min()))
    # illustrative pairs of the Results cosine paragraph (five-seed mean of per-seed cosines)
    for t, a, b in (("Crohn-UC", "MONDO:0005011", "MONDO:0005101"), ("DI-T2D", "MONDO:0004782", "MONDO:0005148")):
        q[f"illus_cos|{t}"] = seed_mean_cos(a, b, "gene_anchored")
    P = _tsv("fig5_clustering", "ixidn_matched_pair_scores.tsv"); T = P[P.pool == "test"]
    q["median_cos_test_pos"] = float(T[T.y == 1].p.median()); q["median_cos_test_neg"] = float(T[T.y == 0].p.median())
    # rare-common recovery: mean percentile of the co-investigated partner among all diseases of the other class
    ids, ix, M = cosine_matrix(); r = disease_map().set_index("mondo_id").is_rare.reindex(ids).values
    com, rar = np.where(r == 0)[0], np.where(r == 1)[0]; f, g = [], []
    RC = T[(T.y == 1) & (T.rare_common_pair == "rare-common")]
    for a, b in zip(RC.a, RC.b):
        ra, co = (a, b) if r[ix[a]] == 1 else (b, a)
        row = np.asarray(M[ix[ra]]); f.append((row[com] < row[ix[co]]).sum() / (len(com) - 1))
        row = np.asarray(M[ix[co]]); g.append((row[rar] < row[ix[ra]]).sum() / (len(rar) - 1))
    q.update(rc_test_pairs=len(f), rc_mean_pct_rare_anchor=float(np.mean(f)), rc_mean_pct_common_anchor=float(np.mean(g)))
    return q


# ---------------------------------------------------------------- IxIDN and IxIDN-Matched (Figure 3, Methods)
def q_ixidn():
    F = _tsv("fig3_ixidn", "funnel_counts.tsv")
    S1 = _tsv("fig3_ixidn", "stage1_ixidn_pos.tsv.gz")
    X = _tsv("fig3_ixidn", "ixidn_matched_pairs.tsv")
    V = _tsv("fig3_ixidn", "matched_candidate_verdicts.tsv").verdict.str.upper().value_counts()
    SD = pd.read_csv(os.path.join(DD, "fig3_ixidn", "shortcut_diagnostics.csv"))
    U = _tsv("fig3_ixidn", "stage2_universe_diseases.tsv.gz")
    m = disease_map().set_index("mondo_id").is_rare
    dz = pd.unique(X[["a", "b"]].values.ravel()); common = int((m.reindex(dz) == 0).sum())
    te = X[X.pool == "test"]; ted = pd.unique(te[["a", "b"]].values.ravel())
    sc = SD[SD.matched_by_designA.astype(str) == "True"].expanded_auroc
    q = dict(ixidn_links=len(S1), ixidn_diseases=int(pd.unique(S1.values.ravel()).size), universe_diseases=len(U),
             matched_candidates=int(V.sum()), verdict_strong=int(V.get("STRONG", 0)), verdict_plausible=int(V.get("PLAUSIBLE", 0)),
             bench_pairs=len(X), bench_pos=int(X.label.sum()), bench_neg=int((X.label == 0).sum()),
             bench_diseases=len(dz), bench_common=common, bench_pct_common=100 * common / len(dz),
             test_split_pairs=len(te), test_split_pos=int(te.label.sum()), test_split_neg=int((te.label == 0).sum()),
             test_half_diseases=int((_tsv("fig3_ixidn", "disease_split.tsv").split == "test").sum()),
             test_pair_diseases=len(ted), shortcut_auroc_min=float(sc.min()), shortcut_auroc_max=float(sc.max()))
    q.update({f"shortcut_auroc|{f}": float(v) for f, v in zip(SD.feature, SD.expanded_auroc)})
    q["ixidn_diseases_in_disatlas"] = int(np.isin(pd.unique(S1.values.ravel()), disease_map().mondo_id).sum())
    q["funnel_pinned_pos"] = int(F[(F.stage == "3_pinned_universe") & (F["class"] == "IxIDN+")].pairs.iloc[0])
    return q


# ---------------------------------------------------------------- cross-system comparison (Figure 6, Table 1)
@functools.lru_cache(None)
def per_run():
    return _tsv("fig6_comparison", "per_run_scores.tsv")


def _system_mean(R, prefix):
    cols = [c for c in R.columns if c.startswith(prefix + "|")]
    return R[cols].mean(axis=1, skipna=True), cols


def q_comparison():
    R = per_run(); y = R.y.values
    ND = _tsv("fig6_comparison", "node_disjoint_per_seed.tsv").set_index("blind_id").loc[R.blind_id]
    nd = ND.filter(like="DisAtlas node-disjoint|").mean(axis=1).values
    full = R["DisAtlas"].values; ok = ~np.isnan(nd); te = (R.pool == "test").values
    q = dict(cmp_pairs=len(R), cmp_pos=int(y.sum()), nd_unscored=int((~ok).sum()), nd_scored=int(ok.sum()), nd_scored_test=int((ok & te).sum()),
             auroc_full=auc(y, full), auroc_nd=auc(y[ok], nd[ok]), auroc_full_on_nd_pairs=auc(y[ok], full[ok]),
             auroc_full_test=auc(y[te], full[te]), test_pairs=int(te.sum()),
             topk=int(round(0.1 * len(R))), jaccard_random=0.1 / (2 - 0.1))
    q["delta_nd_minus_full"] = q["auroc_nd"] - q["auroc_full_on_nd_pairs"]
    for nm, pre in [("claude_chat", "Claude-Chat Opus 5.5"), ("claude_sci", "Claude-Sci Opus 5.5"), ("gemini", "Gemini Pro 3.1"),
                    ("phylo", "Phylo"), ("txgnn", "TxGNN published ckpt")]:
        s, cols = _system_mean(R, pre); q[f"runs|{nm}"] = len(cols); q[f"unscored|{nm}"] = int(s.isna().sum())
    q["runs|disatlas_seeds"] = len([c for c in R.columns if c.startswith("DisAtlas|s")])
    # disease-cluster bootstrap: each benchmark disease resampled with all its pairs, pair weight w_a*w_b
    nodes = pd.unique(R[["a", "b"]].values.ravel()); ix = {d: i for i, d in enumerate(nodes)}
    ai, bi = R.a.map(ix).values, R.b.map(ix).values
    rg = np.random.default_rng(4707); bs = {k: [] for k in ("full", "nd", "delta", "full_test")}
    for _ in range(B_BOOT):
        w = np.bincount(rg.integers(0, len(nodes), len(nodes)), minlength=len(nodes)); ww = (w[ai] * w[bi]).astype(float)
        bs["full"].append(auc(y, full, ww))
        k = ok & (ww > 0); a_nd = auc(y[k], nd[k], ww[k]); a_fu = auc(y[k], full[k], ww[k])
        bs["nd"].append(a_nd); bs["delta"].append(a_nd - a_fu)
    # test split: same scheme, restricted to the test pairs and their own diseases (seed 11)
    yt, ft, tn = y[te], full[te], pd.unique(R[te][["a", "b"]].values.ravel()); tix = {d: i for i, d in enumerate(tn)}
    ta, tb = R[te].a.map(tix).values, R[te].b.map(tix).values; rg = np.random.default_rng(11)
    for _ in range(B_BOOT):
        w = np.bincount(rg.integers(0, len(tn), len(tn)), minlength=len(tn)); ww = (w[ta] * w[tb]).astype(float); k = ww > 0
        if len(set(yt[k])) == 2: bs["full_test"].append(auc(yt[k], ft[k], ww[k]))
    for k, v in bs.items():
        q[f"{k}_lo"], q[f"{k}_hi"] = (float(x) for x in np.percentile(v, [2.5, 97.5]))
    # sensitivity: resample the 55 DisAtlas disease clusters instead of single diseases
    cl = disease_map().set_index("mondo_id").cluster; ca, cb = R.a.map(cl).values, R.b.map(cl).values
    C = np.sort(cl.unique()); cpos = {c: i for i, c in enumerate(C)}; ca = np.array([cpos[c] for c in ca]); cb = np.array([cpos[c] for c in cb])
    rg = np.random.default_rng(4707); cs = {k: [] for k in ("full", "nd", "delta")}
    for _ in range(B_BOOT):
        w = np.bincount(rg.integers(0, len(C), len(C)), minlength=len(C)); ww = (w[ca] * w[cb]).astype(float)
        if len(set(y[ww > 0])) < 2: continue
        cs["full"].append(auc(y, full, ww)); k = ok & (ww > 0)
        a_nd = auc(y[k], nd[k], ww[k]); cs["nd"].append(a_nd); cs["delta"].append(a_nd - auc(y[k], full[k], ww[k]))
    for k, v in cs.items():
        q[f"clboot_{k}_lo"], q[f"clboot_{k}_hi"] = (float(x) for x in np.percentile(v, [2.5, 97.5]))
    return q


def q_run_agreement():
    """Median pairwise Spearman correlation between repeated runs of each system on IxIDN-Matched."""
    A = pd.read_csv(os.path.join(DD, "fig6_comparison", "all_methods_run_agreement.csv"))
    return {f"runagree_median|{k}": float(v) for k, v in A.groupby("method").spearman.median().items()}



def q_bootstrap_ci():
    """Disease-cluster bootstrap AUROC and 95% CI per scope and system (Figure 6, Supplementary)."""
    B = _tsv("fig6_comparison", "disease_bootstrap_auroc_ci.tsv"); q = {}
    for _, r in B.iterrows():
        sy = "nd" if "node-disjoint" in r.system else ("full" if r.system == "DisAtlas (full)" else r.system)
        for c in ("auroc", "lo", "hi", "n_pairs", "n_pos"):
            q[f"boot_{c}|{r.scope}|{sy}"] = float(r[c])
    return q

# ---------------------------------------------------------------- Figure 5 callout pairs (eculizumab inset, panel b)
FIG5_CALLOUTS = {  # key: (disease a, disease b)
    "PNH-NMOSD": ("MONDO:0100244", "MONDO:0019100"),
    "PNH-aHUS": ("MONDO:0100244", "MONDO:0016244"),
    "CESD-FH": ("MONDO:0019149", "MONDO:0005439"),
    "CN2-Gilbert": ("MONDO:0011725", "MONDO:0007745"),
    "CPVT1-SADS": ("MONDO:0010975", "MONDO:0054866"),
    "HD-PD": ("MONDO:0007739", "MONDO:0005180"),
    "CFEOM-ptosis": ("MONDO:0007614", "MONDO:0000728"),
    "GLUT1DS2-IGE": ("MONDO:0012805", "MONDO:0005579"),
}


def q_callouts():
    """Five-seed mean cosine of each Figure 5 callout pair and its percentile among all disease pairs."""
    q = {}
    for t, (a, b) in FIG5_CALLOUTS.items():
        c = pair_cos(a, b); q[f"callout_p|{t}"] = c; q[f"callout_pct|{t}"] = pct_all_pairs(c)
    q["n_all_pairs"] = len(all_pair_cosines())
    return q


def q_prospective():
    """Figure 7 and the event-time Results: conditional-logit odds ratios per 1 SD of score (test pool, popularity
    spline, Efron ties, drug-clustered CIs) and the eculizumab/NMO case (decision year 2007)."""
    import json
    ET = os.path.join(DD, "fig7_prospective")
    A = pd.read_csv(os.path.join(ET, "or_clogit.tsv"), sep="\t", dtype={"set": str})

    def o(model, s="pooled"):
        r = A[(A.set == s) & (A.model == model)]
        assert len(r) == 1, (model, s); return r.iloc[0]
    q = {}
    for k, m in (("dis", "M2"), ("dis_unadj", "M0"), ("perm", "PERM_M2"), ("gene", "M2_gene_jaccard"), ("gobp", "M2_gobp_simgic"),
                 ("hpo", "M2_hpo_simgic"), ("aa", "M2_adamic_adar_w"), ("strength", "M2_strength")):
        r = o(m); q[f"or_{k}"], q[f"or_{k}_lo"], q[f"or_{k}_hi"] = float(r.OR), float(r.lo95), float(r.hi95)
    p = o("M2"); q["or_events"], q["or_drugs"] = int(p.n_cases), int(p.n_drugs)
    Y = A[(A.model == "M2") & A.set.str.isdigit()].assign(y=lambda d: d.set.astype(int))
    Y = Y[(Y.y >= 2000) & (Y.y <= 2021)].set_index("y").sort_index()
    assert list(Y.index) == list(range(2000, 2022)), list(Y.index)
    assert (Y.lo95 > 1).all(), "a yearly odds ratio interval reaches 1"            # 'stays above 1 in each decision year'
    q["or_last_year"] = int(Y.index.max()); q["or_min_year"] = int(Y.OR.idxmin())
    for y in Y.index:
        q[f"or_{y}"], q[f"or_{y}_lo"], q[f"or_{y}_hi"] = float(Y.loc[y].OR), float(Y.loc[y].lo95), float(Y.loc[y].hi95)
    CJ = json.load(open(os.path.join(ET, "nmo_case_clogit.json")))
    q.update(nmo_year=int(CJ["year"]), nmo_rank=int(CJ["disatlas_rank"]), nmo_N=int(CJ["N"]), nmo_top_pct=float(CJ["disatlas_pct_rank"]),
             nmo_pop_rank=int(CJ["popularity_rank"]), nmo_pop_top_pct=float(CJ["popularity_pct_rank"]))
    H = json.load(open(os.path.join(ET, "nmo_hops.json")))["stripped_training_graph"]
    assert H["direct_edge"] is False; q["nmo_hops"] = int(H["hops"])
    CS = pd.read_csv(os.path.join(ET, "nmo_candidate_sets.tsv"), sep="\t")
    hp = CS[CS.pool == "held-out pool"].set_index("set").disatlas_rank
    assert hp["any distance (current)"] == hp["<=2 hops of PNH"] == hp["<=3 hops of PNH"] == q["nmo_rank"]   # 'rank is the same'
    al = CS[(CS.pool == "all embedded diseases") & (CS.set == "any distance (current)")].iloc[0]
    q.update(nmo_all_rank=int(al.disatlas_rank), nmo_all_N=int(al.N), nmo_all_top_pct=100 * al.disatlas_rank / al.N)
    return q


# ---------------------------------------------------------------- architecture
def q_architecture():
    """Embedding dimension of the final (gene-anchored) model, read from its shipped embeddings; pairs are scored by
    the cosine of the embeddings, so there are no further model parameters."""
    dims = {np.load(f)["X"].shape[1] for f in glob.glob(os.path.join(ROOT, "data", "zenodo", "model_outputs", "embeddings", "gene_anchored", "embedding_s*.npz"))}
    assert len(dims) == 1, dims
    return dict(embedding_dim=int(dims.pop()))


def q_txgnn_mapped_and_run_ranges():
    """S9: TxGNN on the pairs its published checkpoint can map, and run-to-run Spearman ranges on the 4,357 pairs.

    Interval: disease-cluster bootstrap over the diseases of the mapped pairs only (own node set, seed 4707),
    the same scheme as the validation block of figures/main/fig6_method_comparison.py.
    """
    from scipy.stats import spearmanr
    R = per_run(); tx = _system_mean(R, "TxGNN published ckpt")[0].values; ok = ~np.isnan(tx)
    Rv = R[ok].reset_index(drop=True); y = Rv.y.values; t = tx[ok]; f = Rv["DisAtlas"].values
    q = dict(txmap_pairs=int(ok.sum()), txmap_auroc_txgnn=auc(y, t), txmap_auroc_full=auc(y, f))
    q["txmap_diff_full_minus_txgnn"] = q["txmap_auroc_full"] - q["txmap_auroc_txgnn"]
    nv = pd.unique(Rv[["a", "b"]].values.ravel()); vix = {d: i for i, d in enumerate(nv)}
    a, b = Rv.a.map(vix).values, Rv.b.map(vix).values; rg = np.random.default_rng(4707); bs = {"txgnn": [], "full": [], "diff": []}
    for _ in range(B_BOOT):
        w = np.bincount(rg.integers(0, len(nv), len(nv)), minlength=len(nv)); ww = (w[a] * w[b]).astype(float); k = ww > 0
        a1, a2 = auc(y[k], t[k], ww[k]), auc(y[k], f[k], ww[k]); bs["txgnn"].append(a1); bs["full"].append(a2); bs["diff"].append(a2 - a1)
    for k, v in bs.items():
        q[f"txmap_{k}_lo"], q[f"txmap_{k}_hi"] = (float(x) for x in np.percentile(v, [2.5, 97.5]))
    # run-to-run agreement on the 4,357 pairs scored by both DisAtlas models; pairs scored by both runs of each pair of runs
    ND = _tsv("fig6_comparison", "node_disjoint_per_seed.tsv").set_index("blind_id").loc[R.blind_id].filter(like="DisAtlas node-disjoint|")
    okn = ND.notna().all(axis=1).values; R4 = R[okn]; N4 = ND[okn]
    groups = {"DisAtlas node-disjoint": N4, "DisAtlas": R4.filter(like="DisAtlas|s")}
    for m in ("Claude-Chat Opus 5.5", "Claude-Sci Opus 5.5", "Gemini Pro 3.1", "Phylo"):
        groups[m] = R4[[c for c in R4.columns if c.startswith(m + "|")]]
    q["ra4357|pairs"] = int(okn.sum())
    for m, X in groups.items():
        cols = list(X.columns)
        rs = [spearmanr(X[[c1, c2]].dropna()[c1], X[[c1, c2]].dropna()[c2])[0] for i, c1 in enumerate(cols) for c2 in cols[i + 1:]]
        q[f"ra4357|{m}|median"], q[f"ra4357|{m}|min"], q[f"ra4357|{m}|max"] = float(np.median(rs)), float(min(rs)), float(max(rs))
        q[f"ra4357|{m}|n_all_runs"] = int(X.notna().all(axis=1).sum())
    return q




def q_release_text():
    """Numbers added in the final revised text: development-pair AUROC, HPO coverage, degree medians,
    TF-IDF graph-SVD baseline difference, and the release / rebuild edge counts."""
    import json
    P = _tsv("fig5_clustering", "ixidn_matched_pair_scores.tsv"); D = P[P.pool == "tune"]
    q = {"dev_pairs": len(D), "auroc_dev": auc(D.y, D.p)}
    O = pd.read_csv(os.path.join(DD, "ontology_baseline", "ontology_baselines.csv"))
    h = O[(O.scope == "all") & (O.column == "hpo_simgic")].iloc[0]
    q.update(hpo_scorable_pairs=float(h.n_scorable), hpo_scorable_pct=100 * float(h.n_scorable) / len(P))
    S = _tsv("degree_distributions", "degree_summary.tsv"); S = S[S.graph == "combined"].set_index("group")
    q.update(median_degree_rare=float(S.loc["rare", "median_degree"]), median_degree_common=float(S.loc["common", "median_degree"]))
    B = _tsv("si_baselines", "baselines_auroc_B2000.tsv")
    b = B[(B.block == "in_graph") & (B.subset == "all") & (B.score == "tfidf_graph_svd_cosine__full_graph")].iloc[0]
    q.update({"tfidf_svd_ingraph|" + c: float(b[c]) for c in ("auc", "lo", "hi", "d", "dlo", "dhi")})
    sig = lambda f: json.load(open(os.path.join(ROOT, "tests", f)))
    q["edges_pfull_release"] = float(sig("pfull_graph_signature.json")["edges"])
    q["edges_pfull_rebuild"] = float(sig("pfull_rebuild_signature.json")["warm"]["edges"])  # rebuild with the shipped broad graph
    return q


def all_quantities():
    q = {}
    for f in (q_disease_set, q_cocluster, q_strata, q_controls, q_callouts, q_run_agreement, q_bootstrap_ci, q_ixidn, q_comparison, q_prospective, q_architecture,
              q_txgnn_mapped_and_run_ranges, q_release_text):
        q.update(f())
    return q


if __name__ == "__main__":
    import json, sys
    q = all_quantities()
    out = sys.argv[1] if len(sys.argv) > 1 else None
    small = {k: v for k, v in q.items() if not k.startswith(("cluster|", "cluster_size|"))}
    if out: json.dump(q, open(out, "w"), indent=0, default=float)
    for k, v in small.items(): print(k, round(v, 4) if isinstance(v, float) else v)
