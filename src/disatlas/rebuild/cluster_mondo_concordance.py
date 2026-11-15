"""Rebuild the "m3_external" block of data/derived/si_cluster/sifu_cluster_13021_values.json (Figure S16 c:
agreement of the K = 55 clusters with the 19 MONDO body-system classes).

    python -m disatlas.rebuild.cluster_mondo_concordance OUTJSON

Inputs: data/raw/mondo/{mondo_body_system_13021.tsv, body_system_classes.tsv} (MONDO 2026-02-06, staged by
scripts/maintainer/stage_mondo_body_system.py) and the K = 55 labels in
data/zenodo/model_outputs/clusters/disease_clusters_k55.parquet. The body below is the authors' code verbatim
(r3minor_analyses.py, section m3), with its random generator numpy default_rng(20261011).
Test: tests/test_rebuild_cluster_mondo_concordance.py.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[3]


def compute():
    ids = json.load(open(ROOT / "data/zenodo/model_outputs/cosine_matrix/ids_13021.json")); N = len(ids)
    idx = {d: i for i, d in enumerate(ids)}
    K = pd.read_parquet(ROOT / "data/zenodo/model_outputs/clusters/disease_clusters_k55.parquet").set_index("disease_id")
    lab = K.loc[ids, "cluster"].to_numpy(int)
    M = pd.read_csv(ROOT / "data/raw/mondo/mondo_body_system_13021.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert list(M.disease_id) == ids
    cats = {d: frozenset(c.split("|")) if c else frozenset() for d, c in zip(M.disease_id, M.categories)}
    C = pd.read_csv(ROOT / "data/raw/mondo/body_system_classes.tsv", sep="\t")
    terms = {m: {"name": n} for m, n in zip(C.mondo_id, C.name)}
    RNG = np.random.default_rng(20261011); LOG = {}
    log = lambda *x: None
    catlist = sorted({c for v in cats.values() for c in v}); ci = {c: i for i, c in enumerate(catlist)}
    Y = np.zeros((N, len(catlist)), bool)
    for d, v in cats.items():
        for c in v: Y[idx[d], ci[c]] = True
    has = Y.any(1)
    glob_frac = Y[has].mean(0)
    def enrich(labels):
        res = []
        for c in np.unique(labels):
            m = (labels == c) & has; n = m.sum()
            if n == 0: continue
            k = Y[m].sum(0)
            p = hypergeom.sf(k - 1, has.sum(), Y[has].sum(0), n)
            res.append((c, n, k, p))
        return res
    def bh(p):
        p = np.asarray(p); o = np.argsort(p); q = p[o] * len(p) / (np.arange(len(p)) + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out
    def n_enriched(labels, min_k=5, fold=2.0):
        res = enrich(labels); ps = []; keys = []
        for c, n, k, p in res:
            for j in range(len(catlist)):
                if k[j] >= min_k and (k[j] / n) / glob_frac[j] >= fold: ps.append(p[j]); keys.append((c, j))
                else: ps.append(1.0); keys.append((c, j))
        q = bh(ps); good = {c for (c, j), qq in zip(keys, q) if qq < 0.05}
        return len(good), res, dict(zip(keys, q))
    def share_pair_rate(labels, npairs=400000, rng=None):
        rng = rng or RNG; hi = np.where(has)[0]
        i = rng.choice(hi, npairs); j = rng.choice(hi, npairs); keep = i != j; i, j = i[keep], j[keep]
        same = labels[i] == labels[j]; share = (Y[i] & Y[j]).any(1)
        return share[same].mean(), share[~same].mean()
    def within_cluster_pairs(labels, nper=400000, rng=None):
        rng = rng or RNG; hi = np.where(has)[0]; out = []
        lab_h = labels[hi]; order = np.argsort(lab_h); lab_s = lab_h[order]; hi_s = hi[order]
        starts = np.searchsorted(lab_s, np.unique(lab_s)); ends = np.r_[starts[1:], len(lab_s)]
        w = (ends - starts).astype(float); w = w * (w - 1); w /= w.sum()
        g = rng.choice(len(starts), nper, p=w)
        i = hi_s[starts[g] + (rng.random(nper) * (ends[g] - starts[g])).astype(int)]
        j = hi_s[starts[g] + (rng.random(nper) * (ends[g] - starts[g])).astype(int)]
        keep = i != j; return (Y[i[keep]] & Y[j[keep]]).any(1).mean()
    hi = np.where(has)[0]
    ri = RNG.choice(hi, 400000); rj = RNG.choice(hi, 400000); kk = ri != rj
    base_share = float((Y[ri[kk]] & Y[rj[kk]]).any(1).mean())
    obs_within = float(within_cluster_pairs(lab))
    n_enr_obs, res_obs, q_obs = n_enriched(lab)
    null_within, null_enr = [], []
    for r in range(200):
        pl = RNG.permutation(lab)
        null_within.append(within_cluster_pairs(pl, nper=100000))
        if r < 50: null_enr.append(n_enriched(pl)[0])
    LOG["m3_external"] = dict(categories=[terms[c]["name"] for c in catlist], category_counts=[int(x) for x in Y.sum(0)],
        criterion="MONDO body-system classes (direct is_a children of MONDO:7770006 'disease by body system or component'), ancestors via is_a",
        n_categories=len(catlist), diseases_with_category=int(has.sum()),
        pair_share_random=base_share, pair_share_within_cluster=obs_within,
        pair_share_within_cluster_null_mean=float(np.mean(null_within)),
        pair_share_within_cluster_null_95=[float(np.quantile(null_within, .025)), float(np.quantile(null_within, .975))],
        perm_p_within=float((1 + np.sum(np.array(null_within) >= obs_within)) / (1 + len(null_within))),
        clusters_with_enriched_category=n_enr_obs, clusters_with_enriched_category_null_max=int(max(null_enr)),
        clusters_with_enriched_category_null_mean=float(np.mean(null_enr)), n_null_perm_enrichment=len(null_enr), n_null_perm_pairs=len(null_within),
        enrichment_rule="hypergeometric, k>=5 members and >=2-fold over background, BH q<0.05 over all clusters x categories")
    log("m3", LOG["m3_external"])

    # adjusted mutual information with each disease's single most specific top category is ill-defined
    # (multi-label); report it on diseases with exactly one top-level category.
    from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
    one = Y.sum(1) == 1
    ami = adjusted_mutual_info_score(lab[one], Y[one].argmax(1)); ari = adjusted_rand_score(lab[one], Y[one].argmax(1))
    ami_null = [adjusted_mutual_info_score(RNG.permutation(lab[one]), Y[one].argmax(1)) for _ in range(100)]
    LOG["m3_external"].update(single_category_diseases=int(one.sum()), ami=float(ami), ari=float(ari),
                              ami_null_max=float(max(ami_null)), ami_null_mean=float(np.mean(ami_null)))
    return LOG["m3_external"]


if __name__ == "__main__":
    json.dump(dict(m3_external=compute()), open(sys.argv[1], "w"), indent=1)
