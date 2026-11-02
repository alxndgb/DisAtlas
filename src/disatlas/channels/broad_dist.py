"""Broad channel distances: four disease-disease distances for every candidate pair.

GO  : disease -> gene -> GO term; gene-gene = min term-term distance / d_max; disease-disease = best-match
      average (BMA) over the two gene sets; BP + 0.25 MF + 0.1 CC over the namespaces present.
HPO : disease -> phenotype term; BMA over the two term sets / d_max.
HPA : gene-gene (1 - cosine) / 2 of cell-type profiles; DOM: gene-gene 1 - Jaccard of Pfam sets; BMA.
Candidates: pairs sharing a gene with GO <= 0.2, or sharing a phenotype with HPO <= 0.2.
"""
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd

NS = {"biological_process": 1.0, "molecular_function": 0.25, "cellular_component": 0.1}
CUT = 0.2
HUB_FRAC = 0.02


def pairs_sharing(df, col):
    inv = defaultdict(set)
    for d, e in zip(df.disease_id, df[col]):
        inv[e].add(d)
    return sorted({p for ds in inv.values() if len(ds) > 1 for p in combinations(sorted(ds), 2)})


def bma(gg):
    """Best-match average; rows/columns that are all NaN are ignored."""
    if np.isnan(gg).all():
        return np.nan
    r, c = ~np.isnan(gg).all(1), ~np.isnan(gg).all(0)
    return 0.5 * (np.nanmin(gg[r], 1).mean() + np.nanmin(gg[:, c], 0).mean())


def gene_terms(go, ns, idx, genes):
    """{gene: sorted term rows}; terms annotated to > HUB_FRAC of the genes are masked unless a gene has no other."""
    s = go[(go.namespace == ns) & go.go_id.isin(idx.keys())].drop_duplicates(["gene", "go_id"])
    frac = s[s.gene.isin(genes)].groupby("go_id").gene.size() / s[s.gene.isin(genes)].gene.nunique()
    hub = set(frac.index[frac > HUB_FRAC])
    out, only_hub = {}, {}
    for g, t in zip(s.gene, s.go_id):
        (only_hub if t in hub else out).setdefault(g, []).append(idx[t])
    for g, v in only_hub.items():
        out.setdefault(g, v)
    return {g: np.array(sorted(set(v)), dtype=np.int64) for g, v in out.items()}


def go_distance(pairs, dg, go, mats):
    """mats: {namespace: (term list, matrix, d_max)} -> dist_GO per pair."""
    genes = dg.groupby("disease_id").target.apply(lambda s: sorted(set(s))).to_dict()
    total = np.zeros(len(pairs))
    seen = np.zeros(len(pairs), dtype=bool)
    for ns, w in NS.items():
        terms, M, dmax = mats[ns]
        rows = gene_terms(go, ns, {t: i for i, t in enumerate(terms)}, set(dg.target))
        cache = {}
        for k, (a, b) in enumerate(pairs):
            g1 = [g for g in genes[a] if g in rows]
            g2 = [g for g in genes[b] if g in rows]
            if not g1 or not g2:
                continue
            gg = np.empty((len(g1), len(g2)))
            for i, x in enumerate(g1):
                for j, y in enumerate(g2):
                    key = (x, y) if x <= y else (y, x)
                    if key not in cache:
                        cache[key] = float(M[np.ix_(rows[x], rows[y])].min())
                    gg[i, j] = cache[key] / dmax
            v = bma(gg)
            if not np.isnan(v):
                total[k] += w * v
                seen[k] = True
    total[~seen] = np.nan
    return total


def hpo_distance(pairs, dp, terms, M, dmax):
    idx = {t: i for i, t in enumerate(terms)}
    rows = defaultdict(set)
    for d, t in zip(dp.disease_id, dp.hpo_id):
        if t in idx:
            rows[d].add(idx[t])
    rows = {d: np.array(sorted(v), dtype=np.int64) for d, v in rows.items()}
    out = np.full(len(pairs), np.nan)
    for k, (a, b) in enumerate(pairs):
        if a in rows and b in rows:
            s = M[np.ix_(rows[a], rows[b])]
            out[k] = 0.5 * (s.min(1).mean() + s.min(0).mean()) / dmax
    return out


def px_distance(pairs, dg, hpa, pfam):
    """(px_HPA, px_DOM) per pair; a gene without the annotation is ignored by the BMA."""
    genes = sorted(set(dg.target))
    gi = {g: i for i, g in enumerate(genes)}
    row = {g: i for i, g in enumerate(hpa.index)}
    X = hpa.to_numpy()
    has = np.array([g in row for g in genes])
    V = np.zeros((len(genes), X.shape[1]))
    V[has] = X[[row[g] for g in genes if g in row]]
    H = (1.0 - V @ V.T) / 2
    H[~has], H[:, ~has] = np.nan, np.nan
    import scipy.sparse as sp
    sub = pfam[pfam.gene.isin(gi)]
    dom = {d: i for i, d in enumerate(sorted(set(sub.pfam_id)))}
    B = sp.csr_matrix((np.ones(len(sub)), (sub.gene.map(gi), sub.pfam_id.map(dom))), shape=(len(genes), len(dom)))
    inter = (B @ B.T).toarray()
    n = np.asarray(B.sum(1)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        P = 1.0 - inter / (n[:, None] + n[None, :] - inter)
    P[n == 0], P[:, n == 0] = np.nan, np.nan
    dgi = {d: np.array([gi[g] for g in sorted(set(s))]) for d, s in dg.groupby("disease_id").target}
    out = np.full((len(pairs), 2), np.nan)
    for k, (a, b) in enumerate(pairs):
        ia, ib = dgi.get(a), dgi.get(b)
        if ia is not None and ib is not None:
            out[k] = bma(H[np.ix_(ia, ib)]), bma(P[np.ix_(ia, ib)])
    return out[:, 0], out[:, 1]


def candidates(dg, dp, go, hpa, pfam, go_mats, hpo_mat):
    """Candidate frame: pair ('a__b', sorted), dist_HPO, dist_GO, px_HPA, px_DOM."""
    gp = pairs_sharing(dg, "target")
    g = pd.DataFrame({"pair": [f"{a}__{b}" for a, b in gp], "dist_GO": go_distance(gp, dg, go, go_mats)})
    hp = pairs_sharing(dp, "hpo_id")
    h = pd.DataFrame({"pair": [f"{a}__{b}" for a, b in hp], "dist_HPO": hpo_distance(hp, dp, *hpo_mat)})
    g, h = g[g.dist_GO <= CUT], h[h.dist_HPO <= CUT]
    m = pd.DataFrame({"pair": pd.Index(g.pair).union(pd.Index(h.pair))})
    m = m.merge(h, on="pair", how="left").merge(g, on="pair", how="left")
    ab = m.pair.str.split("__", n=1, expand=True)
    m["px_HPA"], m["px_DOM"] = px_distance(list(zip(ab[0], ab[1])), dg, hpa, pfam)
    return m

