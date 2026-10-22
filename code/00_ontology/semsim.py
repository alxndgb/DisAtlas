"""HPO semantic-similarity baselines (Resnik / Lin BMA, SimGIC) and the direct MechDist HPO distance.

IC(t) = -log p(t), p = share of annotated diseases whose HPO closure contains t (full corpus).
Resnik = IC of the most informative common ancestor; Lin = 2 Resnik / (IC1 + IC2);
set level = symmetric best-match average; SimGIC = IC-weighted Jaccard of the closures.
All matrices are float32, NaN where a disease has no annotation.
"""

import numpy as np


class Onto:
    """Ancestor-or-self closure over an ancestor table [node, ancestor, dist]."""

    def __init__(self, anc, remap=None):
        self.remap = remap or {}
        self.terms = sorted(set(anc["node"]) | set(anc["ancestor"]))
        self.idx = {t: i for i, t in enumerate(self.terms)}
        n = anc["node"].map(self.idx).to_numpy()
        a = anc["ancestor"].map(self.idx).to_numpy()
        lists = [[] for _ in self.terms]
        for i, j in zip(n, a):
            lists[i].append(j)
        self.anc = [np.unique(np.r_[v, i]).astype(np.int64) for i, v in enumerate(lists)]

    def rows(self, ids):
        out = [self.idx.get(self.remap.get(t, t)) for t in ids]
        return np.unique(np.array([i for i in out if i is not None], np.int64))

    def closure(self, rows):
        return np.unique(np.concatenate([self.anc[i] for i in rows])) if len(rows) else np.empty(0, np.int64)

    def corpus_ic(self, entity_rows):
        cnt = np.zeros(len(self.terms))
        for rows in entity_rows.values():
            cnt[self.closure(rows)] += 1
        ic = np.zeros(len(self.terms))
        ok = cnt > 0
        ic[ok] = -np.log(cnt[ok] / len(entity_rows))
        return ic


def resnik_matrix(onto, ic, U):
    """(|U| x |U|) float32 IC of the most informative common ancestor of term rows U."""
    members = {}
    for pos, u in enumerate(U):
        for k in onto.anc[u]:
            members.setdefault(int(k), []).append(pos)
    R = np.zeros((len(U), len(U)), dtype=np.float32)
    for k in sorted(members, key=lambda k: ic[k]):
        if ic[k] > 0:
            m = np.asarray(members[k])
            R[np.ix_(m, m)] = np.float32(ic[k])
    return R


def lin_from_resnik(R, icU):
    s = icU[:, None] + icU[None, :]
    return np.where(s > 0, 2 * R / np.maximum(s, np.float32(1e-12)), np.float32(0)).astype(np.float32)


def bma(S, ent_rows):
    """(n x n) float32 symmetric best-match average of item similarity S over entity item sets."""
    n = len(ent_rows)
    has = np.array([len(r) > 0 for r in ent_rows])
    D = np.full((n, n), np.nan, dtype=np.float32)
    owner = np.concatenate([np.full(len(r), e) for e, r in enumerate(ent_rows)])
    items = np.concatenate(ent_rows)
    order = np.argsort(owner, kind="stable")
    starts = np.searchsorted(owner[order], np.flatnonzero(has))
    for e in np.flatnonzero(has):
        best = np.maximum.reduceat(S[ent_rows[e]][:, items[order]], starts, axis=1)   # (|e| x n_has)
        D[e, has] = best.mean(axis=0, dtype=np.float32)  # D[a, b] = mean over a of best match in b
    return (np.float32(0.5) * (D + D.T)).astype(np.float32)


def simgic(onto, ic, ent_rows):
    cl = [onto.closure(r) for r in ent_rows]
    allc = np.unique(np.concatenate([c for c in cl if len(c)]))
    pos = {t: i for i, t in enumerate(allc)}
    C = np.zeros((len(ent_rows), len(allc)), dtype=np.float32)
    for e, c in enumerate(cl):
        C[e, [pos[t] for t in c]] = 1
    w = ic[allc].astype(np.float32)
    inter = (C * w) @ C.T
    s = C @ w
    G = inter / np.maximum(s[:, None] + s[None, :] - inter, np.float32(1e-12))
    empty = np.array([len(r) == 0 for r in ent_rows])
    G[empty] = np.nan
    G[:, empty] = np.nan
    return G.astype(np.float32)


def hpo_scores(ids, dis_symptom, anc, remap, dist_matrix, matrix_idx, d_max):
    """{method: (n x n) float32} for diseases ``ids``: hpo_resnik_bma, hpo_lin_bma, hpo_simgic, mechdist_hpo."""
    H = Onto(anc, remap)
    corpus = {d: H.rows(v) for d, v in dis_symptom.groupby("disease_id")["hpo_id"].apply(list).items()}
    ic = H.corpus_ic(corpus)
    et = [corpus.get(d, np.empty(0, np.int64)) for d in ids]
    U = np.unique(np.concatenate([r for r in et if len(r)]))
    p = {t: i for i, t in enumerate(U)}
    er = [np.array([p[t] for t in r], np.int64) for r in et]
    R = resnik_matrix(H, ic, U)
    out = {"hpo_resnik_bma": bma(R, er),
           "hpo_lin_bma": bma(lin_from_resnik(R, ic[U].astype(np.float32)), er),
           "hpo_simgic": simgic(H, ic, et)}
    mrows = [np.array([matrix_idx[H.terms[t]] for t in r if H.terms[t] in matrix_idx], np.int64) for r in et]
    Um = np.unique(np.concatenate([r for r in mrows if len(r)]))
    q = {t: i for i, t in enumerate(Um)}
    S = -(dist_matrix[np.ix_(Um, Um)].astype(np.float32) / np.float32(d_max))
    out["mechdist_hpo"] = bma(S, [np.array([q[t] for t in r], np.int64) for r in mrows])
    for v in out.values():
        np.fill_diagonal(v, np.nan)
    return out
