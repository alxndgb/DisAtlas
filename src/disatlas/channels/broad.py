"""Broad channel: fuse four disease-disease distances and select the broad edges.

Input: one row per candidate pair (GO or HPO distance <= 0.2) with dist_GO, dist_HPO, px_HPA, px_DOM
(NaN = channel absent). The gene-based distances ship in data/zenodo/channels_raw/broad/pair_gene_distances.parquet;
dist_HPO ships in pair_hpo_distances.parquet (warm) or is recomputed from disease_phenotype + hp.obo (--cold).
Output: broad_edges.parquet (u, v) on the largest connected component.

usage: python -m disatlas.channels.broad --out DIR [--cold]
"""
import argparse
import os

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

TOP_N, TOP_K, BRIDGES, FLOOR = 200_000, 1_500_000, 2, 3
PRIOR, NUM_BIN, NUM_ITER = 0.05, 100, 10
PX = ["px_HPA", "px_DOM"]


def _smooth(x):
    s, changed = x.copy(), 1
    while changed:
        p = s.copy()
        for i in range(1, len(x) - 1):
            s[i] = np.median(p[i - 1:i + 2])
        s[0] = np.median([p[0], s[1], 3 * s[1] - 2 * s[2]])
        s[-1] = np.median([p[-1], s[-2], 3 * s[-2] - 2 * s[-3]])
        changed = np.sum(s != p)
    return s


def _rank(a):
    """Average ranks of -a, NaN all tied at rank n."""
    r = np.full(a.size, float(a.size))
    ok = ~np.isnan(a)
    r[ok] = stats.rankdata(-a[ok])
    return r


def crank(data):
    """CRank (Zitnik et al. 2018) with NaN in the worst bin and exact ties broken by mean rank."""
    nr, nc = data.shape
    nrp = int(np.floor(nr * PRIOR))
    rd = np.column_stack([_rank(data[:, i]) for i in range(nc)]) / nr
    binned = np.ceil(rd * NUM_BIN).astype(int)
    bf, bd = np.zeros((NUM_BIN, nc)), np.zeros((nr, nc))
    guess, cprev = rd.mean(1), 0
    for _ in range(NUM_ITER):
        if 1 - cprev < 1e-15:
            break
        last = guess.copy()
        cut = np.partition(guess, nrp - 1)[nrp - 1]
        better, tied = guess < cut, guess == cut
        pos = better.astype(float)
        pos[tied] = (nrp - better.sum()) / tied.sum()
        for i in range(nc):
            for b in range(1, NUM_BIN + 1):
                sel = binned[:, i] <= b
                bf[b - 1, i] = np.log((pos[sel].sum() + 1) / ((1 - pos)[sel].sum() + 1) / (PRIOR / (1 - PRIOR)))
            bf[:, i] = np.maximum.accumulate(_smooth(bf[:, i])[::-1])[::-1]
            bd[:, i] = bf[binned[:, i] - 1, i]
        guess = stats.rankdata(-bd.sum(1))
        cprev = stats.pearsonr(guess, last)[0]
    s = stats.rankdata(-bd.sum(1), method="dense").astype(np.int64)
    m = stats.rankdata(rd.mean(1), method="dense").astype(np.int64)
    return stats.rankdata(s * (m.max() + 1) + m)


def channel_class(m):
    has = {"go": m.dist_GO.notna().to_numpy(), "hpo": m.dist_HPO.notna().to_numpy(),
           "px": m[PX].notna().any(axis=1).to_numpy()}
    n = sum(has.values())
    out = np.full(len(m), "none", dtype=object)
    out[n == 3] = "all"
    for k, v in has.items():
        out[v & (n == 1)] = k + "_only"
    for a, b in [("go", "hpo"), ("go", "px"), ("hpo", "px")]:
        out[has[a] & has[b] & (n == 2)] = a + "+" + b
    return out


def fuse(m):
    """Quotas from the GO+HPO CRank, then rank-mean of the four channels (CRank order breaks ties).

    The rank-mean order (not the CRank score) is what the manuscript graph uses; CRank only fixes ties.
    """
    m = m.reset_index(drop=True)
    S = 1 - m[["dist_GO", "dist_HPO"] + PX].to_numpy(float)
    m["cls"] = channel_class(m)
    o = m.assign(c=crank(S[:, :2])).sort_values("c", kind="mergesort")
    sel = o[o.dist_GO.isna() & o.dist_HPO.notna()].head(TOP_N)
    n_px = int(sel[PX].notna().any(axis=1).sum())
    quotas = {c: TOP_N for c in m.cls.unique()}
    quotas["hpo+px"], quotas["hpo_only"] = n_px, len(sel) - n_px
    m = m.assign(c=crank(S)).sort_values("c", kind="mergesort").reset_index(drop=True)
    S = 1 - m[["dist_GO", "dist_HPO"] + PX].to_numpy(float)
    rm = np.column_stack([_rank(S[:, i]) for i in range(4)]).mean(1)
    order = np.argsort(rm, kind="stable")
    return m.iloc[order].reset_index(drop=True), quotas


def select(m, quotas):
    """Class quotas, per-disease floor, island repair, largest component."""
    xy = m.pair.str.split("__", n=1, expand=True)
    m = m.assign(x=xy[0].to_numpy(), y=xy[1].to_numpy())
    keep = set()
    for c, g in m.groupby("cls", sort=True):
        keep |= set(g.index[:quotas.get(c, 0)])
    long = pd.DataFrame({"i": np.repeat(m.index, 2), "c": np.repeat(m.cls.to_numpy(), 2),
                         "n": np.column_stack([m.x, m.y]).ravel()})
    keep |= set(long.i[long.groupby(["c", "n"]).cumcount() < FLOOR])
    head = m.loc[sorted(keep)]
    G = nx.from_edgelist(zip(head.x, head.y))
    for top_k in (TOP_K, len(m)):
        comps = sorted(nx.connected_components(G), key=len, reverse=True)
        if len(comps) == 1:
            break
        if _repair(G, comps, m.iloc[TOP_N:top_k]) == 0:
            break
    big = max(nx.connected_components(G), key=len)
    e = pd.DataFrame([tuple(sorted(p)) for p in G.subgraph(big).edges()], columns=["u", "v"])
    return e.sort_values(["u", "v"]).reset_index(drop=True)


def _repair(G, comps, band):
    """Attach each island by up to BRIDGES best band pairs to the attached graph; returns #unbridged."""
    attached, isl = set(comps[0]), {n: i for i, c in enumerate(comps[1:]) for n in c}
    cand = {}
    for u, v in zip(band.x, band.y):
        iu, iv = isl.get(u), isl.get(v)
        if iu is not None and iu != iv:
            cand.setdefault(iu, []).append((u, v))
        if iv is not None and iu != iv:
            cand.setdefault(iv, []).append((v, u))
    first = {k: j for j, k in enumerate(cand)}  # dict order = best band rank first
    unbridged = 0
    for i in sorted(range(len(comps) - 1), key=lambda i: first.get(i, np.inf)):
        added = 0
        for a, b in cand.get(i, []):
            if added < BRIDGES and b in attached and a not in attached:
                G.add_edge(a, b)
                added += 1
        if added:
            attached |= comps[i + 1]
        else:
            unbridged += 1
    return unbridged


RAW = os.path.join("data", "raw")
ZEN = os.path.join("data", "zenodo", "channels_raw", "broad")
GENE_COLS = ["dist_GO", "px_HPA", "px_DOM"]
CUT = 0.2


def _with_pair(df):
    return df.assign(pair=df.disease_a + "__" + df.disease_b).drop(columns=["disease_a", "disease_b"])


def phenotype_distances(obo, dp):
    """dist_HPO (<= CUT) for every pair of diseases sharing a phenotype term, from hp.obo + disease_phenotype."""
    from disatlas.channels import broad_dist, ontology
    (nodes, M), _ = ontology.hpo_matrix(obo)
    pairs = broad_dist.pairs_sharing(dp, "hpo_id")
    d = broad_dist.hpo_distance(pairs, dp, nodes, M, float(M.max()))
    h = pd.DataFrame({"pair": [f"{a}__{b}" for a, b in pairs], "dist_HPO": d})
    return h[h.dist_HPO <= CUT]


def candidate_frame(cold=False, raw=RAW, zen=ZEN):
    """Candidate pairs with the four distances; gene-based ones are shipped, dist_HPO shipped (warm) or recomputed (cold)."""
    m = _with_pair(pd.read_parquet(os.path.join(zen, "pair_gene_distances.parquet")))
    if cold:
        dp = pd.read_parquet(os.path.join(raw, "channels", "broad", "disease_phenotype.parquet"), columns=["disease_id", "hpo_id"])
        h = phenotype_distances(os.path.join(raw, "ontology", "hp.obo"), dp)
    else:
        h = _with_pair(pd.read_parquet(os.path.join(zen, "pair_hpo_distances.parquet")))
    m = m.merge(h, on="pair", how="left", validate="1:1")  # keeps the shipped (sorted) pair order
    return m[["pair", "dist_HPO", "dist_GO"] + PX]


def build(cold=False, raw=RAW, zen=ZEN):
    m, quotas = fuse(candidate_frame(cold, raw, zen))
    return select(m, quotas), quotas


def main():
    ap = argparse.ArgumentParser(description="Rebuild the broad channel (broad_edges.parquet).")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cold", action="store_true", help="recompute dist_HPO from disease_phenotype + hp.obo")
    a = ap.parse_args()
    e, q = build(a.cold)
    os.makedirs(a.out, exist_ok=True)
    e.to_parquet(os.path.join(a.out, "broad_edges.parquet"), index=False)
    print({"nodes": len(set(e.u) | set(e.v)), "edges": len(e), "quotas": q})


if __name__ == "__main__":
    main()
