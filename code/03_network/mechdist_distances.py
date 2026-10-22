"""MechDist stages S1-S5 from the raw tables in data/00_raw/ (Methods, mechanism distance).

S1 GO distance (gene-set BMA of min-ancestral-path term distance / d_max, BP+0.25 MF+0.1 CC, GO terms
annotated to >2% of disease genes masked at distance time); S2 HPO distance (symptom-set BMA);
S3 cut-off 0.2 on either channel; S4 candidate-pair union + proximal channels px_HPA ((1 - cosine)/2 of
cell-type profiles) and px_DOM (1 - Jaccard of Pfam sets); S5 CRank order and per-class edge budget.
S1/S2 call code/00_ontology/disease_distance.py (same tables as data/00_rebuilt_ontology/).

Reads data/03_rebuilt/frozen/disease_sort_rank.parquet (stored tie order).
Writes data/03_rebuilt/work/{disease_go_distance,disease_hpo_distance,dis_dist_all_df,
dis_dist_all_df_with_crank}.parquet and data/03_rebuilt/{S3_thresholds,S5_px_budget}.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "code" / "00_ontology"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import disease_distance as DD  # noqa: E402
import crank  # noqa: E402

RAW = REPO_ROOT / "data" / "00_raw"
OUT = REPO_ROOT / "data" / "03_rebuilt"
WORK = OUT / "work"

CUT_GO = CUT_HPO = 0.2
PX_CHANNELS = ("HPA", "DOM")
CRANK_PRIOR, CRANK_NUM_BIN = 0.05, 100
CRANK_NAN_POLICY = "order_free"
TOP_N = 200_000


def load_raw():
    dg = pd.read_parquet(RAW / "associations" / "disease_gene.parquet").rename(columns={"gene_id": "gene_symbol"})
    dh = pd.read_parquet(RAW / "associations" / "disease_symptom.parquet")
    g2g = pd.read_parquet(RAW / "genes" / "gene_go.parquet").rename(columns={"gene_id": "gene_symbol", "go_id": "GO_ID"})
    hpa = pd.read_parquet(RAW / "genes" / "gene_hpa.parquet").set_index("gene_id")
    dom = pd.read_parquet(RAW / "genes" / "gene_pfam.parquet")
    return dg, dh, g2g, {"HPA": hpa, "DOM": dom}


def _bma_from_gene_matrix(gg):
    if np.isnan(gg).all():
        return np.nan
    row_ok = ~np.isnan(gg).all(axis=1)
    col_ok = ~np.isnan(gg).all(axis=0)
    row_min = np.nanmin(gg[row_ok], axis=1)
    col_min = np.nanmin(gg[:, col_ok], axis=0)
    return 0.5 * (row_min.mean() + col_min.mean())


def px_gene_matrix(channel, ann, genes):
    if channel == 'HPA':
        row = {g: i for i, g in enumerate(ann.index)}
        Xv = ann.to_numpy()
        has = np.array([g in row for g in genes])
        V = np.zeros((len(genes), Xv.shape[1]))
        V[has] = Xv[[row[g] for g, h in zip(genes, has) if h]]
        D = (1.0 - V @ V.T) / 2
    elif channel == 'DOM':
        import scipy.sparse as sp
        gi = {g: i for i, g in enumerate(genes)}
        sub = ann[ann.gene_id.isin(gi)]
        doms = {d: i for i, d in enumerate(sorted(set(sub.pfam_id)))}
        B = sp.csr_matrix((np.ones(len(sub)), (sub.gene_id.map(gi).to_numpy(), sub.pfam_id.map(doms).to_numpy())), shape=(len(genes), len(doms)))
        inter = (B @ B.T).toarray()
        n = np.asarray(B.sum(1)).ravel()
        with np.errstate(invalid='ignore', divide='ignore'):
            D = 1.0 - inter / (n[:, None] + n[None, :] - inter)
        has = n > 0
    else:
        raise ValueError(f'unknown px channel {channel!r} (HPA | DOM)')
    D[~has] = np.nan
    D[:, ~has] = np.nan
    return D


def px_disease_distances(pairs, dis_gene, channel, ann, gene_col='gene_symbol'):
    dgenes = dis_gene.drop_duplicates(['disease_id', gene_col]).groupby('disease_id')[gene_col].apply(sorted).to_dict()
    genes = sorted(set(dis_gene[gene_col]))
    gidx = {g: i for i, g in enumerate(genes)}
    D = px_gene_matrix(channel, ann, genes)
    gi = {d: np.array([gidx[g] for g in gs], dtype=np.int64) for d, gs in dgenes.items()}
    out = np.full(len(pairs), np.nan)
    for k, (a, b) in enumerate(pairs):
        ia, ib = (gi.get(a), gi.get(b))
        if ia is not None and ib is not None:
            out[k] = _bma_from_gene_matrix(D[np.ix_(ia, ib)])
    return out


def _key(df):
    a = df["disease_a"].astype(str)
    b = df["disease_b"].astype(str)
    lo = np.where(a < b, a, b)
    hi = np.where(a < b, b, a)
    return pd.Series([f"{x}__{y}" for x, y in zip(lo, hi)], index=df.index)


def tie_order(m, pairs):
    """Row order of the candidate frame = production order (pairs sorted by the original disease-ID
    strings). Blinded IDs sort differently and CRank output has exact ties, so the stable crank sort
    needs this order; ``frozen/disease_sort_rank.parquet`` stores each disease's original-string rank."""
    rk = pd.read_parquet(OUT / "frozen" / "disease_sort_rank.parquet").set_index("disease_id")["sort_rank"]
    a = rk.reindex([x for x, _ in pairs]).to_numpy()
    b = rk.reindex([y for _, y in pairs]).to_numpy()
    assert not (np.isnan(a.astype(float)).any() or np.isnan(b.astype(float)).any()), "disease missing from sort rank"
    o = np.lexsort((np.maximum(a, b), np.minimum(a, b)))
    return m.iloc[o].reset_index(drop=True), [pairs[i] for i in o]


def channel_class(df):
    dist_cols = sorted(c for c in df.columns if c.startswith("dist_"))
    avail = {c: df[c].notna().to_numpy() for c in dist_cols}
    n_avail = np.sum(list(avail.values()), axis=0)
    out = np.full(len(df), "none", dtype=object)
    out[n_avail == len(dist_cols)] = "all"
    for c, mask in avail.items():
        out[mask & (n_avail == 1)] = c[len("dist_"):].lower() + "_only"
    mixed = (n_avail > 1) & (n_avail < len(dist_cols))
    if mixed.any():
        out[mixed] = np.array(["+".join(sorted(k[len("dist_"):].lower() for k, m in avail.items() if m[i]))
                               for i in np.flatnonzero(mixed)], dtype=object)
    return out


def s1(dg, g2g):
    dg = dg.rename(columns={"gene_symbol": "gene_id"})
    g2g = g2g.rename(columns={"gene_symbol": "gene_id", "GO_ID": "go_id"})
    go = DD.go_disease_distances(DD.pairs_sharing_entity(dg, "gene_id"), dg, g2g)
    print(f"S1 GO: {len(go):,} pairs")
    return go


def s2(dh):
    hpo = DD.hpo_disease_distances(DD.pairs_sharing_entity(dh, "hpo_id"), dh)
    print(f"S2 HPO: {len(hpo):,} pairs")
    return hpo


def s3(go, hpo):
    g = go["dist_GO"].dropna().to_numpy()
    h = hpo["dist_HPO"].dropna().to_numpy()
    return {"cut_off_dist_GO_BP": CUT_GO, "cut_off_dist_HPO": CUT_HPO,
            "go": {"n": len(g), "mean": float(g.mean()), "frac_le_cut": float((g <= CUT_GO).mean())},
            "hpo": {"n": len(h), "mean": float(h.mean()), "frac_le_cut": float((h <= CUT_HPO).mean())}}


def s4(go, hpo, dg, ann):
    go = go.assign(pair=_key(go))
    hpo = hpo.assign(pair=_key(hpo))
    E1 = go.loc[go["dist_GO"] <= CUT_GO, ["pair", "dist_GO"]].dropna().drop_duplicates("pair")
    E2 = hpo.loc[hpo["dist_HPO"] <= CUT_HPO, ["pair", "dist_HPO"]].dropna().drop_duplicates("pair")
    allp = pd.DataFrame({"pair": pd.Index(E1["pair"]).union(pd.Index(E2["pair"]))})
    m = allp.merge(E2, on="pair", how="left").merge(E1, on="pair", how="left")[["pair", "dist_HPO", "dist_GO"]]
    ab = m["pair"].str.split("__", n=1, expand=True)
    pairs = list(zip(ab[0], ab[1]))
    m, pairs = tie_order(m, pairs)
    for ch in PX_CHANNELS:
        m["px_" + ch] = px_disease_distances(pairs, dg, ch, ann[ch])
    print(f"S4 candidate pairs: {len(m):,} (GO {len(E1):,}, HPO {len(E2):,})")
    return m


def px_budget(m, px):
    sel = m[m.dist_GO.isna() & m.dist_HPO.notna()].head(TOP_N)
    n_px = int(sel[px].notna().any(axis=1).sum())
    tmp = m[["pair", "dist_GO", "dist_HPO"]].assign(dist_PX=m[px].mean(axis=1, skipna=True))
    classes = pd.Series(channel_class(tmp)).value_counts().to_dict()
    quotas = {c: TOP_N for c in classes}
    quotas["hpo+px"], quotas["hpo_only"] = n_px, len(sel) - n_px
    return {"quotas": quotas, "hpo_only_slots": int(len(sel)),
            "class_pool": {k: int(v) for k, v in classes.items()}, "px_channels": list(PX_CHANNELS), "top_N": TOP_N}


def s5(m):
    px = ["px_" + c for c in PX_CHANNELS]
    cols = ["dist_GO", "dist_HPO"] + px
    for c in cols:
        m["1 - " + c] = 1 - m[c]
    guess, _, _ = crank.crank_rra(m[["1 - dist_GO", "1 - dist_HPO"]].values, prior=CRANK_PRIOR,
                                  num_bin=CRANK_NUM_BIN, nan_policy=CRANK_NAN_POLICY)
    budget = px_budget(m.assign(crank=guess).sort_values("crank", kind="mergesort"), px)
    guess, _, _ = crank.crank_rra(m[["1 - " + c for c in cols]].values, prior=CRANK_PRIOR,
                                  num_bin=CRANK_NUM_BIN, nan_policy=CRANK_NAN_POLICY)
    m["crank"] = guess
    m = m.sort_values("crank", kind="mergesort").reset_index(drop=True)
    m["dist_PX"] = m[px].mean(axis=1, skipna=True)
    return m, budget


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    dg, dh, g2g, ann = load_raw()
    go = s1(dg, g2g)
    go.to_parquet(WORK / "disease_go_distance.parquet")
    hpo = s2(dh)
    hpo.to_parquet(WORK / "disease_hpo_distance.parquet")
    json.dump(s3(go, hpo), open(OUT / "S3_thresholds.json", "w"), indent=2)
    m = s4(go, hpo, dg, ann)
    m.to_parquet(WORK / "dis_dist_all_df.parquet")
    m, budget = s5(m)
    json.dump(budget, open(OUT / "S5_px_budget.json", "w"), indent=2)
    m.drop(columns=[c for c in m.columns if c.startswith("1 - ")]).to_parquet(WORK / "dis_dist_all_df_with_crank.parquet")
    print("S5 quotas", budget["quotas"])


if __name__ == "__main__":
    main()
