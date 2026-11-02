"""Rebuild the DisAtlas v3.2 training graph from the shipped post-threshold channels.

Steps (Methods, "MechDist network"):
  1. specific channel: kept disease-module links -> drop modules linked to < 2 diseases ->
     x_dm = sqrt(e_dm / (ln(1 + n_m) L_d)), e_dm = n_fam + ln(1 + n_doc) -> w_spec = X X^T ->
     per-disease top-50 partners, symmetrised -> drop edges leaving the broad disease set.
  2. union with the broad channel; FDA repurposing control pairs removed.
  3. walk weight w = 1[broad] + alpha * f * w_spec / mean(w_spec), alpha = 2, with the mechanism-class
     balancing factor f (classes with less walk mass than the complement class raised by at most 10x).
  4. training graph = union minus every edge joining a full IxIDN+ or IxIDN- pair.

usage: python -m disatlas.graph.build [--data data] [--out data/zenodo/graph] [--check]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os

import numpy as np
import pandas as pd
from scipy import sparse

ALPHA, K_TOP, C_MAX = 2.0, 50, 10.0


def pair_key(a, b) -> pd.Series:
    a = np.asarray(a, dtype=object).astype(str)
    b = np.asarray(b, dtype=object).astype(str)
    lo, hi = np.where(a < b, a, b), np.where(a < b, b, a)
    return pd.Series([f"{x}__{y}" for x, y in zip(lo, hi)], dtype=object)


def specific_channel(links: pd.DataFrame, broad_nodes: set, k: int = K_TOP):
    sp = links.copy()
    sp = sp[sp.groupby("module_id").disease_id.transform("size") >= 2].copy()
    sp["n_m"] = sp.groupby("module_id").disease_id.transform("size")
    sp["e"] = sp.n_fam.astype(float) + np.log1p(sp.docs_max.astype(float))
    dis = sorted(sp.disease_id.unique())
    di = {d: i for i, d in enumerate(dis)}
    mods = sp.drop_duplicates("module_id").sort_values("module_order").module_id.tolist()
    mi = {m: i for i, m in enumerate(mods)}
    Ld = sp.groupby("disease_id").size()
    val = np.sqrt((1.0 / np.log1p(sp.n_m)) * sp.e / sp.disease_id.map(Ld))
    X = sparse.csr_matrix((val.values, (sp.disease_id.map(di).values, sp.module_id.map(mi).values)),
                          shape=(len(dis), len(mods)))
    S = (X @ X.T).tocsr()
    S.setdiag(0)
    S.eliminate_zeros()
    ed = {}
    for i in range(len(dis)):
        a0, a1 = S.indptr[i], S.indptr[i + 1]
        if a1 == a0:
            continue
        d_, j_ = S.data[a0:a1], S.indices[a0:a1]
        sel = np.argsort(-d_, kind="stable")[:k]
        for j, v in zip(j_[sel].tolist(), d_[sel].tolist()):
            ed[(min(i, j), max(i, j))] = v
    spec = pd.DataFrame([(dis[i], dis[j], v) for (i, j), v in ed.items()], columns=["u", "v", "w"])
    stats = dict(links=int(len(sp)), diseases=len(dis), modules=len(mods), raw_edges=int(len(spec)))
    keep = spec.u.isin(broad_nodes) & spec.v.isin(broad_nodes)
    stats["dropped_outside_broad"] = int((~keep).sum())
    spec = spec[keep].reset_index(drop=True)
    stats.update(edges=int(len(spec)), nodes=int(len(set(spec.u) | set(spec.v))))
    return spec, sp, stats


def union_graph(broad: pd.DataFrame, spec: pd.DataFrame, fda_pairs: pd.DataFrame, alpha: float = ALPHA):
    bro = broad[["u", "v"]].assign(in_broad=True)
    un = bro.merge(spec.assign(in_spec=True), on=["u", "v"], how="outer")
    un["in_broad"] = un.in_broad.astype("boolean").fillna(False).astype(bool)
    un["in_spec"] = un.in_spec.astype("boolean").fillna(False).astype(bool)
    stats = dict(union_edges=int(len(un)))
    ck = set(pair_key(fda_pairs.disease_a, fda_pairs.disease_b))
    hit = pair_key(un.u, un.v).isin(ck).values
    stats["fda_pairs"] = len(ck)
    stats["fda_edges_removed"] = int(hit.sum())
    un = un[~hit].reset_index(drop=True)
    un["w_spec"] = un.w.fillna(0.0)
    un = un.drop(columns="w")
    meanw = float(un.loc[un.in_spec, "w_spec"].mean())
    un["w_base"] = un.in_broad.astype(float) + np.where(un.in_spec, alpha * un.w_spec / meanw, 0.0)
    stats.update(edges=int(len(un)), nodes=int(len(set(un.u) | set(un.v))), broad=int(un.in_broad.sum()),
                 specific=int(un.in_spec.sum()), both=int((un.in_broad & un.in_spec).sum()), mean_w_spec=meanw)
    return un, stats


def arm_weights(un: pd.DataFrame, kept_links: pd.DataFrame, classes: pd.DataFrame,
                alpha: float = ALPHA, cmax: float = C_MAX):
    C = classes.assign(module_id="MECHNODE:" + classes.node)
    k = kept_links[kept_links.module_id.str.startswith("MECHNODE:")][["disease_id", "module_id"]].drop_duplicates()
    pr = k.merge(k, on="module_id")
    pr = pr[pr.disease_id_x < pr.disease_id_y].rename(columns={"disease_id_x": "u", "disease_id_y": "v"})
    pr = pr.merge(C[["module_id", "mechanism_class"]], on="module_id")
    ec = pr[["u", "v", "mechanism_class"]].drop_duplicates()
    sp = un[un.in_spec]
    meanw = float(sp.w_spec.mean())
    ec = ec.merge(sp[["u", "v", "w_spec"]], on=["u", "v"])
    ec["mass"] = alpha * ec.w_spec / meanw
    W = ec.groupby("mechanism_class").mass.sum()
    Wstar = float(W["complement"])
    F = pd.DataFrame(dict(W_m=W))
    F["W_star"] = Wstar
    F["factor"] = np.where(F.W_m < Wstar, np.minimum(Wstar / F.W_m, cmax), 1.0)
    F["capped"] = (Wstar / F.W_m) > cmax
    F.loc[F.index.str.startswith("complement"), ["factor", "capped"]] = [1.0, False]
    ec["f"] = ec.mechanism_class.map(F.factor)
    ef = ec.groupby(["u", "v"]).f.max().reset_index()
    out = un.merge(ef, on=["u", "v"], how="left")
    out["f"] = np.where(out.in_spec, out.f.fillna(1.0), 1.0)
    out["w_arm"] = out.in_broad.astype(float) + np.where(out.in_spec, alpha * out.f * out.w_spec / meanw, 0.0)
    assert (out.w_arm > 0).all()
    stats = dict(n_classes=int(len(F)), W_star=Wstar, n_capped=int(F.capped.sum()),
                 edges_raised=int((out.f > 1).sum()))
    return out, F, stats


def strip_benchmark(un: pd.DataFrame, pos: pd.DataFrame, neg: pd.DataFrame):
    key = pair_key(un.u, un.v)
    kp = set(pair_key(pos.disease_a, pos.disease_b))
    kn = set(pair_key(neg.disease_a, neg.disease_b))
    hp, hn = key.isin(kp).values, key.isin(kn).values
    tr = un[~(hp | hn)].reset_index(drop=True)
    nodes = set(un.u) | set(un.v)
    isolated = nodes - (set(tr.u) | set(tr.v))
    stats = dict(ixidn_plus_pairs=len(kp), ixidn_minus_pairs=len(kn), removed=int((hp | hn).sum()),
                 removed_plus=int(hp.sum()), removed_minus=int(hn.sum()), train_edges=int(len(tr)),
                 train_nodes=len(nodes), new_isolates=len(isolated))
    return tr, stats


def edgeset_sha256(u, v) -> str:
    return hashlib.sha256("\n".join(sorted("\t".join(sorted((str(a), str(b)))) for a, b in zip(u, v)))
                          .encode()).hexdigest()


def build(data: str = "data"):
    ch = os.path.join(data, "raw", "channels")
    bm = os.path.join(data, "raw", "benchmarks")
    broad = pd.read_parquet(os.path.join(ch, "broad", "broad_edges.parquet"))
    bnodes = set(pd.read_parquet(os.path.join(ch, "broad", "broad_nodes.parquet")).disease_id)
    links = pd.read_parquet(os.path.join(ch, "specific", "disease_module_links.parquet"))
    classes = pd.read_csv(os.path.join(ch, "specific", "mechanism_node_classes.tsv"), sep="\t")
    fda = pd.read_csv(os.path.join(bm, "fda_repurposing", "fda_repurposing_pairs.tsv"), sep="\t")
    pos = pd.read_parquet(os.path.join(bm, "ixidn", "ixidn_plus_pairs.parquet"))
    neg = pd.read_parquet(os.path.join(bm, "ixidn", "ixidn_minus_pairs.parquet"))
    spec, kept, s_spec = specific_channel(links, bnodes)
    un, s_un = union_graph(broad, spec, fda)
    arm, F, s_arm = arm_weights(un, kept, classes)
    tr, s_tr = strip_benchmark(arm, pos, neg)
    stats = dict(broad=dict(edges=int(len(broad)), nodes=len(bnodes)), specific=s_spec, union=s_un,
                 arm=s_arm, strip=s_tr, train_edgeset_sha256=edgeset_sha256(tr.u, tr.v))
    return arm, tr, F, stats


def check_against(arm: pd.DataFrame, tr: pd.DataFrame, graph_dir: str):
    """Assert rebuilt union and training graphs equal the shipped ones (same edges and flags, weights within 1e-12)."""
    for name, new in (("union_graph_edges.parquet", arm), ("training_graph_edges.parquet", tr)):
        ref = pd.read_parquet(os.path.join(graph_dir, name))
        canon = lambda d: d.sort_values(["u", "v"]).reset_index(drop=True)  # noqa: E731
        r, n = canon(ref), canon(new)
        assert len(r) == len(n), (name, len(r), len(n))
        assert (r.u.values == n.u.values).all() and (r.v.values == n.v.values).all(), name
        for col in ("in_broad", "in_spec"):
            if col in r and col in n:
                assert (r[col].values == n[col].values).all(), (name, col)
        # weights agree to machine precision; floating-point summation order differs across CPU types
        for col in ("w_spec", "f", "w_arm"):
            assert np.allclose(r[col].values, n[col].values, rtol=0, atol=1e-12), (name, col)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default=None, help="write training_graph_edges.parquet etc. here")
    ap.add_argument("--check", action="store_true", help="compare with data/zenodo/graph/training_graph_edges.parquet")
    a = ap.parse_args()
    arm, tr, F, stats = build(a.data)
    cols = ["u", "v", "in_broad", "in_spec", "w_spec", "f", "w_arm"]
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        arm[cols].to_parquet(os.path.join(a.out, "union_graph_edges.parquet"), index=False)
        tr[cols].to_parquet(os.path.join(a.out, "training_graph_edges.parquet"), index=False)
        F.to_csv(os.path.join(a.out, "mechanism_class_balancing.tsv"), sep="\t")
        json.dump(stats, open(os.path.join(a.out, "graph_build_stats.json"), "w"), indent=1)
    if a.check:
        check_against(arm, tr, os.path.join(a.data, "zenodo", "graph"))
        print("CHECK rebuilt graphs match shipped (identical edges and flags; weights within 1e-12)")
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
