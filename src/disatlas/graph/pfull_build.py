"""Build the P_full training graph: Habcd broad channel + released specific channel (disease-disease edges only).

Inputs (repository paths by default)
  specific  data/raw/channels/specific/specific_edges.parquet: u, v, w (u < v). w is the specific-channel term of the
            walk weight after mechanism-class balancing (2 f w_spec / mean(w_spec)); only the edges are released, not
            the evidence, modules, classes or factors behind them.
  broad     the Habcd broad graph (u, v): data/zenodo/graph/graph_edges.parquet, or the output of
            `python -m disatlas.graph.raw_build` (cold).
  tier1     the tier-1 rows of the released disease-gene table (disatlas.graph.raw_build.gene_edges: public sources plus
            the 8,343 released edges), used by the gene-anchored step.
  fda       data/raw/benchmarks/fda_repurposing/fda_repurposing_pairs.tsv: pairs removed from the graph.

Steps  1 outer union of broad and specific edges (lexicographic key order); 2 FDA strip; 3 w_arm = 1[broad] + w;
       4 gene-anchored step: pairs of diseases that share a tier-1 gene linked to <= 28 diseases, with at least one
       endpoint having <= 3 tier-1 genes, and that are not already specific, are added or re-weighted with
       w_arm = 1[broad] + C * W_MED (C = 2 / mean(w_spec), W_MED = median w_spec of the specific channel).
Output u, v, in_broad, in_spec, w_arm (row order = training order; `python -m disatlas.retrain` trains on w_arm and
       removes the IxIDN benchmark pairs itself) and the build statistics.

    python -m disatlas.graph.pfull_build --out build/graph            # broad graph from data/zenodo/graph
    python -m disatlas.graph.pfull_build --broad build/rawgraph/graph_edges.parquet --out build/graph   # cold broad graph
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os

import numpy as np
import pandas as pd

from disatlas import paths

C = float.fromhex("0x1.5444d7455994fp+3")      # 10.633403430413322
W_MED = float.fromhex("0x1.40a20738818afp-3")  # 0.1565590442771216
GENE_MAX_DIS, END_MAX_GENES = 28, 3
SPEC = paths.RAW / "channels" / "specific" / "specific_edges.parquet"
BROAD = paths.ZENODO / "graph" / "graph_edges.parquet"
FDA = paths.RAW / "benchmarks" / "fda_repurposing" / "fda_repurposing_pairs.tsv"
# Graph of record: training graph of the released final model (tests/pfull_graph_signature.json). Its tier-1 table
# includes licence-restricted disease-gene pairs that are not redistributed, so it is documented, not rebuilt.
RECORD = dict(broad_edgeset="85e8163d97a9f7ad065bb5a7caf3b348f9083040eaf5c7faa4b5751b8b0db8ed",
              tier1_set="980916169147b568f77145553439fe56f5bdc74975a9b57ee92196c2c6bc2582",
              out_edgeset="f6030d6d5356f36d70a4540edf1e109eee5c88e796f3274096ba97bd3a5c80b5")
# Rebuild from the released inputs (tests/pfull_rebuild_signature.json): output edge set per broad graph
# (shipped broad graph, or the broad graph rebuilt by disatlas.graph.raw_build).
EXPECT = dict(tier1_set="6abf5961f85732a5a423884d1a20354d5c2abfe675c067785e665ef054c0bd1a",
              out_edgeset={"85e8163d97a9f7ad065bb5a7caf3b348f9083040eaf5c7faa4b5751b8b0db8ed": "60364203caf3fd95a705f8fdc41b895e15de5dd02d69e1bb39bdbff069afdcb0",
                           "8be8ab613feabbd929bab1d2443673b0a4454bbd9a19af84dabbf6c381328daa": "1431e8489ed43512569864f82d2d86c9082a9bf99f8cdfe6b8d2186c4bc1e81a"})
COLS = ["u", "v", "in_broad", "in_spec", "w_arm"]


def pair_key(a, b):
    a, b = np.asarray(a, dtype=str), np.asarray(b, dtype=str)
    return np.where(a < b, np.char.add(np.char.add(a, "__"), b), np.char.add(np.char.add(b, "__"), a))


def edgeset_sha256(u, v):
    """Order-free signature of an undirected edge set."""
    return hashlib.sha256("\n".join(sorted("\t".join(sorted((str(a), str(b)))) for a, b in zip(u, v))).encode()).hexdigest()


def w_arm_sha256(g):
    """Signature of the walk weights in row order (float64 bytes): equal only for a bit-identical w_arm."""
    return hashlib.sha256(g.w_arm.to_numpy(np.float64).tobytes()).hexdigest()


def load_broad(p):
    b = pd.read_parquet(p, columns=["u", "v"]).astype(str)
    return pd.DataFrame(sorted({tuple(sorted(x)) for x in zip(b.u, b.v)}), columns=["u", "v"])


def tier1_table(path=None):
    if path is None:
        from disatlas.graph.raw_build import gene_edges
        g = gene_edges()
    else:
        g = pd.read_parquet(path, columns=["disease_id", "gene_symbol", "evidence_tier"])
    return g[g.evidence_tier == 1][["disease_id", "gene_symbol"]].astype(str).drop_duplicates()


def set_sha(t):
    return hashlib.sha256("\n".join(sorted(t.disease_id + "\t" + t.gene_symbol)).encode()).hexdigest()


def build(spec=SPEC, broad=BROAD, tier1=None, fda=FDA, check=True):
    S = pd.read_parquet(spec)
    assert list(S.columns) == ["u", "v", "w"] and S.w.dtype == np.float64
    assert (S.u < S.v).all() and not S.duplicated(["u", "v"]).any() and (S.w > 0).all()
    B = load_broad(broad)
    T = tier1_table(tier1)
    st = dict(inputs=dict(spec_edges=len(S), broad_edges=len(B), broad_edgeset_sha256=edgeset_sha256(B.u, B.v),
                          tier1_rows=len(T), tier1_set_sha256=set_sha(T)))
    if check:
        assert st["inputs"]["broad_edgeset_sha256"] in EXPECT["out_edgeset"], "broad graph is neither the shipped nor the rebuilt Habcd broad graph"
        assert st["inputs"]["tier1_set_sha256"] == EXPECT["tier1_set"], "tier-1 disease-gene set differs from the Habcd table"
    un = B.assign(in_broad=True).merge(S.assign(in_spec=True), on=["u", "v"], how="outer")
    for c in ("in_broad", "in_spec"):
        un[c] = un[c].astype("boolean").fillna(False).astype(bool)
    F = pd.read_csv(fda, sep="\t", usecols=["disease_a", "disease_b"])
    hit = np.isin(pair_key(un.u, un.v), pair_key(F.disease_a, F.disease_b))
    st.update(union_edges=len(un), fda_edges_removed=int(hit.sum()), fda_spec_edges_removed=int((hit & un.in_spec.values).sum()))
    un = un[~hit].reset_index(drop=True)
    un["w_arm"] = un.in_broad.astype(float) + np.where(un.in_spec, un.w.fillna(0.0), 0.0)
    G0 = un[COLS].copy()
    G0["pair"] = pair_key(G0.u, G0.v)
    N = set(G0.u) | set(G0.v)
    g = T[T.disease_id.isin(N)]
    gd, ng = g.groupby("gene_symbol").disease_id.nunique(), g.groupby("disease_id").gene_symbol.nunique()
    m = g.merge(g, on="gene_symbol")
    m = m[m.disease_id_x < m.disease_id_y]
    m = m[(m.gene_symbol.map(gd) <= GENE_MAX_DIS)
          & ((m.disease_id_x.map(ng) <= END_MAX_GENES) | (m.disease_id_y.map(ng) <= END_MAX_GENES))]
    Q = sorted(set(pair_key(m.disease_id_x, m.disease_id_y)) - set(G0.pair[G0.in_spec]))
    G2 = G0.set_index("pair")
    ex = [p for p in Q if p in G2.index]
    nw = [p for p in Q if p not in G2.index]
    ab = pd.Series(nw, dtype=str).str.split("__", n=1, expand=True)
    A2 = pd.DataFrame(dict(u=ab[0].values, v=ab[1].values, in_broad=False, in_spec=True, w_arm=0.0),
                      index=pd.Index(nw, name="pair"))
    G2.loc[ex, "in_spec"] = True
    G2 = pd.concat([G2, A2])
    sel = G2.index.isin(Q)
    for c in ("in_spec", "in_broad"):
        G2[c] = G2[c].astype(bool)
    G2.loc[sel, "w_arm"] = G2.loc[sel, "in_broad"].astype(float) + C * 1.0 * W_MED
    G2 = G2.reset_index(drop=True)[COLS]
    assert (G2.w_arm > 0).all() and not G2.duplicated(["u", "v"]).any()
    st.update(edges=len(G2), nodes=len(set(G2.u) | set(G2.v)), broad=int(G2.in_broad.sum()), spec=int(G2.in_spec.sum()),
              both=int((G2.in_broad & G2.in_spec).sum()), ga_added=len(nw), ga_reweighted=len(ex),
              edgeset_sha256=edgeset_sha256(G2.u, G2.v), w_arm_sha256=w_arm_sha256(G2))
    if check:
        assert st["edgeset_sha256"] == EXPECT["out_edgeset"][st["inputs"]["broad_edgeset_sha256"]], "P_full edge set differs from the release rebuild"
    return G2, st


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--spec", default=str(SPEC))
    ap.add_argument("--broad", default=str(BROAD))
    ap.add_argument("--tier1", default=None, help="parquet with disease_id, gene_symbol, evidence_tier (default: rebuilt from data/raw)")
    ap.add_argument("--fda", default=str(FDA))
    ap.add_argument("--out", required=True, help="output folder (pfull_graph_edges.parquet, pfull_graph_stats.json)")
    ap.add_argument("--no-check", action="store_true", help="skip the input/output signature checks (ablations)")
    a = ap.parse_args()
    G, st = build(a.spec, a.broad, a.tier1, a.fda, check=not a.no_check)
    os.makedirs(a.out, exist_ok=True)
    G.to_parquet(os.path.join(a.out, "pfull_graph_edges.parquet"), index=False)
    json.dump(st, open(os.path.join(a.out, "pfull_graph_stats.json"), "w"), indent=1)
    print(json.dumps(st))


if __name__ == "__main__":
    main()

