"""Step 3 — IxIDN- negative-candidate pairs.

Manuscript: Results "IxIDN-Matched" and Methods "IxIDN-Matched benchmark": negatives are IxIDN
disease pairs at least 3 steps apart in the IxIDN+ co-investigation graph that share no trial
drug with a known mechanism of action, counting drugs trialled for a direct MONDO parent or
child of either disease (co-investigation veto). Pairs in different components are excluded.
All candidates are kept (no sampling).

Inputs: data/01_rebuilt/ixidn_pos_edges.parquet, disease_coinvestigated_drugs.parquet,
dd_nodes.parquet; data/00_raw/frozen/ixidn_eval_pair_universe_frozen.tsv.gz.
Outputs (data/01_rebuilt/):
  ixidn_neg_edges.parquet — disease_1 < disease_2, label_1, label_2, d_sp (hops in IxIDN+),
      degree_cell ("(low, high)" IxIDN-degree quintile bins), ixidn_degree_1, ixidn_degree_2
  ixidn_eval_pairs.parquet — a, b, class (pos / neg): IxIDN+ plus vetoed IxIDN-
  ixidn_train_strip_pairs.parquet — a, b, class, source: IxIDN+ plus all >= 3-hop pairs before the
      veto (source "ixidn"), plus the remaining pairs of the frozen evaluation-pair universe
      (source "eval_universe"); every listed pair is removed from the MechDist training graph
      (Methods, E_IxIDN held out: 6,245 edges stripped, 325,837 training edges)
"""
import os
import sys

import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ixidn_lib as L  # noqa: E402


def main():
    k = L.KNOBS
    pos = L.read_out("ixidn_pos_edges.parquet")
    G = nx.Graph()
    G.add_edges_from(zip(pos.disease_1, pos.disease_2))
    nodes = sorted(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    D = np.full((len(nodes), len(nodes)), -1, dtype=np.int16)
    for s, dd in nx.all_pairs_shortest_path_length(G):
        for t, d in dd.items():
            D[idx[s], idx[t]] = d
    iu, ju = np.triu_indices(len(nodes), k=1)
    arr = np.array(nodes, dtype=object)
    sp = pd.DataFrame({"u": arr[iu], "v": arr[ju], "d": D[iu, ju]})
    cand = sp[sp.d >= k["neg_d_min"]]                      # unreachable (-1) excluded

    cv = L.read_out("disease_coinvestigated_drugs.parquet")
    S = {n: frozenset(g) for n, g in cv.groupby("disease_id").drug_id}
    empty = frozenset()
    veto = np.fromiter((bool(S.get(u, empty) & S.get(v, empty)) for u, v in zip(cand.u, cand.v)),
                       dtype=bool, count=len(cand))
    sel = cand[~veto].reset_index(drop=True)

    deg = pd.Series(dict(G.degree())).sort_index()
    edges = np.unique(np.quantile(deg.values, np.linspace(0, 1, k["neg_degree_bins"] + 1)))
    b = np.clip(np.searchsorted(edges, deg.values, side="right") - 1, 0, len(edges) - 2)
    bins = dict(zip(deg.index, b.astype(int)))
    label = L.read_out("dd_nodes.parquet").set_index("node_id").label
    out = pd.DataFrame({
        "disease_1": sel.u, "disease_2": sel.v,
        "label_1": sel.u.map(label), "label_2": sel.v.map(label),
        "d_sp": sel.d.astype(int),
        "degree_cell": [str((min(bins[u], bins[v]), max(bins[u], bins[v]))) for u, v in zip(sel.u, sel.v)],
        "ixidn_degree_1": sel.u.map(deg).astype(int), "ixidn_degree_2": sel.v.map(deg).astype(int)})
    out.to_parquet(L.out_path("ixidn_neg_edges.parquet"), index=False)

    # pair sets used downstream: evaluation pairs (IxIDN+ and vetoed IxIDN-) and the training-graph
    # strip (IxIDN+ and every >= 3-hop pair before the veto, removed from MechDist training)
    P = pd.DataFrame({"a": pos.disease_1, "b": pos.disease_2, "class": "pos"})
    ev = pd.concat([P, pd.DataFrame({"a": out.disease_1, "b": out.disease_2, "class": "neg"})], ignore_index=True)
    st = pd.concat([P, pd.DataFrame({"a": cand.u.values, "b": cand.v.values, "class": "neg"})], ignore_index=True)
    st["source"] = "ixidn"
    # plus every pair of the frozen evaluation-pair universe not already listed
    u = pd.read_csv(L.EVAL_UNIVERSE, sep="\t", usecols=["a", "b", "label"])
    have = set(map(frozenset, zip(st.a, st.b)))
    u = u[[frozenset(x) not in have for x in zip(u.a, u.b)]]
    st = pd.concat([st, pd.DataFrame({"a": u.a.values, "b": u.b.values,
                                      "class": np.where(u.label.values == 1, "pos", "neg"),
                                      "source": "eval_universe"})], ignore_index=True)
    ev.to_parquet(L.out_path("ixidn_eval_pairs.parquet"), index=False)
    st.to_parquet(L.out_path("ixidn_train_strip_pairs.parquet"), index=False)
    print(f"IxIDN-: {len(cand)} pairs >= {k['neg_d_min']} hops, {int(veto.sum())} vetoed, "
          f"{len(out)} kept over {pd.concat([out.disease_1, out.disease_2]).nunique()} diseases")


if __name__ == "__main__":
    main()
