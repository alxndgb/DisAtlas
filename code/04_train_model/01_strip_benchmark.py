"""Strip benchmark (IxIDN+/IxIDN-) disease-pairs out of the MechDist training graph.

Any pair that is an IxIDN+ positive, an IxIDN- negative candidate, or a pair of the IxIDN
evaluation universe (data/00_raw/frozen/ixidn_eval_pair_universe_frozen.tsv.gz) AND exists as an
edge of the MechDist network is removed before embedding/decoder training, so no benchmark label
can leak into training. With the default inputs this gives the frozen model's training graph
(325,837 edges; FROZEN_MECHDIST_MODEL.json expected_train_edgeset_sha256_blinded_ids).

Isolated nodes created by the strip are KEPT (reported, not dropped) -- node2vec simply cannot
embed a degree-0 node in the next stage.
"""
import argparse
import json
import pickle

import pandas as pd
import networkx as nx


def load_mechdist_graph(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)


def benchmark_pairs_as_ids(edges_df, id_col_1="disease_1", id_col_2="disease_2"):
    return {frozenset((row[id_col_1], row[id_col_2])) for _, row in edges_df.iterrows()}


def strip_for_training(G, pos_pairs, neg_pairs, eval_pairs=frozenset()):
    H = G.copy()
    graph_edges = {frozenset(e) for e in H.edges()}
    present_pos = pos_pairs & graph_edges
    present_neg = neg_pairs & graph_edges
    present_eval = eval_pairs & graph_edges
    to_strip = present_pos | present_neg | present_eval
    for pair in to_strip:
        u, v = tuple(pair)
        if H.has_edge(u, v):
            H.remove_edge(u, v)
    new_isolates = [n for n in H.nodes() if H.degree(n) == 0]
    report = {
        "edges_before": G.number_of_edges(), "edges_after": H.number_of_edges(),
        "n_present_ixidn_pos": len(present_pos), "n_present_ixidn_neg": len(present_neg),
        "n_present_eval_universe": len(present_eval),
        "n_stripped_total": len(to_strip), "n_new_isolates": len(new_isolates),
        "n_components_after": nx.number_connected_components(H),
    }
    return H, report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", default="../../data/zenodo/network/mechdist_graph_v2.1.0.pickle")
    ap.add_argument("--ixidn-pos", default="../../data/02S_build_IxIDN/ixidn_pos_edges.parquet.gz")
    ap.add_argument("--ixidn-neg", default="../../data/02S_build_IxIDN_neg/ixidn_neg_pinned_edges.parquet.gz")
    ap.add_argument("--eval-universe", default="../../data/00_raw/frozen/ixidn_eval_pair_universe_frozen.tsv.gz")
    ap.add_argument("--output", default="../../data/zenodo/network/G_train.pickle")
    args = ap.parse_args()

    G = load_mechdist_graph(args.graph)
    # graph nodes are relabeled ints with the original disease id stored under node attr 'ID'
    id_of = nx.get_node_attributes(G, "ID")
    G = nx.relabel_nodes(G, id_of)

    pos_df = pd.read_parquet(args.ixidn_pos)
    neg_df = pd.read_parquet(args.ixidn_neg)
    pos_pairs = benchmark_pairs_as_ids(pos_df)
    neg_pairs = benchmark_pairs_as_ids(neg_df)

    ev = pd.read_csv(args.eval_universe, sep="\t", usecols=["a", "b"])
    eval_pairs = {frozenset(p) for p in zip(ev.a, ev.b)}

    G_train, report = strip_for_training(G, pos_pairs, neg_pairs, eval_pairs)
    print(json.dumps(report, indent=1))

    with open(args.output, "wb") as fh:
        pickle.dump(G_train, fh)
    print("wrote", args.output, "nodes:", G_train.number_of_nodes(), "edges:", G_train.number_of_edges())


if __name__ == "__main__":
    main()
