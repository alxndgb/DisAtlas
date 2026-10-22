"""MechDist network (Methods, mechanism distance network; Results, MechDist network).

Stages S6-S7 on the S5 frame written by mechdist_distances.py: rank-mean fusion of the GO, HPO,
px_HPA and px_DOM similarities (missing channel = last rank; exact ties keep CRank order),
stratified per-channel-class top-N with the S5 budget, per-disease floor of 3, connectivity repair
(up to 2 bridges per island from ranks (top-N, top-K], escalating to the full pool), largest component.

Run from the repo root:
    python code/03_network/mechdist_distances.py
    python code/03_network/build_mechdist_network.py
Writes data/03_rebuilt/mechdist_network_{edges,nodes}.parquet.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mechdist_distances import OUT, WORK, PX_CHANNELS, channel_class  # noqa: E402

TOP_N = 200_000
TOP_K = 1_500_000
MAX_BRIDGES_M = 2
PER_DISEASE_FLOOR = 3


def _split_pairs(pair_series):
    if len(pair_series) == 0:
        return pd.DataFrame({'x': [], 'y': []})
    xy = pair_series.str.split('__', n=1, expand=True)
    return pd.DataFrame({'x': xy[0].values, 'y': xy[1].values})

def rankmean_fuse(df):
    cols = ['dist_GO', 'dist_HPO', 'px_HPA', 'px_DOM']
    sim = np.column_stack([1 - df[c].to_numpy(dtype=float) for c in cols])
    n = len(df)
    ranks = np.column_stack([np.where(np.isnan(sim[:, i]), n, rankdata(-np.nan_to_num(sim[:, i], nan=-np.inf))) for i in range(sim.shape[1])])
    tieb = ranks.mean(1)
    crank = np.empty(n)
    crank[np.lexsort((tieb, tieb))] = np.arange(1, n + 1)
    return crank

def assemble_topn(df, top_n=TOP_N, per_disease_floor=PER_DISEASE_FLOOR, quotas=None):
    work = df.copy()
    if PX_CHANNELS:
        px_cols = ['px_' + c for c in PX_CHANNELS]
        work['dist_PX'] = work[px_cols].mean(axis=1, skipna=True)
    work['_cls'] = channel_class(work)
    work = work.sort_values('crank', kind='mergesort')
    present = list(work['_cls'].unique())
    quotas = dict(quotas or {})
    if not quotas:
        quotas = {c: top_n for c in present}
    parts = [g.head(quotas.get(c, 0)) for c, g in work.groupby('_cls', sort=True) if quotas.get(c, 0) > 0]
    head = pd.concat(parts, ignore_index=True) if parts else work.head(0)
    floor_pairs = set()
    if per_disease_floor:
        ends = _split_pairs(work['pair'])
        long = pd.concat([pd.DataFrame({'p': work['pair'].to_numpy(), 'c': work['_cls'].to_numpy(), 'n': ends[k].to_numpy()}) for k in ('x', 'y')], ignore_index=True)
        long = long.iloc[np.argsort(np.r_[np.arange(len(work)), np.arange(len(work))], kind='mergesort')]
        floor = set(long.p[long.groupby(['c', 'n']).cumcount() < per_disease_floor])
        floor_pairs = floor - set(head['pair'])
        head = work[work['pair'].isin(floor | set(head['pair']))]
    xy = _split_pairs(head['pair'])
    xy['weight'] = 1.0
    xy['channel_class'] = head['_cls'].to_numpy()
    xy['crank'] = head['crank'].to_numpy()
    G = nx.from_pandas_edgelist(xy, 'x', 'y', edge_attr=['weight', 'channel_class', 'crank'])
    for _, _, d in G.edges(data=True):
        d['origin'] = 'topN'
    if floor_pairs:
        fx = _split_pairs(pd.Series(sorted(floor_pairs)))
        for u, v in zip(fx['x'], fx['y']):
            G.edges[u, v]['origin'] = 'floor'
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    return (G, comps)

def _candidate_bridges(df, top_n, top_k):
    band = df.iloc[top_n:top_k]
    xy = _split_pairs(band['pair'])
    return list(zip(xy['x'].values, xy['y'].values, band['crank'].values))

def repair_connectivity(G, comps, df, top_n=TOP_N, top_k=TOP_K, max_bridges_m=MAX_BRIDGES_M):
    df = df.sort_values('crank', kind='mergesort').reset_index(drop=True)
    lcc = set(comps[0])
    islands = [set(cmp) for cmp in comps[1:]]
    attached = set(lcc)
    island_of = {n: i for i, isl in enumerate(islands) for n in isl}
    cands = _candidate_bridges(df, top_n, top_k)
    from collections import defaultdict
    by_island = defaultdict(list)
    for u, v, cr in cands:
        iu, iv = (island_of.get(u), island_of.get(v))
        if iu is not None and iv is None:
            by_island[iu].append((cr, u, v))
        elif iv is not None and iu is None:
            by_island[iv].append((cr, v, u))
        elif iu is not None and iv is not None and (iu != iv):
            by_island[iu].append((cr, u, v))
            by_island[iv].append((cr, v, u))
    for k in by_island:
        by_island[k].sort(key=lambda t: t[0])
    low_data_nodes = set()
    bridges_added = []
    unbridged = []
    order = sorted(range(len(islands)), key=lambda i: by_island[i][0][0] if by_island[i] else np.inf)
    for i in order:
        isl = islands[i]
        added = 0
        for cr, isl_node, other in by_island.get(i, []):
            if added >= max_bridges_m:
                break
            if other in attached and isl_node not in attached:
                G.add_edge(isl_node, other, origin='repair', bridge_rank=float(cr), weight=1.0, channel_class='repair')
                added += 1
        if added:
            attached |= isl
            low_data_nodes |= isl
            bridges_added.append({'island_size': len(isl), 'bridges': added})
        else:
            unbridged.append(sorted(isl))
    report = {'top_n': top_n, 'top_k': top_k, 'max_bridges_m': max_bridges_m, 'n_islands': len(islands), 'n_islands_bridged': len(bridges_added), 'n_bridges_added': int(sum((b['bridges'] for b in bridges_added))), 'n_low_data_nodes': len(low_data_nodes), 'n_unbridged_islands': len(unbridged), 'unbridged_island_nodes': sum((len(u) for u in unbridged))}
    return (G, low_data_nodes, report)

def assemble_repair(df, top_n=TOP_N, top_k=TOP_K, max_bridges_m=MAX_BRIDGES_M, per_disease_floor=PER_DISEASE_FLOOR, quotas=None):
    G, comps = assemble_topn(df, top_n=top_n, per_disease_floor=per_disease_floor, quotas=quotas)
    print(f'assembled: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, {len(comps)} components')
    G, low, rep = repair_connectivity(G, comps, df, top_n, top_k, max_bridges_m=max_bridges_m)
    if rep['n_unbridged_islands']:
        comps2 = sorted(nx.connected_components(G), key=len, reverse=True)
        G, low2, rep2 = repair_connectivity(G, comps2, df, top_n, len(df), max_bridges_m=max_bridges_m)
        low |= low2
        rep['escalated'] = rep2
        if rep2['n_unbridged_islands']:
            print(f"LOG: {rep2['n_unbridged_islands']} islands ({rep2['unbridged_island_nodes']} nodes) unbridgeable at any rank -- dropped by largest_cc().")
    return (G, low, rep)

def largest_cc(G):
    if G.number_of_nodes() == 0:
        return G.copy()
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    return G.subgraph(comps[0]).copy()

def finalize(G, low_data_nodes):
    for n in G.nodes():
        G.nodes[n]['kind'] = 'disease'
        G.nodes[n]['data_support'] = 'low' if n in low_data_nodes else 'full'
    G_int = nx.convert_node_labels_to_integers(G, label_attribute='ID')
    return G_int


def main():
    df = pd.read_parquet(WORK / "dis_dist_all_df_with_crank.parquet")
    df["crank"] = rankmean_fuse(df)
    df = df.sort_values("crank", kind="mergesort").reset_index(drop=True)
    quotas = json.load(open(OUT / "S5_px_budget.json"))["quotas"]
    G, low, rep = assemble_repair(df, quotas=quotas)
    G = largest_cc(finalize(G, low))
    edges = pd.DataFrame([(min(G.nodes[u]["ID"], G.nodes[v]["ID"]), max(G.nodes[u]["ID"], G.nodes[v]["ID"]),
                           d["channel_class"], d.get("crank", d.get("bridge_rank")), d["origin"], d["weight"])
                          for u, v, d in G.edges(data=True)],
                         columns=["disease_a", "disease_b", "channel_class", "crank", "origin", "weight"])
    edges = edges.sort_values(["disease_a", "disease_b"]).reset_index(drop=True)
    nodes = pd.DataFrame([(d["ID"], d["data_support"]) for _, d in G.nodes(data=True)],
                         columns=["disease_id", "data_support"]).sort_values("disease_id").reset_index(drop=True)
    edges.to_parquet(OUT / "mechdist_network_edges.parquet", index=False)
    nodes.to_parquet(OUT / "mechdist_network_nodes.parquet", index=False)
    json.dump(rep, open(OUT / "S7_repair_report.json", "w"), indent=2, default=str)
    print(f"MechDist: {len(nodes):,} nodes, {len(edges):,} edges, "
          f"{int((nodes.data_support == 'low').sum()):,} low-data nodes")


if __name__ == "__main__":
    main()
