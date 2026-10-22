"""Step 2 — project the DD network to IxIDN+ (disease–disease co-investigation graph).

Manuscript: Methods, "Indication-based disease network (IxIDN)", layers 2-5: each shared
drug adds 1/deg(drug); drugs above the 99.9th degree percentile do not score; a pair with
>= 2 shared drugs is kept if it is in the top-40 partners of either disease (hub cap 100);
rare-rare pairs (>= 1 shared scoring drug) are ranked among rare partners only (top-100,
cap 300); a pair sharing a drug with rare-disease-registry trials on both diseases is always
kept. Also writes the per-disease co-investigated drug sets (own trials of any size, plus
direct MONDO parents/children) used by the IxIDN- veto.

Inputs: data/01_rebuilt/dd_nodes.parquet, dd_survivors.parquet; mondo.obo.
Outputs (data/01_rebuilt/):
  ixidn_pos_edges.parquet           one row per IxIDN+ edge (disease_1 < disease_2)
  disease_coinvestigated_drugs.parquet  disease_id, drug_id, via (own / hier)
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ixidn_lib as L  # noqa: E402


def offdiag(X):
    X = X.tocoo()
    m = (X.row != X.col) & (X.data != 0)
    X = sp.csr_matrix((X.data[m], (X.row[m], X.col[m])), shape=X.shape)
    X.sum_duplicates()
    return X


def rank_partners(M, N):
    """Per row: partner -> rank (1 = best; ties by more shared drugs, then column index)."""
    ranks = {}
    for i in range(M.shape[0]):
        a, b = M.indptr[i], M.indptr[i + 1]
        cols, sc = M.indices[a:b], M.data[a:b]
        if not len(cols):
            continue
        nn = np.asarray(N[i, cols].todense()).ravel()
        order = np.lexsort((cols, -nn, -sc))
        for r, j in enumerate(cols[order], start=1):
            ranks[(i, int(j))] = r
    return ranks


def backbone(ranks, K, Mx):
    keep = set()
    for (i, j), r in ranks.items():
        if i < j:
            rji = ranks[(j, i)]
            if (r <= K and rji <= Mx) or (rji <= K and r <= Mx):
                keep.add((i, j))
    return keep


def coinvestigation(surv_all, dis_nodes):
    """own: {disease -> drugs in its surviving trials of any size};
    hier: {disease -> drugs of its direct MONDO parents/children in the DD network}."""
    cap = L.KNOBS["coinv_hier_max_children"]
    own = surv_all.groupby("disease_id").drug_id.agg(set).to_dict()
    par = L.mondo_is_a()
    kids = defaultdict(set)
    for c, ps in par.items():
        for p in ps:
            kids[p].add(c)
    present = set(dis_nodes)
    hier = defaultdict(set)
    for n in dis_nodes:
        if not n.startswith("MONDO:"):
            continue
        rel = {p for p in par.get(n, ()) if len(kids.get(p, ())) <= cap}
        if len(kids.get(n, ())) <= cap:
            rel |= kids.get(n, set())
        for r in rel:
            if r != n and r in present and own.get(r):
                hier[n] |= own[r]
    return own, hier


def main():
    k = L.KNOBS
    nodes = L.read_out("dd_nodes.parquet")
    s_all = L.read_out("dd_survivors.parquet")
    dis_nodes = nodes[nodes.kind == "disease"]
    label = dis_nodes.set_index("node_id").label
    rare_of = dis_nodes.set_index("node_id").is_rare.astype(bool)
    deg_of = nodes[nodes.kind == "drug"].set_index("node_id").degree

    # layer 0 trial-size cap, then one record per trial for each disease–drug edge
    s = s_all[s_all.passes_size_cap].drop_duplicates(["disease_id", "drug_id", "trial_key"])
    edges = s.groupby(["drug_id", "disease_id"], sort=True)
    surviving = defaultdict(dict)
    for (drug, d), g in edges:
        surviving[drug][d] = g
    orph = s[s.source_code == L.RARE_SOURCE].groupby("drug_id").disease_id.agg(set).to_dict()

    linking = sorted(c for c, m in surviving.items() if len(m) >= 2)
    degs = [deg_of[c] for c in sorted(surviving)]
    thr = float(np.percentile(degs, k["drug_degree_cap_pct"]))
    capped = {c for c in surviving if deg_of[c] > thr}
    dis_list = sorted({d for c in linking for d in surviving[c]})
    di = {d: i for i, d in enumerate(dis_list)}

    def incidence(pairs):
        pairs = list(pairs)
        return sp.csr_matrix((np.ones(len(pairs)), ([di[d] for d, _ in pairs], [j for _, j in pairs])),
                             shape=(len(dis_list), len(linking)))

    B = incidence((d, j) for j, c in enumerate(linking) if c not in capped for d in surviving[c])
    W = sp.diags([0.0 if c in capped else 1.0 / deg_of[c] for c in linking])
    N = offdiag(B @ B.T)
    S = offdiag(B @ W @ B.T)
    Bo = incidence((d, j) for j, c in enumerate(linking) for d in orph.get(c, ()) if d in di)
    O = offdiag(Bo @ Bo.T)
    cand = offdiag(N >= k["min_shared_drugs"])
    M = offdiag(S.multiply(cand))
    ranks = rank_partners(M, N)
    keep = backbone(ranks, k["top_k"], k["max_partners"])
    rare = np.array([bool(rare_of[d]) for d in dis_list])
    Rc = offdiag(N >= k["rare_rare_min_shared"]).tocoo()
    m = rare[Rc.row] & rare[Rc.col]
    Rr = sp.csr_matrix((np.ones(int(m.sum())), (Rc.row[m], Rc.col[m])), shape=N.shape)
    Mr = offdiag(S.multiply(Rr))
    rr = backbone(rank_partners(Mr, N), k["rare_top_k"], k["rare_max_partners"])
    Oc = sp.triu(O, k=1).tocoo()
    orph_pairs = set(zip(Oc.row.tolist(), Oc.col.tolist()))
    keep |= rr | orph_pairs

    drugs_of = defaultdict(set)
    for c in linking:
        for d in surviving[c]:
            drugs_of[d].add(c)
    rows = []
    for i, j in sorted(keep):
        d1, d2 = dis_list[i], dis_list[j]
        shared = sorted(drugs_of[d1] & drugs_of[d2])
        g1 = pd.concat([surviving[c][d1] for c in shared])
        g2 = pd.concat([surviving[c][d2] for c in shared])
        ph = pd.concat([g1.max_phase, g2.max_phase]).dropna()
        ph = ph[ph > 0]
        rows.append({
            "disease_1": d1, "disease_2": d2, "label_1": label.get(d1), "label_2": label.get(d2),
            "n_drugs": len(shared), "n_shared_drugs": int(N[i, j]), "score": round(float(S[i, j]), 6),
            "rank_1": ranks.get((i, j)), "rank_2": ranks.get((j, i)),
            "rare_registry_supported": (i, j) in orph_pairs, "rare_rare_rule": (i, j) in rr,
            "n_trials_d1": g1.trial_key.nunique(), "n_trials_d2": g2.trial_key.nunique(),
            "max_phase_int": float(ph.max()) if len(ph) else np.nan,
            "sources": json.dumps(sorted(set(g1.source_code) | set(g2.source_code))),
            "shared_drugs": json.dumps(shared)})
    out = pd.DataFrame(rows)
    for c in ("rank_1", "rank_2"):
        out[c] = out[c].astype("Int64")
    out.to_parquet(L.out_path("ixidn_pos_edges.parquet"), index=False)

    own, hier = coinvestigation(s_all, list(dis_nodes.node_id[dis_nodes.is_disease_indication]))
    ixn = sorted(set(out.disease_1) | set(out.disease_2))
    cv = [(n, c, "own") for n in ixn for c in sorted(own.get(n, ()))]
    cv += [(n, c, "hier") for n in ixn for c in sorted(hier.get(n, set()) - own.get(n, set()))]
    pd.DataFrame(cv, columns=["disease_id", "drug_id", "via"]).to_parquet(
        L.out_path("disease_coinvestigated_drugs.parquet"), index=False)
    print(f"IxIDN+: {len(ixn)} diseases, {len(out)} edges "
          f"({len(rr)} rare-rule, {len(orph_pairs)} registry-supported); "
          f"degree cap {thr:.1f}, {len(capped)} capped drugs; {len(cv)} co-investigated disease–drug rows")


if __name__ == "__main__":
    main()
