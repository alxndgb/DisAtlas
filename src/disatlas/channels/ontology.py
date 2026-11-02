"""Term-term distances from an OBO file: d(x, y) = min over common ancestors c of d(c, x) + d(c, y).

GO: one graph per namespace, is_a/part_of weight 1, regulates family weight 2. HPO: is_a only.
"""
from collections import defaultdict

import numpy as np
import pandas as pd

GO_ROOTS = {"biological_process": "GO:0008150", "molecular_function": "GO:0003674",
            "cellular_component": "GO:0005575"}
WEIGHT = {"is_a": 1, "part_of": 1, "regulates": 2, "positively_regulates": 2, "negatively_regulates": 2}
BIG = np.iinfo(np.int16).max


def read_obo(path):
    """{term: {'ns', 'obsolete', 'parents': [(parent, rel)], 'alt', 'replaced_by', 'consider'}}."""
    terms, t = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith("["):
            t = None
            if line == "[Term]":
                t = {"parents": [], "alt": [], "replaced_by": [], "consider": [], "obsolete": False, "ns": None}
            continue
        if t is None or ": " not in line:
            continue
        k, v = line.split(": ", 1)
        v = v.split(" ! ")[0].strip()
        if k == "id":
            terms[v] = t
        elif k == "namespace":
            t["ns"] = v
        elif k == "is_obsolete":
            t["obsolete"] = v == "true"
        elif k == "is_a":
            t["parents"].append((v.split()[0], "is_a"))
        elif k == "relationship":
            rel, p = v.split()[:2]
            t["parents"].append((p, rel))
        elif k == "alt_id":
            t["alt"].append(v)
        elif k in ("replaced_by", "consider"):
            t[k].append(v)
    return terms


def remap(terms, same_ns=True):
    """{alt/obsolete id: active id}: alt_id owner, else single active replaced_by, else single consider."""
    active = {k for k, v in terms.items() if not v["obsolete"]}
    owner = {a: k for k, v in terms.items() for a in v["alt"]}
    out = {}
    for src in sorted(set(owner) | (set(terms) - active)):
        if owner.get(src) in active:
            out[src] = owner[src]
            continue
        t = terms.get(src)
        if t is None:
            continue
        for key in ("replaced_by", "consider"):
            ok = [x for x in t[key] if x in active and (not same_ns or terms[x]["ns"] == t["ns"])]
            if len(ok) == 1:
                out[src] = ok[0]
                break
            if len(ok) > 1:
                break
    return out


def distance_matrix(terms, ns=None):
    """(sorted active term ids, int16 distance matrix) for one namespace (None = all terms)."""
    nodes = sorted(k for k, v in terms.items() if not v["obsolete"] and (ns is None or v["ns"] == ns))
    pos = {n: i for i, n in enumerate(nodes)}
    par = {n: {} for n in nodes}
    for n in nodes:
        for p, rel in terms[n]["parents"]:
            if p in pos and rel in WEIGHT:
                par[n][p] = min(WEIGHT[rel], par[n].get(p, 99))
    anc = {}

    def ancestors(n):  # weighted shortest up-distances, memoised over the DAG
        if n not in anc:
            d = {n: 0}
            for p, w in par[n].items():
                for a, x in ancestors(p).items():
                    if x + w < d.get(a, 99):
                        d[a] = x + w
            anc[n] = d
        return anc[n]

    import sys
    sys.setrecursionlimit(100000)
    desc = defaultdict(list)
    for n in nodes:
        for a, x in ancestors(n).items():
            desc[pos[a]].append((pos[n], x))
    M = np.full((len(nodes), len(nodes)), BIG, dtype=np.int16)
    for c, lst in desc.items():
        s = np.array([i for i, _ in lst])
        d = np.array([x for _, x in lst], dtype=np.int16)
        ix = np.ix_(s, s)
        M[ix] = np.minimum(M[ix], d[:, None] + d[None, :])
    assert (M < BIG).all()
    return nodes, M


def go_matrices(path):
    terms = read_obo(path)
    return {ns: distance_matrix(terms, ns) for ns in GO_ROOTS}, remap(terms)


def hpo_matrix(path):
    terms = read_obo(path)
    return distance_matrix(terms), remap(terms, same_ns=False)


def to_frame(nodes, M, keep=None):
    """Square DataFrame, optionally restricted to the terms in ``keep``."""
    df = pd.DataFrame(M, index=nodes, columns=nodes)
    if keep is not None:
        k = sorted(set(keep) & set(nodes))
        df = df.loc[k, k]
    return df

