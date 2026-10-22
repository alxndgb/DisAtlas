"""ORDO and MONDO max-LCA disease distances (appendix: Dependence between ORDO and HPO; Negative-set realism).

Distance definition (both ontologies)::

    dist(x, y) = min over common ancestors c of [ longest(c -> x) + longest(c -> y) ]

i.e. the lowest-common-ancestor distance with longest-path (deepest classification) legs.

ORDO: the classification graph is fixed by ``data/00_raw/ontologies/ordo/ordo_ancestors.parquet``
(longest-path ancestor closure of the Orphanet disorder classes, ``part_of`` + ``subClassOf`` edges with
typing metaclasses removed, plus a pseudo-root ``-1`` above the top categories). A MONDO disease is keyed to
ORDO through ``data/00_raw/diseases/orpha_mondo_crosswalk.parquet``: the first crosswalk row (file order)
whose Orphanet code is an ORDO node wins.

MONDO: longest-path ``is_a`` ancestor closure of ``mondo.obo``; only diseases on the evaluation axis
(``data/00_raw/ontologies/ordo/mondo_lca_axis.csv``) receive a distance.

Unmapped pairs get -1.
"""

import functools
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

import ontology as O

ORDO_DIR = O.ONTO_DIR / "ordo"
CROSSWALK = O.REPO_ROOT / "data" / "00_raw" / "diseases" / "orpha_mondo_crosswalk.parquet"


def lca_distance(ax, ay):
    """min over common ancestors of ax[c] + ay[c]; None if no common ancestor."""
    if len(ax) > len(ay):
        ax, ay = ay, ax
    vals = [d + ay[c] for c, d in ax.items() if c in ay]
    return min(vals) if vals else None


def ordo_ancestors():
    """{orpha code (int): {ancestor: longest-path distance}} for the ORDO disorder nodes."""
    anc = pd.read_parquet(ORDO_DIR / "ordo_ancestors.parquet")
    out = defaultdict(dict)
    for n, a, d in anc.itertuples(index=False):
        out[int(n)][int(a)] = int(d)
    return dict(out)


def mondo_to_ordo(ordo_nodes):
    """{MONDO id: orpha code}, first crosswalk row whose Orphanet code is an ORDO node."""
    xw = pd.read_parquet(CROSSWALK)
    out = {}
    for orpha, mondo in zip(xw.orpha_id, xw.mondo_id):
        code = int(orpha.split(":")[1])
        if code in ordo_nodes:
            out.setdefault(mondo, code)
    return out


def mondo_ancestors(path=O.MONDO_OBO):
    """Callable MONDO id -> {ancestor: longest is_a path length} (memoised)."""
    terms, edges = O.read_obo(path)
    up = O.parent_graph(terms, edges, rels={"is_a"})
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 100000))

    @functools.lru_cache(maxsize=None)
    def anc(n):
        d = {n: 0}
        for p in up.get(n, {}):
            for a, k in anc(p).items():
                if k + 1 > d.get(a, -1):
                    d[a] = k + 1
        return d

    return anc


def mondo_axis():
    return set(pd.read_csv(ORDO_DIR / "mondo_lca_axis.csv").mondo_id)


def pair_distances(a, b):
    """DataFrame [ordo_dist, mondo_dist] (int, -1 = unmapped) for MONDO-keyed disease pairs."""
    oa = ordo_ancestors()
    key = mondo_to_ordo(set(oa))
    ma, axis = mondo_ancestors(), mondo_axis()
    od, md = [], []
    for x, y in zip(a, b):
        ox, oy = key.get(x), key.get(y)
        d = lca_distance(oa[ox], oa[oy]) if ox is not None and oy is not None else None
        od.append(-1 if d is None else d)
        d = lca_distance(ma(x), ma(y)) if x in axis and y in axis else None
        md.append(-1 if d is None else d)
    return pd.DataFrame({"ordo_dist": np.array(od, int), "mondo_dist": np.array(md, int)})


def random_ordo_pairs(rng, n=50000):
    """ORDO distances of random disease pairs (reference distribution, Negative-set realism panel a).

    Draws ``n`` codes twice (with replacement) from the ordered disease list in
    ``ordo/random_pair_diseases.csv`` restricted to ORDO nodes, drops self-pairs. Consumes ``rng``
    exactly as the manuscript run did, so the later bootstrap draws line up.
    """
    oa = ordo_ancestors()
    dis = [int(d) for d in pd.read_csv(ORDO_DIR / "random_pair_diseases.csv").orpha_code if int(d) in oa]
    ra = rng.choice(dis, n)
    rb = rng.choice(dis, n)
    m = ra != rb
    return np.array([lca_distance(oa[x], oa[y]) for x, y in zip(ra[m], rb[m])], float)
