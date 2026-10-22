"""OBO parsing, term remap, ancestor closure and min-ancestral-path term distances."""

import heapq
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
ONTO_DIR = REPO_ROOT / "data" / "00_raw" / "ontologies"
ZENODO_ONTO_DIR = REPO_ROOT / "data" / "zenodo" / "00_raw" / "ontologies"
HP_OBO = ONTO_DIR / "hp.obo"
GO_OBO = ZENODO_ONTO_DIR / "go-basic.obo"
MONDO_OBO = ZENODO_ONTO_DIR / "mondo.obo"
ANCESTORS_DIR = ONTO_DIR / "ancestors"

GO_ROOTS = {"biological_process": "GO:0008150", "molecular_function": "GO:0003674",
            "cellular_component": "GO:0005575"}
HPO_ROOT = "HP:0000001"
STRUCTURAL_REL = {"is_a", "part_of"}
REG_REL = {"regulates", "positively_regulates", "negatively_regulates"}
# variant key -> (namespace, regulates-family weight); is_a / part_of always weight 1.
GO_VARIANTS = {"go_bp_semi": ("biological_process", 2), "go_mf": ("molecular_function", 2),
               "go_cc": ("cellular_component", 2)}

_TAG = re.compile(r"^(?P<tag>.+?):\s*(?P<value>.+?)\s*(?P<mod>(?<!\\)\{.*?(?<!\\)\})?\s*(?P<c>(?<!\\)!.*?)?$")


def read_obo(path):
    """Parse [Term] stanzas -> ({id: {tag: [values]}}, [(child, relation, parent)])."""
    terms, cur = {}, None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("["):
                if cur and "id" in cur:
                    terms[cur["id"][0]] = cur
                cur = {} if line == "[Term]" else None
            elif cur is not None and line and not line.startswith("!"):
                m = _TAG.match(line)
                if m:
                    cur.setdefault(m.group("tag"), []).append(m.group("value"))
    if cur and "id" in cur:
        terms[cur["id"][0]] = cur
    edges = []
    for tid, t in terms.items():
        edges += [(tid, "is_a", p.split()[0]) for p in t.get("is_a", [])]
        edges += [(tid, *r.split()[:2]) for r in t.get("relationship", [])]
    return terms, edges


def node_table(terms, edges):
    """(all graph nodes incl. edge targets, {node: namespace}, obsolete set)."""
    nodes = set(terms) | {p for _, _, p in edges}
    ns = {n: terms.get(n, {}).get("namespace", [None])[0] for n in nodes}
    obsolete = {n for n in nodes if terms.get(n, {}).get("is_obsolete", ["false"])[0] == "true"}
    return nodes, ns, obsolete


def term_remap(terms, edges, same_namespace):
    """Obsolete/alt id -> active canonical id: alt_id, then single replaced_by, then single consider.

    ``same_namespace`` (GO) requires replaced_by / consider targets to share the source namespace.
    Returns (remap [source_id, canonical_id, reason], unmappable [id, reason]).
    """
    nodes, ns, obsolete = node_table(terms, edges)
    active = nodes - obsolete
    alt_owner = {a: n for n, t in terms.items() for a in t.get("alt_id", [])}
    remap, unmap = [], []
    for src in sorted(obsolete | set(alt_owner)):
        canon = reason = None
        if src in alt_owner and alt_owner[src] in active:
            canon, reason = alt_owner[src], "alt_id"
        multi = None
        for tag, rname in (("replaced_by", "replaced_by"), ("consider", "consider_single")):
            if canon is not None or multi:
                break
            hits = [t for t in terms.get(src, {}).get(tag, [])
                    if t in active and (not same_namespace or ns.get(t) == ns.get(src))]
            if len(hits) == 1:
                canon, reason = hits[0], rname
            elif len(hits) > 1:
                multi = tag + "_multi"
        if multi:
            unmap.append((src, multi))
        elif canon is not None:
            remap.append((src, canon, reason))
        elif src in obsolete:
            unmap.append((src, "no_successor"))
    return (pd.DataFrame(remap, columns=["source_id", "canonical_id", "reason"]),
            pd.DataFrame(unmap, columns=["id", "reason"]))


def parent_graph(terms, edges, namespace=None, reg_weight=None, rels=None):
    """Active-term adjacency {child: {parent: weight}} (smallest weight over parallel relations).

    GO: is_a / part_of weight 1, regulates family ``reg_weight`` (None drops it), one namespace.
    HPO / MONDO: ``rels={'is_a'}`` (weight 1).
    """
    nodes, ns, obsolete = node_table(terms, edges)
    keep = {n for n in nodes - obsolete if namespace is None or ns[n] == namespace}
    up = {n: {} for n in keep}
    for c, rel, p in edges:
        if c not in keep or p not in keep:
            continue
        if rels is not None:
            w = 1 if rel in rels else None
        else:
            w = 1 if rel in STRUCTURAL_REL else (reg_weight if rel in REG_REL else None)
        if w is not None and w < up[c].get(p, np.inf):
            up[c][p] = w
    return up


def _dijkstra(adj, src):
    dist, heap = {src: 0}, [(0, src)]
    while heap:
        d, u = heapq.heappop(heap)
        if d > dist[u]:
            continue
        for v, w in adj.get(u, {}).items():
            if d + w < dist.get(v, np.inf):
                dist[v] = d + w
                heapq.heappush(heap, (d + w, v))
    return dist


def ancestor_closure(up):
    """Long DataFrame [node, ancestor, dist int32]: shortest weighted up-distance, self at 0."""
    rows = [(n, a, d) for n in sorted(up) for a, d in _dijkstra(up, n).items()]
    df = pd.DataFrame(rows, columns=["node", "ancestor", "dist"])
    df["dist"] = df["dist"].astype(np.int32)
    return df


def depth_from_root(up, root):
    """{term: shortest weighted down-distance from ``root``}."""
    down = {}
    for c, ps in up.items():
        for p, w in ps.items():
            down.setdefault(p, {})[c] = w
    return _dijkstra(down, root)


def build_go(path=GO_OBO):
    """{variant key: (sorted nodes, ancestors, remap)} for go_bp_semi / go_mf / go_cc."""
    terms, edges = read_obo(path)
    remap, _ = term_remap(terms, edges, same_namespace=True)
    _, ns, _ = node_table(terms, edges)
    out = {}
    for key, (namespace, w) in GO_VARIANTS.items():
        up = parent_graph(terms, edges, namespace=namespace, reg_weight=w)
        rm = remap[remap["canonical_id"].map(ns) == namespace].reset_index(drop=True)
        out[key] = (sorted(up), ancestor_closure(up), rm)
    return out


def build_hpo(path=HP_OBO):
    """(sorted nodes, ancestors, remap) for HPO, is_a only."""
    terms, edges = read_obo(path)
    remap, _ = term_remap(terms, edges, same_namespace=False)
    up = parent_graph(terms, edges, rels={"is_a"})
    return sorted(up), ancestor_closure(up), remap


def load_ancestors(key, anc_dir=ANCESTORS_DIR):
    """Shipped (nodes, ancestors DataFrame, {source_id: canonical_id}) for go_bp_semi | go_mf | go_cc | hpo."""
    d = Path(anc_dir)
    nodes = pd.read_csv(d / f"{key}_nodes.csv").iloc[:, 0].tolist()
    anc = pd.read_parquet(d / f"{key}_ancestors.parquet")
    rm = pd.read_csv(d / f"{key}_term_remap.csv")
    return nodes, anc, dict(zip(rm["source_id"], rm["canonical_id"]))


def distance_matrix(nodes, anc):
    """Square int16 matrix, dist(x, y) = min over common ancestors c of d(x, c) + d(y, c)."""
    pos = {t: i for i, t in enumerate(nodes)}
    ni = anc["node"].map(pos).to_numpy()
    ai = anc["ancestor"].map(pos).to_numpy()
    di = anc["dist"].to_numpy().astype(np.int16)
    n = len(nodes)
    by_anc = np.argsort(ai, kind="stable")
    a_bounds = np.searchsorted(ai[by_anc], np.arange(n + 1))
    by_node = np.argsort(ni, kind="stable")
    n_bounds = np.searchsorted(ni[by_node], np.arange(n + 1))
    big = np.iinfo(np.int16).max
    M = np.full((n, n), big, dtype=np.int16)
    for x in range(n):
        for k in by_node[n_bounds[x]:n_bounds[x + 1]]:
            sel = by_anc[a_bounds[ai[k]]:a_bounds[ai[k] + 1]]
            np.minimum.at(M[x], ni[sel], di[sel] + di[k])
    assert not (M == big).any(), "unreachable term pair"
    return M


def lca_pair(ax, ay):
    """(min-ancestral-path distance, sorted tied LCAs) for two {ancestor: dist} dicts."""
    common = ax.keys() & ay.keys()
    if not common:
        return None, []
    s = {c: ax[c] + ay[c] for c in common}
    best = min(s.values())
    return best, sorted(c for c, v in s.items() if v == best)
