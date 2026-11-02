"""MONDO terms of the embedded nodes, as used by the compositionality analysis (si_compositionality_recompute.py).

The full ontology file is not shipped; data/derived/si_compositionality/mondo_nodes.tsv (written by
scripts/maintainer/build_mondo_nodes.py) holds the per-node slice the analysis needs.  load_terms() returns the same
{mondo_id: {"name", "is_a"}} and node -> MONDO mapping from either source, so results are identical.
"""
from collections import defaultdict
from pathlib import Path
import json
import pandas as pd

_D = Path(__file__).resolve().parents[2] / "data" / "derived" / "si_compositionality"
TABLE, SOURCE = _D / "mondo_nodes.tsv", _D / "mondo_nodes.source.json"


def parse_obo(path):
    terms, cur = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line == "[Term]": cur = {"is_a": [], "xref": []}; continue
        if line.startswith("[") and line.endswith("]"): cur = None; continue
        if cur is None or ": " not in line: continue
        k, v = line.split(": ", 1)
        if k == "id": cur["id"] = v; terms[v] = cur
        elif k == "name": cur["name"] = v
        elif k == "is_a": cur["is_a"].append(v.split(" ")[0])
        elif k == "xref": cur["xref"].append(v.split(" ")[0])
        elif k == "is_obsolete" and v == "true": cur["obsolete"] = True
    return {k: v for k, v in terms.items() if k.startswith("MONDO:") and not v.get("obsolete")}


def _mondo_of(terms, nodes):
    x2m = defaultdict(set)
    for k, v in terms.items():
        for x in v["xref"]: x2m[x].add(k)

    def to_mondo(d):
        if d.startswith("MONDO:"): return d if d in terms else None
        c = x2m.get(d, set()); return sorted(c)[0] if len(c) == 1 else None
    return {d: to_mondo(d) for d in nodes}


def node_table(terms, nodes):
    m_of, inn = _mondo_of(terms, nodes), set(nodes)
    rows = []
    for d in nodes:
        m = m_of[d]; t = terms.get(m, {}) if m else {}
        rows.append(dict(node_id=d, mondo_id=m, name=t.get("name"),
                         is_a_embedded=";".join(p for p in t.get("is_a", []) if p in inn) if m else None))
    return pd.DataFrame(rows)


def load_terms(nodes, obo=None):
    """-> (terms, mondo_of, source_md5 or None). From the ontology file when given, else from the shipped table."""
    if obo:
        terms = parse_obo(obo); return terms, _mondo_of(terms, nodes), None
    T = pd.read_csv(TABLE, sep="\t", dtype=str, keep_default_na=False)
    assert T.node_id.tolist() == list(nodes), "mondo_nodes.tsv does not match the model's nodes.txt"
    terms, m_of = {}, {}
    for r in T.itertuples():
        m_of[r.node_id] = r.mondo_id or None
        if r.mondo_id:
            t = terms.setdefault(r.mondo_id, {"is_a": [p for p in r.is_a_embedded.split(";") if p]})
            if r.name: t["name"] = r.name
    return terms, m_of, json.load(open(SOURCE))["mondo.obo"]
