"""Maintainer: stage data/raw/mondo/ from MONDO release 2026-02-06 (mondo.obo of the authors' working tree).

    python scripts/maintainer/stage_mondo_body_system.py /path/to/MONDO_060226/mondo.obo

Writes data/raw/mondo/mondo_body_system_13021.tsv (disease_id, mondo_id = the MONDO term the disease maps to,
categories = '|'-joined MONDO ids of its body-system classes) and data/raw/mondo/body_system_classes.tsv (the 19
direct is_a children of MONDO:7770006 'disease by body system or component', id and name). The parsing and mapping
block below is the authors' code verbatim (r3minor_analyses.py, section MONDO).
"""
import json, sys
from collections import defaultdict
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[2]


def main(OBO):
    ids = json.load(open(REPO / "data/zenodo/model_outputs/cosine_matrix/ids_13021.json"))
    LOG = {}
    log = print
    terms = {}; cur = None
    for line in open(OBO, encoding="utf-8"):
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
    terms = {k: v for k, v in terms.items() if k.startswith("MONDO:") and not v.get("obsolete")}
    xref2mondo = defaultdict(set)
    for k, v in terms.items():
        for x in v["xref"]: xref2mondo[x].add(k)
    ROOT = "MONDO:7770006"   # disease by body system or component (19 classes)
    top = sorted(k for k, v in terms.items() if ROOT in v["is_a"])
    anc_cache = {}
    def ancestors(m):
        if m in anc_cache: return anc_cache[m]
        out = set(); st = list(terms.get(m, {}).get("is_a", []))
        while st:
            p = st.pop()
            if p in out or p not in terms: continue
            out.add(p); st.extend(terms[p]["is_a"])
        anc_cache[m] = out; return out
    def to_mondo(d):
        if d.startswith("MONDO:"): return d if d in terms else None
        c = xref2mondo.get(d, set())
        return sorted(c)[0] if len(c) == 1 else None
    mondo_of = {d: to_mondo(d) for d in ids}
    topset = set(top)
    cats = {}
    for d in ids:
        m = mondo_of[d]
        cats[d] = frozenset((ancestors(m) | {m}) & topset) if m else frozenset()
    LOG["mondo"] = dict(n_terms=len(terms), n_top_categories=len(top), embedded_mapped=int(sum(v is not None for v in mondo_of.values())),
                        embedded_with_category=int(sum(len(v) > 0 for v in cats.values())))
    log("mondo", LOG["mondo"])

    out = REPO / "data/raw/mondo"; out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dict(disease_id=ids, mondo_id=[mondo_of[d] or "" for d in ids],
                      categories=["|".join(sorted(cats[d])) for d in ids])).to_csv(out / "mondo_body_system_13021.tsv", sep="\t", index=False)
    pd.DataFrame(dict(mondo_id=top, name=[terms[c]["name"] for c in top])).to_csv(out / "body_system_classes.tsv", sep="\t", index=False)
    print(LOG["mondo"])


if __name__ == "__main__":
    main(sys.argv[1])
