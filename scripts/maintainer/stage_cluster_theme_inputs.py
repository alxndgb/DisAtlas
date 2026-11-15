"""Maintainer: stage the public ontology inputs of the SI cluster-themes rebuild (Supplementary Table S10).

    python scripts/maintainer/stage_cluster_theme_inputs.py MONDO_OBO GO_OBO

Writes
  data/raw/mondo/mondo_terms_13021.tsv.gz   (disease_id, term_id): is_a ancestors-or-self of the MONDO term each of the
                                            13,021 clustered diseases maps to (data/raw/mondo/mondo_body_system_13021.tsv),
                                            restricted to terms with 10-2000 mapped diseases plus 'cancer or benign tumor';
  data/raw/mondo/mondo_terms.tsv            (term_id, name, under_mech_root) for those terms; under_mech_root marks the three
                                            mechanism roots (disease by molecular mechanism, inborn errors of metabolism,
                                            inborn error of immunity) and their descendants;
  data/raw/go/go_bp_terms.tsv.gz            (term_id, name, namespace, ancestors, organ_level): the Gene Ontology terms
                                            annotated in data/raw/specific_channels/gene2go.parquet and their is_a/part_of
                                            ancestors; ancestors '|'-joined; organ_level = the term or an ancestor is
                                            GO:0048856 (anatomical structure development) or GO:0003008 (system process).
The parsers are the authors' code (anatomy_audit.py, parse_obo / anc_fn) verbatim. No dates, no licence-restricted content.
"""
import json, sys
from collections import defaultdict
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
MECH_ROOTS = {"MONDO:7770011", "MONDO:0019052", "MONDO:0003778"}


def parse_obo(path, prefix):
    out = {}; cur = None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line == "[Term]": cur = {"is_a": [], "xref": [], "ns": None}; continue
        if line.startswith("[") and line.endswith("]"): cur = None; continue
        if cur is None or ": " not in line: continue
        k, v = line.split(": ", 1)
        if k == "id": cur["id"] = v; out[v] = cur
        elif k == "name": cur["name"] = v
        elif k == "namespace": cur["ns"] = v
        elif k == "is_a": cur["is_a"].append(v.split(" ")[0])
        elif k == "xref": cur["xref"].append(v.split(" ")[0])
        elif k == "relationship" and v.startswith("part_of ") and prefix == "GO:": cur["is_a"].append(v.split(" ")[1])
        elif k == "is_obsolete" and v == "true": cur["obs"] = True
    return {k: v for k, v in out.items() if k.startswith(prefix) and not v.get("obs")}


def anc_fn(tt):
    cache = {}
    def f(m):
        if m in cache: return cache[m]
        o = set(); st = list(tt.get(m, {}).get("is_a", []))
        while st:
            p = st.pop()
            if p in o or p not in tt: continue
            o.add(p); st.extend(tt[p]["is_a"])
        cache[m] = o; return o
    return f


def main(MONDO_OBO, GO_OBO):
    ids = json.load(open(REPO / "data/zenodo/model_outputs/cosine_matrix/ids_13021.json"))
    M = pd.read_csv(REPO / "data/raw/mondo/mondo_body_system_13021.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert list(M.disease_id) == ids
    terms = parse_obo(MONDO_OBO, "MONDO:"); anc = anc_fn(terms)
    dset = {d: (anc(m) | {m}) for d, m in zip(M.disease_id, M.mondo_id) if m}
    cnt = defaultdict(int)
    for v in dset.values():
        for t in v: cnt[t] += 1
    neo = [t for t, v in terms.items() if v.get("name") == "cancer or benign tumor"]; assert len(neo) == 1
    keep = {t for t, n in cnt.items() if 10 <= n <= 2000} | {neo[0]}
    rows = [(d, t) for d, v in dset.items() for t in sorted(v & keep)]
    out = REPO / "data/raw/mondo"
    pd.DataFrame(rows, columns=["disease_id", "term_id"]).to_csv(out / "mondo_terms_13021.tsv.gz", sep="\t", index=False)
    pd.DataFrame([(t, terms[t]["name"], bool(t in MECH_ROOTS or (anc(t) & MECH_ROOTS))) for t in sorted(keep)],
                 columns=["term_id", "name", "under_mech_root"]).to_csv(out / "mondo_terms.tsv", sep="\t", index=False)
    print("mondo", len(rows), len(keep))

    GO = parse_obo(GO_OBO, "GO:"); goanc = anc_fn(GO)
    g2 = pd.read_parquet(REPO / "data/raw/specific_channels/gene2go.parquet")
    used = set(g2[(g2.namespace == "biological_process") & g2.go_id.isin(GO)].go_id)
    allt = set(used)
    for t in used: allt |= goanc(t)
    ORG = {"GO:0048856", "GO:0003008"}
    go_out = REPO / "data/raw/go"; go_out.mkdir(exist_ok=True)
    pd.DataFrame([(t, GO[t]["name"], GO[t]["ns"], "|".join(sorted(goanc(t))), bool((goanc(t) | {t}) & ORG)) for t in sorted(allt)],
                 columns=["term_id", "name", "namespace", "ancestors", "organ_level"]).to_csv(go_out / "go_bp_terms.tsv.gz", sep="\t", index=False)
    print("go", len(allt))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

