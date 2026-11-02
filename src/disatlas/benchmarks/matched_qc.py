"""COLD (opt-in): rerun the IxIDN-Matched language-model QC of the 4,302 candidate positive pairs.

WARM runs never call this module: they read the stored verdicts in
data/raw/benchmarks/ixidn_matched/llm_curation_verdicts.tsv. This module rebuilds the QC inputs described in
the Methods (both disease names, the shared trial drugs and their highest trial phase, a four-level scale:
strong, plausible, weak, suspect) and asks a model of your choice for one verdict per pair, in batches of 20.
Pairs without a shared trial drug in the IxIDN graph are assigned NOT_IN_GRAPH without a call.

    ANTHROPIC_API_KEY=... python -m disatlas.benchmarks.matched_qc --out build/qc --model MODEL_ID
    python -m disatlas.benchmarks.matched_qc --out build/qc --dry-run      # write the prompts only

Inputs (shipped): the candidate pairs (pairs of llm_curation_verdicts.tsv), the shared drugs per pair
(data/zenodo/figures_main/ixidn_graph_edges.parquet), disease names (--names TSV id<TAB>name, default
data/raw/ontology/mondo.obo). Drug names and highest trial phase come from the licensed trial source and are
not redistributed; pass them with --drugs TSV (drug_id<TAB>name) and --phases TSV (disease_a, disease_b,
max_phase) if you hold them, otherwise drugs are given by their public identifiers and phase is omitted.
The prompt below follows the Methods description; the authors' original prompt text was not stored, so a rerun
reproduces the stored verdicts approximately, not exactly. Build the benchmark from the rerun with
python -m disatlas.benchmarks.ixidn_matched --verdicts build/qc/llm_curation_verdicts.tsv --out ...
"""
import argparse
import json
import os
import re
import time
import urllib.request

import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
VERDICTS = os.path.join(ROOT, "data", "raw", "benchmarks", "ixidn_matched", "llm_curation_verdicts.tsv")
EDGES = os.path.join(ROOT, "data", "zenodo", "figures_main", "ixidn_graph_edges.parquet")
MONDO = os.path.join(ROOT, "data", "raw", "ontology", "mondo.obo")
SCALE = ("STRONG", "PLAUSIBLE", "WEAK", "SUSPECT")
BATCH = 20

SYSTEM = (
    "You review disease pairs from a clinical-trial co-investigation benchmark. Two diseases form a candidate "
    "positive pair because the same drugs were tested in trials for both. For each pair, judge whether the "
    "shared-drug evidence reflects a genuine therapeutic relationship between the two diseases.\n"
    "STRONG: the shared drugs act on a mechanism relevant to both diseases.\n"
    "PLAUSIBLE: a shared therapeutic rationale is likely but less specific.\n"
    "WEAK: the overlap is explained by broad-use drugs (supportive care, anti-infectives, general "
    "immunosuppression) or by co-occurring conditions rather than a shared disease mechanism.\n"
    "SUSPECT: the trial link is probably an artefact (mis-mapped disease, umbrella or basket trial, "
    "complication or adverse-event indication).\n"
    'Answer with a JSON array only: [{"i": <pair number>, "verdict": "STRONG|PLAUSIBLE|WEAK|SUSPECT"}, ...]')


def obo_names(path):
    names, cur = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("id: "):
                cur = line[4:].strip()
            elif line.startswith("name: ") and cur:
                names.setdefault(cur, line[6:].strip())
            elif line.startswith("[") and not line.startswith("[Term]"):
                cur = None
    return names


def table(path, k, v):
    if not path:
        return {}
    d = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(d[k], d[v])) if isinstance(k, str) else {frozenset(x): p for x, p in zip(zip(*[d[c] for c in k]), d[v])}


def evidence(names_path, drugs_path, phases_path):
    C = pd.read_csv(VERDICTS, sep="\t")[["disease_a", "disease_b"]]
    E = pd.read_parquet(EDGES)
    shared = {frozenset(k): json.loads(s) for k, s in zip(zip(E.disease_1, E.disease_2), E.shared_drugs)}
    names = table(names_path, "id", "name") if names_path else (obo_names(MONDO) if os.path.exists(MONDO) else {})
    if not names_path:  # names of the benchmark diseases also ship with the Figure 3 data
        D = pd.read_csv(os.path.join(ROOT, "data", "derived", "fig3_ixidn", "ixidn_matched_pairs.tsv"), sep="\t")
        for i, n in list(zip(D.a, D.disease_a_name)) + list(zip(D.b, D.disease_b_name)):
            names.setdefault(i, n)
    dn, ph = table(drugs_path, "drug_id", "name"), table(phases_path, ["disease_a", "disease_b"], "max_phase")
    C["name_a"] = C.disease_a.map(names).fillna(C.disease_a)
    C["name_b"] = C.disease_b.map(names).fillna(C.disease_b)
    C["drugs"] = [[dn.get(x, x) for x in shared.get(frozenset(k), [])] for k in zip(C.disease_a, C.disease_b)]
    C["max_phase"] = [ph.get(frozenset(k), "") for k in zip(C.disease_a, C.disease_b)]
    return C


def prompt(batch):
    lines = []
    for i, r in enumerate(batch.itertuples(), 1):
        p = f"; highest trial phase {r.max_phase}" if r.max_phase else ""
        lines.append(f"{i}. {r.name_a} | {r.name_b} | shared trial drugs: {', '.join(r.drugs)}{p}")
    return "\n".join(lines)


def call(model, text, key, retries=4):
    body = json.dumps({"model": model, "max_tokens": 2000, "system": SYSTEM,
                       "messages": [{"role": "user", "content": text}]}).encode()
    req = urllib.request.Request("https://api.anthropic.com/v1/messages", body, {
        "x-api-key": key, "anthropic-version": "2023-06-01", "content-type": "application/json"})
    for t in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return "".join(b.get("text", "") for b in json.load(r)["content"])
        except Exception:
            if t == retries - 1:
                raise
            time.sleep(10 * (t + 1))


def parse(text, n):
    out = {int(d["i"]): str(d["verdict"]).upper() for d in json.loads(re.search(r"\[.*\]", text, re.S).group(0))}
    return [out.get(i) if out.get(i) in SCALE else None for i in range(1, n + 1)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--model", help="model identifier (required unless --dry-run)")
    ap.add_argument("--names"), ap.add_argument("--drugs"), ap.add_argument("--phases")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    C = evidence(a.names, a.drugs, a.phases)
    C["verdict"] = None
    C.loc[C.drugs.str.len() == 0, "verdict"] = "NOT_IN_GRAPH"
    todo = C[C.verdict.isna()]
    batches = [todo.iloc[i:i + BATCH] for i in range(0, len(todo), BATCH)]
    with open(os.path.join(a.out, "qc_prompts.txt"), "w") as fh:
        fh.write(SYSTEM + "\n")
        for b in batches:
            fh.write("\n---\n" + prompt(b) + "\n")
    if a.dry_run:
        print(f"{len(todo)} pairs in {len(batches)} prompts; {int((C.verdict == 'NOT_IN_GRAPH').sum())} NOT_IN_GRAPH")
        return
    key = os.environ.get("ANTHROPIC_API_KEY")
    if not key or not a.model:
        raise SystemExit("--cold-qc needs ANTHROPIC_API_KEY in the environment and --model (or QC_MODEL)")
    for j, b in enumerate(batches):
        v = parse(call(a.model, prompt(b), key), len(b))
        miss = [k for k, x in enumerate(v) if x is None]
        if miss:  # one retry for unparsed pairs
            v2 = parse(call(a.model, prompt(b.iloc[miss]), key), len(miss))
            for k, x in zip(miss, v2):
                v[k] = x
        C.loc[b.index, "verdict"] = v
        print(f"batch {j + 1}/{len(batches)}", flush=True)
    V = C[["disease_a", "disease_b", "verdict"]].assign(rejudged=False, first_pass_verdict="",
                                                       keep=C.verdict.isin(("STRONG", "PLAUSIBLE")))
    V.to_csv(os.path.join(a.out, "llm_curation_verdicts.tsv"), sep="\t", index=False)
    print(V.verdict.value_counts(dropna=False).to_dict())


if __name__ == "__main__":
    main()
