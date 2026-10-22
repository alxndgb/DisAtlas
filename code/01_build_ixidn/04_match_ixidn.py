"""Step 4 — IxIDN-Matched: cell-match judged positives to IxIDN- negatives.

Manuscript: Methods "IxIDN-Matched benchmark" — each kept positive (LLM judge verdict STRONG or
PLAUSIBLE) is matched to a negative in the same cell of trial degree (deciles of the lower and
higher endpoint degree), organ system, clinical area, name similarity, MONDO similarity and gene
coverage; 4,364 pairs (2,182 + 2,182).

Default: reads the stored judge verdicts (data/00_raw/llm_annotations/ixidn_judge_verdicts.parquet).
--rerun-llm: re-judges the candidates with code/02_llm_annotations/call_llm.py judge using YOUR
own API key (only this switch needs a key; LLM output is non-deterministic, so the benchmark can
differ from the published one). Not exercised by the repository tests.

Inputs: frozen matching features/candidates, disease_split.tsv.gz, disease_causal_genes.parquet,
diseases.parquet; data/01_rebuilt/ixidn_neg_edges.parquet.
Output: data/01_rebuilt/ixidn_matched.tsv — pair_uid, a, b, disease_a_name, disease_b_name, label,
pool, split_a, split_b, cellA, deg_a, deg_b, dsum, gmin, gmax, lgp, name_tfidf_cos, mondo_simgic,
genecov, shared_mondo_organsys, shared_ordo_area, verdict, n_drugs, max_phase, min_drug_deg.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ixidn_lib as L  # noqa: E402

KEYS = ["dlo10", "dhi10", "c_org", "c_ordo", "c_name4", "c_sg4", "genecov"]
KEEP_VERDICTS = ("STRONG", "PLAUSIBLE")


def cell_match(sub, keys, rng):
    k = np.zeros(len(sub), bool)
    y = sub.label.values
    cell = sub[keys].astype(str).agg("|".join, axis=1).values
    for cv in np.unique(cell):
        ip = np.where((cell == cv) & (y == 1))[0]
        inn = np.where((cell == cv) & (y == 0))[0]
        m = min(len(ip), len(inn))
        if m:
            k[rng.choice(ip, m, replace=False)] = True
            k[rng.choice(inn, m, replace=False)] = True
    return k


def rerun_verdicts(cand):
    """Untested here: re-judge with the user's own key via code/02_llm_annotations/call_llm.py."""
    tmp = L.out_path("ixidn_judge_rerun_input.csv")
    out = L.out_path("ixidn_judge_rerun.jsonl")
    cand.to_csv(tmp, index=False)
    subprocess.run([sys.executable, os.path.join(L.REPO_ROOT, "code", "02_llm_annotations", "call_llm.py"),
                    "judge", tmp, out], check=True)
    v = pd.read_json(out, lines=True)
    return v.rename(columns={"a": "disease_a", "b": "disease_b"})[["disease_a", "disease_b", "verdict"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rerun-llm", action="store_true",
                    help="re-run the judge prompts with your own API key instead of the stored verdicts")
    args = ap.parse_args()

    F = pd.read_parquet(L.MATCH_FEATURES)
    cg = pd.read_parquet(L.CAUSAL_GENES)
    tier = pd.to_numeric(cg.tier.astype(str).str.extract(r"(\d)")[0], errors="coerce")
    ng = cg[tier <= 2].groupby("disease_id").gene_symbol.nunique()
    na, nb = F.a.map(ng).fillna(0).values, F.b.map(ng).fillna(0).values
    F["gmin"], F["gmax"] = np.minimum(na, nb), np.maximum(na, nb)
    F["lgp"] = np.log1p(na) + np.log1p(nb)
    F["dsum"] = F.deg_a + F.deg_b
    F["k"] = [frozenset(x) for x in zip(F.a, F.b)]

    cand = pd.read_parquet(L.MATCH_CANDIDATES)
    V = rerun_verdicts(cand) if args.rerun_llm else pd.read_parquet(L.JUDGE_VERDICTS)
    vk = dict(zip(map(frozenset, zip(V.disease_a, V.disease_b)), V.verdict))
    cand["verdict"] = [vk.get(frozenset(x)) for x in zip(cand.a, cand.b)]
    assert cand.verdict.notna().all(), "unjudged candidates"
    keep = cand[cand.verdict.isin(KEEP_VERDICTS)]
    keepk = set(map(frozenset, zip(keep.a, keep.b)))

    sub = F[((F.label == 1) & F.k.isin(keepk)) | (F.label == 0)].reset_index(drop=True)
    M = sub[cell_match(sub, KEYS, np.random.default_rng(0))].copy()
    npos, nneg = int((M.label == 1).sum()), int((M.label == 0).sum())
    assert npos == nneg == len(keep), (npos, nneg, len(keep))

    neg = L.read_out("ixidn_neg_edges.parquet")
    okneg = set(map(frozenset, zip(neg.disease_1, neg.disease_2)))
    assert all(x in okneg for x in M.k[M.label == 0]), "a matched negative is not an IxIDN- pair"
    assert not M.duplicated(["a", "b"]).any()

    sp = pd.read_csv(L.DISEASE_SPLIT, sep="\t").set_index("disease_id").split
    M["split_a"], M["split_b"] = M.a.map(sp), M.b.map(sp)
    names = L.load_diseases().set_index("disease_id").label
    M["disease_a_name"], M["disease_b_name"] = M.a.map(names), M.b.map(names)
    assert M.disease_a_name.notna().all() and M.disease_b_name.notna().all()
    ck = cand.assign(k=[frozenset(x) for x in zip(cand.a, cand.b)]).set_index("k")
    for c in ("verdict", "n_drugs", "max_phase", "min_drug_deg"):
        M[c] = M.k.map(ck[c])
    cols = ["a", "b", "disease_a_name", "disease_b_name", "label", "pool", "split_a", "split_b",
            "cellA", "deg_a", "deg_b", "dsum", "gmin", "gmax", "lgp", "name_tfidf_cos",
            "mondo_simgic", "genecov", "shared_mondo_organsys", "shared_ordo_area",
            "verdict", "n_drugs", "max_phase", "min_drug_deg"]
    Mo = M[cols].sort_values(["cellA", "label", "a", "b"]).reset_index(drop=True)
    Mo.insert(0, "pair_uid", ["E%05d" % i for i in range(1, len(Mo) + 1)])
    name = "ixidn_matched_rerun_llm.tsv" if args.rerun_llm else "ixidn_matched.tsv"
    Mo.to_csv(L.out_path(name), sep="\t", index=False)
    print(f"IxIDN-Matched: {npos} positives + {nneg} negatives = {len(Mo)} pairs, "
          f"{pd.concat([Mo.a, Mo.b]).nunique()} diseases; pools {Mo.pool.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
