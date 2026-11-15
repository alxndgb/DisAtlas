"""Maintainer-only: stage the frozen IxIDN-Matched matching table into data/raw/benchmarks/ixidn_matched_qc.
The stored QC verdicts themselves are staged by stage_benchmarks.py (ixidn_matched/llm_curation_verdicts.tsv).
Source locations are passed as KEY=PATH and are never tracked.

usage: python scripts/maintainer/stage_matched_qc.py <data dir> universe=<matching universe parquet>

Rows: every veto-passing IxIDN- pair (y=0) and the 4,302 QC candidate positives (y=1), in the original row order
(the seeded cell match depends on it). cell = the seven matching keys (degree deciles of the lower- and
higher-degree disease, shared organ system, shared Orphanet area, name-similarity and semantic-similarity
bins, gene-coverage flag) joined by '|'. pair_id = the blinded id of the pair in ixidn_matched_pairs.tsv.
"""
import os
import sys

import pandas as pd

data = sys.argv[1]
src = dict(a.split("=", 1) for a in sys.argv[2:])
B = os.path.join(data, "raw", "benchmarks")
out = os.path.join(B, "ixidn_matched_qc")
os.makedirs(out, exist_ok=True)

bm = pd.read_csv(os.path.join(B, "ixidn_matched", "ixidn_matched_pairs.tsv"), sep="\t")
ids = {frozenset(k): p for k, p in zip(zip(bm.disease_a, bm.disease_b), bm.pair_id)}
ori = {frozenset(k): k for k in zip(bm.disease_a, bm.disease_b)}
V = pd.read_csv(os.path.join(B, "ixidn_matched", "llm_curation_verdicts.tsv"), sep="\t")
cand = set(map(frozenset, zip(V.disease_a, V.disease_b)))

F = pd.read_parquet(src["universe"])
K = ["dlo10", "dhi10", "c_org", "c_ordo", "c_name4", "c_sg4", "genecov"]
assert (F[K].astype(str).agg("|".join, axis=1) == F.cellA).all()
F = F[(F.label.values == 0) | [frozenset(k) in cand for k in zip(F.a, F.b)]]
k = [ori.get(frozenset(x), x) for x in zip(F.a, F.b)]  # benchmark pairs keep the benchmark orientation
U = pd.DataFrame({"disease_a": [x[0] for x in k], "disease_b": [x[1] for x in k],
                  "y": F.label.astype("int8").values, "cell": F.cellA.values})
U["pair_id"] = [ids.get(frozenset(x)) for x in k]
assert int(U.y.sum()) == len(V) and U.pair_id.notna().sum() == len(bm)
U.to_parquet(os.path.join(out, "matching_pairs.parquet"), index=False, compression="zstd")
print(len(U), int(U.y.sum()))
