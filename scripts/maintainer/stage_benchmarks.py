"""Maintainer-only: stage the manuscript benchmarks (minimal columns) into data/raw/benchmarks, comparator
and DisAtlas benchmark scores into data/zenodo/scores, and precomputed event-time result tables into
data/derived/event_time. Source locations are passed on the command line as KEY=PATH (never tracked).

usage: python scripts/maintainer/stage_benchmarks.py <data dir> scores4364=... per_run=... nodedisjoint=...
       judge=... labels=... event_time=... ordo_pairs=... ordo_summary=... fda=...
"""
import json
import os
import shutil
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
data = sys.argv[1]
src = dict(a.split("=", 1) for a in sys.argv[2:])
B = os.path.join(data, "raw", "benchmarks")
Z = os.path.join(data, "zenodo", "scores")
for d in (B + "/ixidn_matched", B + "/comparators", Z, os.path.join(data, "derived", "event_time")):
    os.makedirs(d, exist_ok=True)
SEEDS = [4707, 4708, 4709, 4710, 4711]
log = {}

# IxIDN-Matched (4,364 pairs)
S = pd.read_csv(src["scores4364"], sep="\t")
assert len(S) == 4364 and int(S.y.sum()) == 2182
assert (S.split_a.eq("test") & S.split_b.eq("test")).eq(S.pool.eq("test")).all()
bm = pd.DataFrame({"pair_id": S.blind_id, "disease_a": S.a, "disease_b": S.b, "y": S.y.astype(int),
                   "split": S.pool, "pair_class": S.rare_common_pair})
bm.to_csv(B + "/ixidn_matched/ixidn_matched_pairs.tsv", sep="\t", index=False)
log["ixidn_matched"] = dict(pairs=len(bm), pos=int(bm.y.sum()), test=int((bm.split == "test").sum()),
                            test_pos=int(((bm.split == "test") & (bm.y == 1)).sum()),
                            diseases=int(pd.concat([bm.disease_a, bm.disease_b]).nunique()))

# DisAtlas (full) scores on IxIDN-Matched, reproduced by the shipped ensemble
from disatlas.model import load_ensemble  # noqa: E402

m = load_ensemble(os.path.join(data, "zenodo", "model", "disatlas_v3.2_ensemble"))
P = m.score_pairs(bm.disease_a, bm.disease_b, per_seed=True)
diff = float(np.nanmax(np.abs(P - S[[f"p_s{s}" for s in SEEDS]].values.T)))
log["ixidn_matched_full_max_abs_diff_vs_stored"] = diff
assert diff < 1e-5, diff
full = bm[["pair_id"]].copy()
for i, s in enumerate(SEEDS):
    full[f"p_s{s}"] = S[f"p_s{s}"].values
full["p_ensemble"] = S.p.values
full.to_csv(Z + "/ixidn_matched_disatlas_full_scores.tsv", sep="\t", index=False)

# node-disjoint arm
ND = pd.read_csv(src["nodedisjoint"], sep="\t")
assert (ND.blind_id.values == S.blind_id.values).all()
nd = bm[["pair_id"]].copy()
for s in SEEDS:
    nd[f"p_s{s}"] = ND[f"DisAtlas node-disjoint|s{s}"].values
nd["p_ensemble"] = nd[[f"p_s{s}" for s in SEEDS]].mean(axis=1)
nd.to_csv(Z + "/ixidn_matched_disatlas_nodedisjoint_scores.tsv", sep="\t", index=False)

# stored comparator outputs (language-model, agentic and knowledge-graph systems), per run
R = pd.read_csv(src["per_run"], sep="\t")
assert (R.blind_id.values == S.blind_id.values).all()
ren = {}
for c in R.columns:
    if "|" not in c or c.startswith("DisAtlas"):
        continue
    sysname, run = c.split("|", 1)
    run = run.replace(".tsv", "").replace("ixidn_expanded_BLIND_with_DDA_", "")
    ren[c] = f"{sysname}|{run}"
cmp_ = pd.concat([bm[["pair_id"]], R[list(ren)].rename(columns=ren)], axis=1)
cmp_.to_csv(B + "/comparators/ixidn_matched_comparator_scores.tsv", sep="\t", index=False)
log["comparator_columns"] = list(ren.values())

# language-model curation of benchmark positives (verdicts as used to build IxIDN-Matched)
J = pd.read_parquet(src["judge"])
J = J.rename(columns={"prior_verdict": "first_pass_verdict"})
J["rejudged"] = J.verdict_source.eq("expanded_pass_rescored")
J[["disease_a", "disease_b", "verdict", "rejudged", "first_pass_verdict", "keep"]].to_csv(
    B + "/ixidn_matched/llm_curation_verdicts.tsv", sep="\t", index=False)
log["llm_curation"] = dict(candidates=len(J), kept=int(J.keep.sum()), rejudged=int(J.rejudged.sum()),
                           verdicts=J.verdict.value_counts().to_dict())

# disease labels used for rare/common and cancer stratification
L = pd.read_parquet(src["labels"])
L[["disease_id", "is_rare", "is_cancer", "is_malignant", "disease_class"]].to_csv(
    B + "/disease_labels.tsv", sep="\t", index=False)

# ORDO / MONDO ontology distance for each IxIDN-Matched pair (-1 = not defined), and the summary of
# ORDO distance for random disease pairs used as the reference distribution
if "ordo_pairs" in src:
    O = pd.read_csv(src["ordo_pairs"], sep="\t")
    assert (O.blind_id.values == S.blind_id.values).all()
    os.makedirs(B + "/ordo", exist_ok=True)
    O.rename(columns={"blind_id": "pair_id"}).to_csv(B + "/ordo/ixidn_matched_ontology_distance.tsv",
                                                     sep="\t", index=False)
    Q = pd.read_csv(src["ordo_summary"], sep="\t").drop(columns=["source"])
    Q.to_csv(B + "/ordo/ordo_distance_summary.tsv", sep="\t", index=False)
    log["ordo_defined_pairs"] = int((O.ordo_dist >= 0).sum())

# FDA-approved repurposing positive controls (holdout split excluded): stratum, mechanism class, split and FDA label verification
# (the training-graph strip uses only the pair list in fda_repurposing_pairs.tsv)
if "fda" in src:
    F = pd.read_csv(src["fda"], sep="\t").rename(columns={"mondo_a": "disease_a", "mondo_b": "disease_b"})
    F = F[["disease_a", "disease_b", "stratum", "pair_class", "split", "fda_verified_pair", "in_model"]]
    F = F[F.split.ne("holdout")]  # the holdout is sealed: no holdout rows or labels are released
    F.to_csv(B + "/fda_repurposing/fda_repurposing_controls.tsv", sep="\t", index=False)
    log["fda"] = dict(pairs=len(F), split=F.split.value_counts().to_dict(),
                      dev_verified_in_model=int((F.split.eq("dev") & F.fda_verified_pair & F.in_model).sum()))

# event-time evaluation: precomputed result tables (inputs are proprietary trial evidence)
ET = src["event_time"]
for f in ("coverage.tsv", "coverage_by_year.tsv", "deltas.tsv",
          "merge_and_matching_census.tsv", "per_year_vs_wAA.tsv", "point_estimates_all.tsv",
          "point_estimates_test.tsv", "table_pooled_E1_E2.tsv", "table_pooled_top10_precision.tsv"):
    shutil.copyfile(os.path.join(ET, f), os.path.join(data, "derived", "event_time", f))
print(json.dumps(log, indent=1, default=str))
