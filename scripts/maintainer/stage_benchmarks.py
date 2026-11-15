"""Maintainer-only: stage the manuscript benchmarks (minimal columns) and the stored comparator scores into
data/raw/benchmarks. DisAtlas benchmark scores are the embedding cosines in data/zenodo/model_outputs/scores. Source locations are passed on the command line as KEY=PATH (never tracked).

usage: python scripts/maintainer/stage_benchmarks.py <data dir> scores4364=... per_run=...
       judge=... labels=... ordo_pairs=... ordo_summary=... fda=...
"""
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
data = sys.argv[1]
src = dict(a.split("=", 1) for a in sys.argv[2:])
B = os.path.join(data, "raw", "benchmarks")
for d in (B + "/ixidn_matched", B + "/comparators"):
    os.makedirs(d, exist_ok=True)
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

# ORDO / MONDO ontology distances of the IxIDN-Matched pairs are no longer staged: no released code read them (release/MINIMISATION.tsv).
# The ontology baseline (Table S4/S5) reads data/derived/ixidn_robustness/ixidn4364_ontology_distance.tsv.

# FDA-approved repurposing positive controls (holdout split excluded): stratum, mechanism class, split and FDA label verification
# (the training-graph strip uses only the pair list in fda_repurposing_pairs.tsv)
if "fda" in src:
    F = pd.read_csv(src["fda"], sep="\t").rename(columns={"mondo_a": "disease_a", "mondo_b": "disease_b"})
    F = F[["disease_a", "disease_b", "stratum", "pair_class", "split", "fda_verified_pair", "in_model"]]
    F = F[F.split.ne("holdout")]  # the holdout is sealed: no holdout rows or labels are released
    F.to_csv(B + "/fda_repurposing/fda_repurposing_controls.tsv", sep="\t", index=False)
    log["fda"] = dict(pairs=len(F), split=F.split.value_counts().to_dict(),
                      dev_verified_in_model=int((F.split.eq("dev") & F.fda_verified_pair & F.in_model).sum()))

print(json.dumps(log, indent=1, default=str))
