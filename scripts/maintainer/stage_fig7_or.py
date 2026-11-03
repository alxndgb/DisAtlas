"""Maintainer-only: stage the stored inputs of Figure 7 (conditional-logit design), the SI event-time figures and the
SI event-time tables into data/derived/fig7_prospective.  Nothing is recomputed: every file is a column-trimmed copy of a
stored result.  Source locations are passed as KEY=PATH and are never tracked.

usage: python scripts/maintainer/stage_fig7_or.py <data dir> KEY=PATH ...
keys:  clogit_yearly clogit_pooled clogit_pm0      conditional-logit fits (scripts/cold/fig7_clogit, opt-in cold path)
       abs auc_per_drug auc_deltas auc_check auc_box auc_dropped e2 e2_2000 weighted early small_sets
       nmo_case nmo_scatter nmo_sets nmo_hops nmo_ms_slice ecu_trials decision
Licensed trial-registry identifiers and record dates are dropped; the drug of the illustrative case (eculizumab) is named in
the manuscript, so its case files carry no compound identifier.
"""
import json
import os
import sys

import pandas as pd

data = sys.argv[1]
src = dict(a.split("=", 1) for a in sys.argv[2:])
out = os.path.join(data, "derived", "fig7_prospective")
os.makedirs(out, exist_ok=True)
rd = lambda k, **kw: pd.read_csv(src[k], sep="\t", **kw)


def wr(df, name):
    df.to_csv(os.path.join(out, name), sep="\t", index=False)
    print(name, df.shape)


# conditional-logit odds ratios: efron ties, full sample; one row per (set, model, term)
OR_COLS = ["set", "model", "term", "term_type", "OR", "lo95", "hi95", "coef", "se_robust", "p_robust",
           "n_cases", "n_drugs", "n_strata", "n_rows", "converged"]
A = pd.concat([rd(k, dtype={"set": str}) for k in ("clogit_yearly", "clogit_pooled", "clogit_pm0")])
A = A[(A["sample"].astype(str) == "full") & (A["method"] == "efron")].copy()
A["set"] = A["set"].replace({"pooled_2000_2021": "pooled"})
A = A[A.term.astype(str).str.startswith("z_")]
assert not A.duplicated(["set", "model", "term"]).any(), "duplicate odds-ratio rows"
wr(A[OR_COLS], "or_clogit.tsv")

for k, name in [("abs", "abs_percentiles.tsv"), ("auc_per_drug", "auc_per_drug.tsv"), ("auc_deltas", "auc_deltas.tsv"),
                ("auc_check", "auc_check.tsv"), ("auc_box", "auc_check_yearly_box.tsv"), ("auc_dropped", "auc_dropped.tsv"),
                ("e2", "e2_per_drug.tsv"), ("e2_2000", "e2_per_drug_from2000.tsv"), ("small_sets", "small_sets_by_year.tsv"),
                ("nmo_scatter", "nmo_case_scatter.tsv"), ("nmo_sets", "nmo_candidate_sets.tsv")]:
    wr(rd(k, dtype=str, keep_default_na=False), name)
W = rd("weighted")
wr(W[(W.weights.isin(["wex", "w51"])) & (W.trimmed.astype(str) == "False")], "weighted_yearly.tsv")
E = rd("early", dtype=str, keep_default_na=False)
wr(E[["pool", "set", "y", "slice", "n_drugs", "n_events", "events_per_drug", "N_median", "share_N_lt10", "E2", "E2_lo",
      "E2_hi", "chance", "auc", "auc_lo", "auc_hi"]], "early_years_metrics.tsv")
M = rd("nmo_ms_slice", dtype=str, keep_default_na=False)
wr(M, "nmo_ms_by_slice.tsv")
# eculizumab -> NMO spectrum disorder in the trial table: public registry records only (source, NCT id, phase, first year)
T = rd("ecu_trials", dtype=str, keep_default_na=False)
T = T[(T.kind == "nct_present") | (T.source == "orphanet")]
T = T.assign(source=T.source.where(T.kind == "record", "ixidn_trial_table"))[["kind", "nct_id", "source", "phase", "first_year_rule"]]
wr(T, "ecu_nmo_trial_records.tsv")

CJ = json.load(open(src["nmo_case"]))
keep = ["year", "t", "nmo_id", "nmo_in_test_pool", "N", "disatlas_rank", "disatlas_ties", "disatlas_pct_rank",
        "popularity_rank", "popularity_ties", "popularity_pct_rank", "yearly_table_pct_s4", "yearly_table_n_cand"]
json.dump({k: CJ[k] for k in keep}, open(os.path.join(out, "nmo_case_clogit.json"), "w"), indent=1)
H = json.load(open(src["nmo_hops"]))["stripped_training_graph"]
json.dump({"stripped_training_graph": {k: H[k] for k in ("n_nodes", "n_edges", "direct_edge", "hops", "pnh_degree", "nmo_degree",
                                                         "shared_neighbours")}},
          open(os.path.join(out, "nmo_hops.json"), "w"), indent=1)
D = json.load(open(src["decision"]))
json.dump({k: D[k] for k in ("rule", "ci", "passes", "n_pass", "verdict")}, open(os.path.join(out, "decision_event_time.json"), "w"), indent=1)
