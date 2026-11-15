"""Maintainer-only: stage the per-event percentile inputs of the event-time analysis (Figure 7, Figures S19-S22,
Tables S13-S16) from the pre-registered scoring run into data/derived/si_event_time/.
    python scripts/maintainer/stage_event_time.py RANK_DIR OUT_DIR
RANK_DIR holds the et02 outputs etrank_y<year>_<method>_<pool>_sh<shard>.tsv.gz and etprec_*.tsv.gz (pool all / test).
Written files: event_percentiles.parquet (one row per pool x decision year x drug (ChEMBL id) x disease (MONDO id) x method:
score, six percentile analyses, coverage flag, merge/matching counters) and event_top10_precision.parquet (per drug x pool x
method x population stratum top-10% precision). No trial ids and no trial dates are carried; the year is the decision
(first-trial) year of the cleared triplet."""
import sys, glob
import pandas as pd
RD, OUT = sys.argv[1], sys.argv[2]
AN = ["pct_p1", "pct_s1", "pct_s2", "pct_s3", "pct_s4", "pct_s5"]
COLS = ["y", "compound_id", "event_v3", "method", "score", "covered", "n_events_drug_year", "p1_merge_level", "p1_small", "s2_widen", "s2_short", "sd_s2", "n_cand"] + AN
parts, prec = [], []
for pool in ("all", "test"):
    for f in sorted(glob.glob(f"{RD}/etrank_y*_{pool}_sh*.tsv.gz")):
        d = pd.read_csv(f, sep="\t", usecols=COLS); d.insert(0, "pool", pool); parts.append(d)
    for f in sorted(glob.glob(f"{RD}/etprec_y*_{pool}_sh*.tsv.gz")):
        d = pd.read_csv(f, sep="\t", dtype={"pop_stratum": str}, usecols=["compound_id", "method", "pop_stratum", "prec10", "y"]); d.insert(0, "pool", pool); prec.append(d)
X = pd.concat(parts, ignore_index=True); PR = pd.concat(prec, ignore_index=True)
for c in ("covered", "p1_small"): X[c] = X[c].astype(bool)
for c in ("p1_merge_level", "s2_widen", "s2_short", "n_events_drug_year", "n_cand", "y"): X[c] = X[c].astype("int32")
X = X.sort_values(["pool", "y", "compound_id", "event_v3", "method"], kind="stable").reset_index(drop=True)
PR = PR.sort_values(["pool", "y", "compound_id", "method", "pop_stratum"], kind="stable").reset_index(drop=True)
print(X.shape, PR.shape, round(X.memory_usage().sum() / 1e6), "MB")
X.to_parquet(f"{OUT}/event_percentiles.parquet", compression="zstd", index=False)
PR.to_parquet(f"{OUT}/event_top10_precision.parquet", compression="zstd", index=False)
