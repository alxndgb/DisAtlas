# Figure 7 and SI event-time odds-ratio summaries

Small precomputed tables behind Figure 7, Supplementary Figures S18-S20 and Supplementary Tables S12-S17. They are staged from the event-time analysis by `scripts/maintainer/stage_fig7_or.py`; the conditional-logit fits were run on the licensed per-candidate trial rows with `scripts/maintainer/fig7_clogit_fit.R` and ship as fitted odds ratios only. Each line gives the content and the scripts that read the file.

- `or_clogit.tsv`: conditional-logit odds ratios per 1 SD of each score, by decision year and pooled, per model. Read by: fig7_prospective.py, si_event_time_or_tables.py.
- `nmo_case_clogit.json`, `nmo_case_scatter.tsv`, `nmo_candidate_sets.tsv`, `nmo_hops.json`: eculizumab case (decision year 2007): NMO rank, candidate pools and graph distance. Read by: fig7_prospective.py, si_event_time_or_tables.py (not the scatter).
- `nmo_ms_by_slice.tsv`, `nmo_ms_yearly_all_methods.tsv`, `fig7_nmo_ms_yearly_ens5_test.tsv`: NMO and multiple sclerosis ranks by evidence slice. Read by: si_event_time_or_tables.py (first), fig7_prospective.py (other two).
- `nmosd_timeline.tsv`, `ecu_nmo_trial_records.tsv`: sourced milestones and trial records of the eculizumab case. Read by: fig7_prospective.py, si_event_time_or_tables.py (trial records).
- `abs_percentiles.tsv`, `auc_per_drug.tsv`, `auc_check.tsv`, `auc_check_yearly_box.tsv`, `auc_deltas.tsv`, `auc_dropped.tsv`: median percentiles and per-drug probabilities (AUROC) by method and decision year. Read by: fig7_prospective.py, si_event_time_auc_methods.py, si_event_time_or_tables.py.
- `e2_per_drug.tsv`, `e2_per_drug_from2000.tsv`, `weighted_yearly.tsv`: top-10% share per drug and decision year, unweighted and popularity-weighted. Read by: fig7_prospective.py (`--si hits10`, `--si weighted`), si_event_time_or_tables.py.
- `early_years_metrics.tsv`, `small_sets_by_year.tsv`: early decision years and small comparison sets. Read by: si_event_time_or_tables.py.
- `merge_and_matching_census.tsv`: event merging and candidate-matching census. Read by: fig7_prospective.py.
- `decision_event_time.json`: the pre-specified decision rule and its outcome. Read by: si_event_time_or_tables.py.
