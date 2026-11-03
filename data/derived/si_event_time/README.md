# Event-time (time-sliced) evaluation summaries

Small precomputed tables that cannot be regenerated inside this repository (upstream model runs or benchmark construction). Each line gives the content and the scripts that read the file.

- `coverage.tsv`: coverage of new trials by slice, pool and comparison set. Read by: si_cluster_rarecommon_numbers.py, si_event_time_numbers.py, si_mechanism_coverage_slices.py, si_rare_common_neighbour_coverage.py.
- `deltas.tsv`: DisAtlas minus reference-score differences in median percentile (`dE1`, percentile points) and top-10% share (`dE2`, share of drugs), with drug-cluster bootstrap intervals. Read by: si_event_time_numbers.py, si_event_time_tables.py.
- `merge_and_matching_census.tsv`: event merging and candidate-matching census. Read by: si_event_time_numbers.py.
- `point_estimates_all.tsv`: point estimates, all-disease pool. Read by: si_event_time_numbers.py.
- `point_estimates_test.tsv`: point estimates, test-disease pool. Read by: si_event_time_numbers.py.
- `posthoc_covid_2020_counts.tsv`: post hoc 2020 COVID-19 event counts. Read by: si_event_time_numbers.py.
- `posthoc_covid_2020_dE1.tsv`: post hoc change in median percentile (`dE1`) excluding COVID-19 events. Read by: si_event_time_numbers.py.
- `table_pooled_top10_precision.tsv`: pooled top-10% precision table. Read by: si_event_time_numbers.py.
