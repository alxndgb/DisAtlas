# Input evidence sources (counts only)

Small precomputed table that cannot be regenerated inside this repository: it is counted on the integrated tier-1
disease-gene and disease-symptom association tables of the network build, which contain licensed and internal rows and
are not shipped. No association rows are included, only counts.

- `input_source_counts.json`: per-source association, disease, gene/symptom and single-source counts, the combined
  tier-1 totals after and before the per-disease cap, the MONDO-propagation layer size, and the gene/symptom overlap of
  the 13,080 candidate diseases. Read by: tables/SI/table_input_sources.py, tests/test_sibench_parity.py.
