# GO release-matched sensitivity (post hoc)

Small precomputed tables that cannot be regenerated inside this repository (upstream model runs or benchmark construction). Each line gives the content and the scripts that read the file.

- `grm_annotation_counts.tsv`: gene-GO pairs, annotated genes and GO terms per arm and slice (current vs slice-dated GO annotation files). Read by: si_go_release_table.py.
- `grm_common_node_census.tsv`: diseases, events and drugs on the common-node set per decision year. Read by: si_go_release_table.py.
- `grm_nmo_S2008_summary.tsv`: eculizumab PNH -> NMO matched percentile per arm (2008 slice, five seeds). Read by: si_go_release_table.py.
