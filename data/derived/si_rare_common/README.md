# Rare-common stratified analyses

Small precomputed tables that cannot be regenerated inside this repository (upstream model runs or benchmark construction). Each line gives the content and the scripts that read the file.

- `ccplus_strata.tsv`: IxIDN-Matched AUROC by rare/common stratum. Read by: si_cluster_rarecommon_numbers.py, si_rare_common_transfer.py.
- `common_anchor_rank_common_partners.tsv`: documented-partner ranks by anchor and partner class. Read by: si_rare_common_transfer.py.
- `common_anchor_rank_rare_partners.tsv`: documented-partner ranks by anchor and partner class. Read by: si_rare_common_transfer.py.
- `coverage_by_k_summary.tsv`: share of diseases with an other-class disease among k nearest, observed and three nulls (1,000 permutations). Read by: si_cluster_rarecommon_numbers.py, si_rare_common_neighbour_coverage.py.
- `enrichment_other_class_frac10.tsv`: mean other-class fraction among 10 nearest neighbours, observed vs nulls. Read by: si_cluster_rarecommon_numbers.py.
- `four.tsv`: AUROC by stratum and method. Read by: figstyle_si.py, si_cluster_rarecommon_numbers.py, si_fda_prioritisation_text.py, si_rare_common_transfer.py.
- `headline_test_auroc.tsv`: main-text headline test-split AUROC and disease-cluster bootstrap interval, as quoted in the SI. Read by: si_cluster_rarecommon_numbers.py.
- `ixidn_negative_pool_funnel.tsv`: IxIDN benchmark construction: within-stratum pairs, post-veto negative pool and selected negatives per stratum (model independent). Read by: si_cluster_rarecommon_numbers.py.
- `partner_percentile_per_pair.tsv`: per-pair partner percentiles for the 68 test rare-common pairs (input to the two-way disease bootstrap). Read by: si_cluster_rarecommon_numbers.py.
- `partner_percentile_per_pair_txgnn.tsv`: TxGNN partner percentile for the same 68 pairs and two directions (published checkpoint, scored as in the IxIDN-Matched comparison; pairs it cannot map ranked last; TxGNN_pair_mapped flags the 65 mapped pairs). Written by scripts/maintainer/txgnn_partner_percentile.py. Read by: si_cluster_rarecommon_numbers.py.
- `partner_universe.tsv`: candidate universe of the partner percentiles (disease_id in score-matrix order, disease_class). Written by scripts/maintainer/stage_partner_universe.py. Read by: scripts/maintainer/txgnn_partner_percentile.py.
- `per_disease_first_other_class.tsv.gz`: rank of the first other-class neighbour per disease. Read by: si_rare_common_neighbour_coverage.py.
- `r4mc3.tsv`: stratified method comparison table. Read by: si_rare_common_transfer.py.
- `rare_anchor_rank_common_partners.tsv`: documented-partner ranks by anchor and partner class. Read by: si_rare_common_transfer.py.
- `rare_anchor_rank_rare_partners.tsv`: documented-partner ranks by anchor and partner class. Read by: si_rare_common_transfer.py.
- `rc_coverage_meta.json`: permutation count and median first-common rank. Read by: si_cluster_rarecommon_numbers.py.
