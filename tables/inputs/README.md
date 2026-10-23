# Table inputs

Processed outputs used only by the table scripts (one row per pair or disease; identifiers as in
`data/00_raw/`).

- `comparator_run_scores.tsv` — per-run scores of DisAtlas (per seed and ensemble), the LLM comparators (stored outputs), TxGNN (DistMult) and Phylo (seeds 42-46) on the 4,364 IxIDN-Matched pairs; blank = unscored/unmapped.
- `sensitivity_post_test_scores.tsv` — per-seed probabilities on the 1,094 test-split pairs for the DisAtlas reference and the four sensitivity arms.
- `cluster_labels_by_K.tsv` — seed-0 k-medoids cluster labels for every network disease, K = 10..80.
- `cluster_control_pairs.csv` — the 17 cluster-validation control pairs, their ensemble score, independence flag and percentile among IxIDN-Matched pairs.
- `ixidn_matched_ccplus_scores.tsv` — DisAtlas ensemble scores on the supplementary common-common benchmark (IxIDN-Matched plus 309 added matched negative common-common pairs).
- `rare_common_stratified_bootstrap.tsv`, `rare_common_ccplus_bootstrap.tsv` — stored disease-cluster bootstrap summaries (AUROC, 95% CI, DisAtlas minus matched-shortcut-stack delta).
- `hpo_baselines.csv` — stored HPO-similarity baseline AUROCs on the HPO-scorable pairs with DisAtlas deltas and bootstrap CIs.
- `source_display_labels.csv` — public display label for each opaque source code in `data/00_raw/associations/`.
- `si_fig3a_trials.tsv` — the 90 trials behind the ten disease-drug lines of Figure 3a; columns trial_key (as in `data/00_raw/trials/trial_triples.parquet`), disease_id, drug_id, public_id (public registry identifier, blank when none), phase (as printed), ctgov_phase (phase on the ClinicalTrials.gov record, blank when not looked up), title.
- `umbrella_placement.tsv` — the 444 umbrella disease terms (at least five lexical subtypes): number of subtypes, whether the subtype-average position differs beyond projection noise and lies in another cluster (93), and the released alternative-position flag for the 69 such diseases without own symptom data (subtype-average K = 50 cluster, map shift, projection noise).
