# Tables

Scripts that regenerate every table of the manuscript body and its Appendix from repository data.
Each script computes every cell, writes `<table>.csv` and `<table>.md` beside itself, and compares
every cell with the value printed in the manuscript (a mismatch stops the run).

    python tables/build_all.py          # from the repository root

| Manuscript table | Section | Script | Inputs |
|---|---|---|---|
| tbl:protocols | Methods, Comparator systems | `main/table_protocols.py` | `inputs/comparator_run_scores.tsv` |
| tbl:SI_input_dis_gene_sources | Appendix, input data sources | `SI/table_SI_input_dis_gene_sources.py` | `data/00_raw/associations/`, `data/00_raw/diseases/`, `inputs/source_display_labels.csv` |
| tbl:SI_input_dis_symptom_sources | Appendix, input data sources | `SI/table_SI_input_dis_symptom_sources.py` | as above |
| tbl:SI_ontology_baseline (HPO rows) | Appendix, Ontology-only classifier and concordance checks | `SI/table_SI_ontology_baseline.py` | `data/08S_eval/ixidn4364_scores.tsv`, `data/00_raw/associations/disease_symptom.parquet`, `inputs/hpo_baselines.csv` |
| tbl:SI_cluster_controls (+ co-clustering ratios in the text) | Appendix, The full control list | `SI/table_SI_cluster_controls.py` | `inputs/cluster_labels_by_K.tsv`, `inputs/cluster_control_pairs.csv`, `data/08S_eval/ixidn4364_scores.tsv` |
| tbl:SI_rare_common_strata | Appendix, rare/common Summary table | `SI/table_SI_rare_common_strata.py` | `data/08S_eval/ixidn4364_scores.tsv`, `inputs/ixidn_matched_ccplus_scores.tsv`, `inputs/rare_common_*_bootstrap.tsv` |
| tbl:SI_sensitivity_post_test | Appendix, Sensitivity to settings fixed after test-split results | `SI/table_SI_sensitivity_post_test.py` | `inputs/sensitivity_post_test_scores.tsv`, `data/08S_eval/ixidn4364_scores.tsv` |
| tbl:SI_fig3a_trials | Appendix, Trials behind the Figure 3a example | `SI/table_SI_fig3a_trials.py` | `inputs/si_fig3a_trials.tsv`, `data/00_raw/trials/trial_triples.parquet`, `data/00_raw/diseases/`, `data/00_raw/drugs/` |

No script calls a language model or downloads anything; LLM comparator scores are the stored
outputs of the runs reported in the paper.

## Notes on specific cells

- **Input-source tables.** The column printed as "Diseases exclusively from this source" counts
  tier-1 associations supported by that source alone. Percentages are diseases / 13,080, where
  13,080 is the candidate universe: MechDist-universe diseases that share at least one gene or HPO
  term (any evidence tier) with another disease. The MONDO-propagation row is evidence tier 9 of
  `disease_symptom.parquet`.
- **Ontology baseline.** Only the two HPO rows are rebuilt. The manuscript's DisAtlas cell for the
  HPO-scorable pairs repeats the baseline range; the rebuilt value is 0.749 (DisAtlas higher by
  0.03-0.06, not significant), as stated in the Appendix prose. HPO similarity scores are read from
  the stored baseline output. The ORDO-dependent rows (least-concordant subset, outside-Orphanet
  control) and the Spearman caption sentence are not rebuilt.
- **Bootstrap intervals** (rare/common CIs, shortcut deltas, HPO deltas) are read from stored
  bootstrap outputs; point estimates are recomputed from scores.
- **Cluster controls.** K = 50 labels are the shipped seed-0 clustering. NF1-mutant/NF2 is separated at
  13 of the 14 other K values (the manuscript says "13 of the 15 other K values"; the grid has 15 K values, so 14 others).
