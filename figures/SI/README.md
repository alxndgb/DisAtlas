# Appendix (SI) figures

Each script rebuilds one appendix figure from repo-relative inputs, asserts the values the appendix text reports, and writes `<figure>.png` (300 dpi) and `<figure>.pdf` next to itself. Shared styling is in `_style.py`; figure-level source tables are in `source_data/`.

Run from the repo root, e.g. `python figures/SI/si_uncertainty.py`. Scripts marked `--from-model` additionally recompute their distances from the frozen model package in `data/zenodo/models/` (optional; needs torch and safetensors).

| Figure | Script | Inputs | Appendix section |
|---|---|---|---|
| ablation_aggregation_rule.png | si_ablation_aggregation.py | source_data/ablation_retrained_arms_auroc.tsv, aggregation_edges_kept.json | Aggregation-rule ablation (#sec:SI_ablations_aggregation, #fig:SI_ablation_aggregation) |
| ablation_evidence_sources.png | si_ablation_evidence.py | source_data/ablation_retrained_arms_auroc.tsv, ablation_direct_distance_baselines.csv, ablation_decoder_vs_cosine.csv, ablation_arm_graph_coverage.csv, ablation_coverage.tsv | Component ablations / evidence-source ablation (#sec:SI_ablations_evidence, #fig:SI_ablation_evidence) |
| ablation_go_measure.png | si_ablation_go_measure.py | source_data/go_measure_*.tsv/json | GO gene-pair distance rule (#sec:SI_ablations_gomeasure, #fig:SI_ablation_gomeasure) |
| extended_cluster_validation.png | si_cluster_validation.py | source_data/cluster_k_selection_bootstrap.tsv.gz, cluster_within_between_k50.tsv.gz, cluster_mondo_concordance.json | Cluster validation (#fig:SI_cluster_validation) |
| SI_cluster_map_vs_model_distance.png | si_layout_fidelity.py | figures/main/data/fig3_display_map.tsv.gz, source_data/layout_fidelity_random_pairs.tsv.gz (optional --from-model: data/zenodo/models) | Display layout fidelity (#fig:SI_cluster_validation_layout) |
| figS4_lca_similarity_diameter_revised.png | si_lca_formulas.py | none (formula sheet) | Threshold and distance-metric sensitivity (#fig:SI_MechDist) |
| ixidn_matched_funnel.png | si_ixidn_matched_funnel.py | data/00_raw/frozen/{pinned_train_strip_pairs.parquet, ixidn_eval_pairs.parquet, ixidn_eval_pair_universe_frozen.tsv.gz, ixidn_matched_universe_features.parquet, ixidn_matched_positive_candidates.parquet}, data/00_raw/llm_annotations/ixidn_judge_verdicts.parquet, data/08S_eval/ixidn4364_scores.tsv; writes source_data/ixidn_matched_funnel_counts.tsv, ixidn_matched_hop_distribution.tsv | IxIDN-Matched benchmark construction (#fig:SI_ixidn_matched_funnel) |
| mechanism_inventory_sensitivity.png | si_mechanism_inventory.py | source_data/inventory_*.tsv/json | Mechanism inventory sensitivity (#fig:SI_mechanism_sensitivity) |
| mechdist_degree_distribution.png | si_mechdist_degree_distribution.py | data/03_network/mechdist_network_{edges,nodes}.parquet.gzip | IxIDN-Matched benchmark construction (#fig:SI_mechdist_degree_dist) |
| mechdist_schematic.png | si_mechdist_schematic.py | data/03_network/mechdist_network_edges.parquet.gzip, source_data/mechdist_training_strip_edges.tsv.gz, mechdist_schematic_pair_distance.tsv, mechdist_schematic_example.json (optional --from-model: data/zenodo/models) | From the MechDist network to DisAtlas (#fig:SI_mechdist_embedding_schematic) |
| negative_set_realism.png | si_negative_set_realism.py | source_data/negative_realism_ordo_distance.tsv.gz, negative_realism_hard_negative_auroc.tsv | Negative-set realism (#fig:SI_negrealism) |
| node_disjoint_robustness.png | si_node_disjoint.py | source_data/node_disjoint_auroc.tsv | Node-disjoint evaluation (#fig:SI_node_disjoint) |
| SI_placement_confidence.png | si_placement_confidence.py | source_data/placement_exclusion_sensitivity.tsv | Per-disease placement confidence (#fig:SI_placement_confidence) |
| SI_rare_common_transfer.png | si_rare_common_transfer.py | source_data/rare_common_strata_auroc.tsv, rare_common_anchored_groups.tsv | Rare-common stratified evaluation (#fig:SI_rare_common_stratified) |
| threshold_sensitivity.png | si_threshold_sensitivity.py | source_data/threshold_sweep.tsv, threshold_reference_disatlas_ensemble.tsv | Robustness to the cutoff value (#fig:SI_threshold_sensitivity) |
| uncertainty_quantification.png | si_uncertainty.py | source_data/uncertainty_seed_vs_sample.tsv, uncertainty_disatlas_ci.tsv | Uncertainty quantification (#fig:SI_uncertainty) |

Not script-generated: `LCA.png` (#fig:SI_lca) is a hand-drawn illustration of ontology distance; `rare_disease_stats.png` (#fig:SI_rare_disease_stats) is reproduced from an external source credited in its caption. Both ship as static images with the manuscript.
