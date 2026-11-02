# Disease-cluster validation (K = 55 on 13,021 diseases)

Small precomputed tables that cannot be regenerated inside this repository (upstream model runs or benchmark construction). Each line gives the content and the scripts that read the file.

- `SI_controls_table.tsv`: positive and negative cluster control pairs with best ranks and same-cluster flags. Read by: si_cluster_rarecommon_numbers.py.
- `dmd_replacement_candidates.tsv`: candidate replacement diseases for the Duchenne/DM1 control pair, from the cluster-controls analysis. Read by: si_cluster_rarecommon_numbers.py.
- `k_selection_boot.tsv.gz`: K sweep: silhouette per K over 13,021-disease subsamples (bootstrap of the clustering), used for K selection. Read by: si_cluster_rarecommon_numbers.py, si_cluster_validation.py.
- `raw_vs_learned_silhouette_k55.json`: silhouette at K = 55 in cosine and learned-distance space. Read by: si_cluster_rarecommon_numbers.py.
- `sifu_cluster_13021_values.json`: internal and external cluster-validation statistics on the 13,021-disease clustering. Read by: si_cluster_map_vs_model.py, si_cluster_rarecommon_numbers.py, si_cluster_validation.py.
- `within_between_K55_13021.tsv.gz`: within- versus between-cluster pair scores at K = 55. Read by: si_cluster_validation.py.
