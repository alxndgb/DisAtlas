# si_cluster_themes

Per-cluster annotation of the K = 55 clustering of the 13,021 diseases, used by `tables/SI/table_cluster_themes.py`
(Supplementary Table SI_cluster_themes and the theme counts in the SI cluster-validation section and main text).
Ontology terms are public (MONDO, Gene Ontology); no licensed content.

- `bodysystem_by_cluster_long.parquet`: share and fold enrichment of each MONDO body-system class among the annotated members of each cluster (hypergeometric test, Benjamini–Hochberg `fdr`).
- `mechanism_terms_by_cluster.parquet`: enrichment of each MONDO mechanism class (inborn errors of metabolism, lysosomal storage, inborn errors of immunity and similar) per cluster.
- `go_bp_top5_by_cluster.parquet`: five highest-ranked Gene Ontology biological-process terms per cluster, with `organ_level` flagging organ- or system-level processes.
- `mondo_terms_top5_by_cluster.parquet`: five most enriched MONDO terms per cluster (table labels).
- `cluster_annotation_k55.parquet`: cluster size, direct neoplasm share and within-cluster pair enrichment over random pairs (phenotype, gene and GO biological-process overlap) used for the Mann–Whitney comparison of anatomy and process clusters.
- `theme_reference.tsv`: the theme printed for each cluster in the Supplementary Table; the generator asserts that its rule reproduces every row.
