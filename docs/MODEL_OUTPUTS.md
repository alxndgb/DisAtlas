# data/zenodo/model_outputs

Frozen outputs of the models behind the manuscript. Figures, tables and numbers are computed from these files, so they
reproduce exactly without retraining. Every file is a byte copy, or a row/column subset, of the authors' run outputs.
`MANIFEST.tsv` lists the source role of each file and `SHA256SUMS` lists its checksum.

| path | content |
|---|---|
| `embeddings/gene_anchored/embedding_s<seed>.npz` | final (gene-anchored) disease embedding, seeds 4707-4711; `ids` (13,021 diseases), `X` float32 (unit rows, dim 32) |
| `cosine_matrix/cosine_ens5_13021.npy` | 13,021 x 13,021 float32 mean cosine over the 5 seeds, rows/columns in `ids_13021.json` order (Figure 5 and SI clustering) |
| `scores/ixidn_matched_cosine_gene_anchored.parquet` | IxIDN-Matched pairs (4,364): ensemble and per-seed cosine (`c`, `c_s<seed>`), node-disjoint placement (`nd_c`) |
| `scores/node_disjoint_placement_cosine_gene_anchored.parquet` | node-disjoint placement cosines, fold x seed |
| `scores/ixidn_full_pair_scores_gene_anchored.parquet` | full-pool IxIDN pair scores |
| `clusters/disease_clusters_k55.parquet` | per-disease K = 55 cluster label and Isomap coordinates |
| `clusters/cluster_summary_k55.parquet` | per-cluster theme, organ-system, mechanism and GO enrichment tests |
| `pfull_reference/pfull_ixidn_cosine_5seed.parquet` | the P_full training-graph retrain: IxIDN-Matched pairs (`blind_id`, `a`, `b`, `label`, `pool`) and per-seed cosines (`Pfull_cos_s<seed>`), replayed by `scripts/pfull_parity.py` |

Licensed identifiers are blinded (`blind_id`); no licensed name, drug or rationale string is shipped.

No file carries a decoder score: DisAtlas pair scores are 5-seed mean embedding cosines.
