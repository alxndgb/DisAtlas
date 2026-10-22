# frozen

Frozen benchmark inputs, shipped as used.

- `ixidn_eval_pair_universe_frozen.tsv.gz` — IxIDN evaluation pair universe with matching factors, original row order; columns row, a, b, label, deg_a, deg_b, shared_mondo_organsys, shared_ordo_area, name_tfidf_cos, mondo_simgic, genecov; 224,469 rows.
- `disease_split.tsv.gz` — tune/test assignment of the 985 benchmark diseases; columns disease_id, split, split_node_disjoint.
- `ixidn_matched_universe_features.parquet` — IxIDN-Matched candidate universe with matching cells, original row order; 182,790 rows, 28 columns.
- `ixidn_matched_positive_candidates.parquet` — IxIDN-Matched positive candidates; columns pair_id, a, b, n_drugs, max_phase, min_drug_deg; 4,302 rows.
- `disease_causal_genes.parquet` — disease causal-gene sets for the shortcut audit; columns disease_id, gene_symbol, tier; 16,063 rows.
- `disease_coinvestigated_drugs.parquet` — drugs co-investigated per IxIDN disease, own trials or direct MONDO parent/child; columns disease_id, drug_id, via (own / hier); 285,773 rows.
- `pinned_train_strip_pairs.parquet` — IxIDN pairs removed from the training graph; columns a, b, class (pos / neg); 456,968 rows.
- `ixidn_eval_pairs.parquet` — IxIDN evaluation pairs after the co-investigation veto; columns a, b, class (pos / neg); 378,669 rows.
