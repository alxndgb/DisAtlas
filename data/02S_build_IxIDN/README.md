# IxIDN+ edges (stored reference)

- `ixidn_pos_edges.parquet.gz` — IxIDN+ disease pairs (43,112), the co-investigation graph before
  negative-set construction. Same pair set as `data/01_rebuilt/ixidn_pos_edges.parquet`, rebuilt by
  `code/01_build_ixidn/02_build_ixidn_pos.py`. Read by `code/04_train_model/01_strip_benchmark.py`
  and `04_evaluate_ensemble.py`.
