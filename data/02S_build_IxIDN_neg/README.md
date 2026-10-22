# IxIDN- edges (stored reference)

- `ixidn_neg_veto_edges.parquet.gz` — IxIDN- disease pairs after the co-investigation veto (335,557).
  Same pair set as `data/01_rebuilt/ixidn_neg_edges.parquet`, rebuilt by
  `code/01_build_ixidn/03_build_ixidn_neg.py`.
- `ixidn_neg_pinned_edges.parquet.gz` — the 413,856 negative candidates (shortest-path distance >= 3)
  before the veto. Read by `code/04_train_model/01_strip_benchmark.py` and `04_evaluate_ensemble.py`.
