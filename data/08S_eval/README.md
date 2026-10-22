# DisAtlas scores on IxIDN-Matched

- `ixidn4364_scores.tsv` — frozen DisAtlas ensemble score `p` (and per-seed `p_s*`), label `y` and
  `pool` for the 4,364 IxIDN-Matched pairs, as scored by the frozen model in `data/zenodo/models/`.
  AUROC is 0.765 on all pairs and 0.782 on the held-out test split; `tests/test_ixidn_parity.py`
  and `benchmarks/verify_benchmarks.py` check these values. It is the shared score input of the
  benchmark, ontology-baseline and shortcut-baseline scripts.
- `ixidn4364_ci.tsv` — 95% disease-cluster bootstrap intervals of the ensemble AUROC on all pairs
  (0.729-0.799), the held-out test split (0.721-0.841) and the rare-rare, rare-common and
  common-common strata, as reported in the Results and Supplementary Information. Read by
  `notebooks/04_headline_evaluation.ipynb`.
