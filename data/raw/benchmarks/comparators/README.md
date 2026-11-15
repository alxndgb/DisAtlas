# data/raw/benchmarks/comparators

Inputs of the IxIDN-Matched comparator benchmark (Figure 6, Table 1). Read by `src/disatlas/rebuild/fig6_inputs.py` and
`tests/test_comparator_scores.py`; staged by `scripts/maintainer/stage_benchmarks.py`.

| file | content |
|---|---|
| `ixidn_matched_blind_input.tsv` | the 4,364 IxIDN-Matched pairs as given to every comparator: blind_id, MONDO ids `a`, `b`, disease names; no labels |
| `ixidn_matched_comparator_scores.tsv` | `pair_id` (= blind_id) and one column per comparator run: the `DDA` score of that run on each pair. Column-to-original-file map and NA counts: `data/README.md` ("Column provenance of comparators/ixidn_matched_comparator_scores.tsv"); per-column checksums in `tests/test_comparator_scores.py` |

## TxGNN: 213 NA pairs, kept as unmappable

The TxGNN column (`TxGNN published ckpt|DistMult`, published checkpoint, not retrained) is NA for 213 pairs:

- 187 pairs: at least one of the two diseases has no TxGNN disease node (`mapped` false in the original `txgnn_pair_scores.tsv.gz`).
- 26 pairs, "collapsed": both diseases map to the same TxGNN disease node (`idx_a` = `idx_b`). The original file carries a DistMult
  value for them, but it is a self-score of one node (embedding cosine 1), not a comparison of two diseases, so they are set to NA.

All 213 are ranked last in every analysis (one value below the lowest scored pair; `data/derived/fig5_clustering/fig5c_txgnn_pair_scores.tsv`
stores them as the sentinel -9.929e-10). The 4,151 remaining pairs are the "can map" subset quoted in the Results.
