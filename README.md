# DisAtlas v3.2

Code and data to reproduce the figures, tables and numbers of the revised DisAtlas manuscript and its
Supplementary Information. Data layout, column dictionaries and licensing: `data/README.md`.

Setup: Python 3.11 with `pip install -r requirements.txt`; download the Zenodo bundle into `data/zenodo`
(see `data/README.md`). All commands run from the repository root.

## Warm run (default)

```
./run_pipeline.sh                    # stages: verify graph figures si tests notebooks
./run_pipeline.sh --from figures     # resume at a stage;  --only STAGE runs one stage
```

Reads `data/raw`, `data/derived` and the processed Zenodo data in `data/zenodo` (training graph, frozen
5-seed model, benchmark and all-pairs scores), after checking it against `data/ZENODO_SHA256SUMS.txt`.
The graph stage rebuilds the training graph from `data/raw/channels` and checks it edge for edge against
the Zenodo graph. Every number is **exact**: the tests compare each manuscript and SI number with the
value recomputed from these files.

Measured on one cluster node (16 cores, 160 GB): 14 min wall, 8.7 GB peak memory (tests 1 min); up to 46 min
on a node shared with other jobs.

## Cold run

```
./run_pipeline.sh --cold             # stages: verify channels graph retrain score figures si tests notebooks
./run_pipeline.sh --cold --from graph
```

Builds the processed data from `data/raw` into `build/` (gitignored): channels (`build/channels`), training
graph (`build/graph`), the 5-seed ensemble retrained with the published recipe (`build/model`), and the
benchmark and all-pairs scores (`build/scores`), then compares `build/` with the frozen record in
`data/zenodo`. The channel stage needs `disatlas.channels.broad/specific`; without them, start at
`--from graph`, which uses the shipped `data/raw/channels`. Figures, tables, notebooks and the manuscript
numbers are always produced from the frozen record (they need `data/zenodo`), because a retrained ensemble
reproduces the published scores closely but not bit for bit.

- **Exact** (`tests/test_cold_exact.py`): channel files equal the shipped ones; broad 13,021 diseases /
  332,082 edges, specific 422,804 links / 191,019 edges, union 498,155 and training 488,103 edges, and the
  training edge set equals the Zenodo graph (SHA-256 and edge for edge).
- **Close** (`tests/test_cold_tolerance.py`): the retrain repeats the frozen recipe with the same random
  streams (walks, decoder pairs, validation split and initial weights are identical); only floating-point
  rounding differs (torch thread count, BLAS build), so embeddings agree to about 2e-7 but not bit for bit.
  Benchmark AUROCs of the retrained ensemble must lie within 0.005 of the frozen ensemble, and ensemble pair
  scores must be rank-correlated with the frozen ones (Spearman >= 0.99 on the 4,364 IxIDN-Matched pairs and
  on 10^6 random pairs). Observed on the cluster (seeds 4707-4711, retrains with 32, 8 and 3 torch threads,
  the last from this pipeline on 16 cores): |AUROC delta| <= 0.0004 on IxIDN-Matched, its test, rare-rare and
  rare-common splits and IxIDN+/-, <= 0.0008 on the FDA dev controls; Spearman >= 0.99996 on IxIDN-Matched
  and >= 0.9997 on random pairs (two frozen seeds agree with each other at 0.964-0.967 on IxIDN-Matched).
  A run is deterministic for a fixed thread count: two 32-thread runs of seed 4707 gave bit-identical
  embeddings.
- The warm test suite then runs as in a warm run, so the manuscript numbers stay exact.

Measured on one cluster node (16 cores, 160 GB), `--cold --from graph`: 37 min wall, 10.4 GB peak memory
(graph 1 min; retrain 16 min with the five seeds in parallel, about 2 GB each; score 7 min, mostly the
all-pairs matrix; figures, tables, tests and notebooks 13 min); up to 2 h on a node shared with other jobs.
