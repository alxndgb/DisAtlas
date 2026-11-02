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

Measured on one cluster node (16 cores): 46 min wall, 8.6 GB peak memory (tests 4 min).

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
reproduces the published scores only statistically.

- **Exact** (`tests/test_cold_exact.py`): channel files equal the shipped ones; broad 13,021 diseases /
  332,082 edges, specific 422,804 links / 191,019 edges, union 498,155 and training 488,103 edges, and the
  training edge set equals the Zenodo graph (SHA-256 and edge for edge).
- **Statistical** (`tests/test_cold_tolerance.py`): retraining uses different random streams from the
  frozen run, so scores are not bit-identical. Benchmark AUROCs of the retrained ensemble must lie within
  0.03 of the frozen ensemble (0.04 for the IxIDN-Matched test split), and ensemble pair scores must be
  rank-correlated with the frozen ones (Spearman >= 0.90 on the 4,364 IxIDN-Matched pairs, >= 0.65 on
  10^6 random pairs). Observed on the cluster (seeds 4707-4711): AUROC deltas -0.015 (IxIDN-Matched),
  -0.022 (test), -0.003 (rare-rare), -0.014 (rare-common), +0.014 (IxIDN+/-), -0.007 (FDA dev controls);
  Spearman 0.946 on IxIDN-Matched and 0.726 on random pairs (frozen seeds agree with each other at
  0.75-0.77 on random pairs, retrained seeds at 0.62-0.64). Two independent retrains gave the same deltas
  within 0.001.
- The warm test suite then runs as in a warm run, so the manuscript numbers stay exact.

Measured on one cluster node (16 cores), `--cold --from graph`: 1 h 59 min wall, 10.4 GB peak memory
(graph 3 min; retrain 63 min with the five seeds in parallel, about 1.6 GB each; score 16 min, mostly the
all-pairs matrix; figures, tables, tests and notebooks 38 min). One seed alone takes about 15 min on 8 cores.
