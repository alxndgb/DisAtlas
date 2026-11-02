# DisAtlas

Code and data to reproduce the results of the manuscript

> **A Mechanistic Disease Representation Framework for Therapeutic Knowledge Transfer Across Rare and Common Conditions**

and its Supplementary Information, including the trial co-investigation benchmarks **IxIDN+**, **IxIDN-** and
**IxIDN-Matched**. The scope is the revised manuscript and its Supplementary Information: every main and
Supplementary figure and table that is computed from data, and every number in the main text, is regenerated
by one command and checked by the test suite.

DisAtlas places 13,021 diseases (10,489 rare, 2,532 common) in one embedding space learned from a
disease-disease network (MechDist) with two evidence channels: a **broad** channel (shared gene function, phenotype,
cell-type expression and protein domains) and a **specific** channel (shared mechanism modules). A five-seed
node2vec ensemble with a Hadamard MLP decoder scores any disease pair.

## Contents

1. [Manuscript to repository map](#manuscript-to-repository-map)
2. [Repository structure](#repository-structure)
3. [Quick start](#quick-start)
4. [Run modes](#run-modes)
5. [Data and licensing](#data-and-licensing)
6. [Tests and parity](#tests-and-parity)
7. [Troubleshooting](#troubleshooting)
8. [How to cite](#how-to-cite)
9. [Licence](#licence)

## Manuscript to repository map

Every script runs from the repository root without arguments (`./run_pipeline.sh` calls them in order) and writes
its outputs next to itself. Most scripts also write a `*_values.tsv` / `*_values.json` / `*_numbers.json` file with
every plotted or quoted number; the tests read those files.

### Main text

| Item | Content | Script | Output |
|---|---|---|---|
| Figure 1 | Knowledge transfer between rare and common diseases (illustration) | none (static image) | `figures/main/fig1_knowledge_transfer.png` |
| Figure 2 | DisAtlas pipeline and its two evaluations | `figures/main/fig2_disatlas_pipeline.py` | `figures/main/disatlas_pipeline.{png,pdf}` |
| Figure 3 | The IxIDN-Matched benchmark | `figures/main/fig3_ixidn_benchmark.py` | `figures/main/ixidn_benchmark.{png,pdf}` |
| Figure 4 | Embedding and decoder architecture | `figures/main/fig4_upna_ddal_architecture.py` | `figures/main/UPNA_DDAL_architecture.{png,pdf,svg}` |
| Figure 5 | Disease map, clusters and rarity-stratified performance | `figures/main/fig5_disease_clustering_map.py` | `figures/main/disease_clustering_map.{png,pdf}` |
| Figure 6 | Cross-system comparison on IxIDN-Matched | `figures/main/fig6_method_comparison.py` | `figures/main/fig5_method_comparison.{png,pdf}` (file name kept from the submission) |
| Figure 7 | Simulated prospective (event-time) test | `figures/main/fig7_prospective.py` | `figures/main/fig7_prospective.{png,pdf}` |
| Table 1 | Protocols of the cross-system comparison | `tables/main/table1_protocols.py` | `tables/main/table1_protocols.{md,tsv}` |
| Main-text numbers | All 536 number tokens of the main text, each with its kind and check | `tests/main_quantities.py`, `tests/test_main_parity.py` | registry `tests/main_text_numbers.tsv` |

### Supplementary figures

| Item | Content | Script (`figures/SI/`) | Output (`figures/SI/`) |
|---|---|---|---|
| Figure S1 | Prevalence of rare diseases (reproduced from the literature) | none (third-party image, not redistributed) | not shipped |
| Figure S2 | Detailed DisAtlas pipeline | `si_pipeline_detailed.py` | `si_pipeline_detailed.{png,pdf}` |
| Figure S3 | Distance to lowest common ancestor (schematic) | none | not shipped |
| Figure S4 | From the MechDist network to embeddings (schematic) | none | not shipped |
| Figure S5 | Construction of IxIDN-Matched and its shortcut audit | `si_ixidn_matched_funnel.py` | `si_ixidn_matched_funnel.{png,pdf}` |
| Methods, Figure 3b, Figure S5 | IxIDN-Matched quality control: language-model verdicts on the 4,302 candidate positives (1,285 strong + 897 plausible kept; repeat check 84.7%) and the seeded cell match | `src/disatlas/benchmarks/ixidn_matched.py` (stored verdicts `data/raw/benchmarks/ixidn_matched/llm_curation_verdicts.tsv`), rerun `src/disatlas/benchmarks/matched_qc.py`, `tests/test_ixidn_matched_qc.py` | `data/derived/fig3_ixidn/qc_verdict_summary.tsv` (kept/removed by reason) |
| Figure S6 | Broad and specific evidence channels | `si_mechanism_sensitivity.py` | `si_mechanism_sensitivity.{png,pdf}` |
| Figure S7 | Degree in each graph of the MechDist network | `si_mechdist_degree_by_channel.py` | `si_mechdist_degree_by_channel.{png,pdf}` |
| Figure S8 | Threshold sensitivity | `si_threshold_sensitivity.py` | `si_threshold_sensitivity.{png,pdf}` |
| Figure S9 | Component ablations (31 arms) | `si_ablation_31arm.py` | `si_ablation_31arm.{png,pdf}` |
| Figure S10 | Seed and sampling uncertainty | `si_uncertainty.py` | `si_uncertainty.{png,pdf}` |
| Figure S11 | Per-disease placement confidence | `si_placement_confidence.py` | `si_placement_confidence.{png,pdf}` |
| Figure S12 | Node-disjoint evaluation | `si_node_disjoint_negrealism.py` | `si_node_disjoint.{png,pdf}` |
| Figure S13 | Negative-set realism and hard negatives | `si_node_disjoint_negrealism.py` | `si_negrealism.{png,pdf}` |
| Figure S14 | Extended cluster validation | `si_cluster_validation.py` | `extended_cluster_validation.{png,pdf}` |
| Figure S15 | Display layout fidelity | `si_cluster_map_vs_model.py` | `SI_cluster_map_vs_model_distance.{png,pdf}` |
| Figure S16 | Performance by disease rarity and anchor disease | `si_rare_common_transfer.py` | `SI_rare_common_transfer.{png,pdf}` |
| Figure S17 | Rare-common neighbour coverage | `si_rare_common_neighbour_coverage.py` | `SI_rare_common_neighbour_coverage.{png,pdf}` |
| Figure S18 | Later-trialled diseases within popularity strata | `si_popularity_strata.py` | `si_popularity_strata.{png,pdf}` |
| Figure S19 | Approved-repurposing benchmark (schematic) | `si_fda_schematic.py` | `SI_approved_repurposing_schematic.{png,pdf}` |
| Figure S20 | Evidence channels on FDA-approved repurposing pairs | `si_fda_channels.py` | `SI_fda_channels.{png,pdf}` |
| Figure S21 | Mechanism coverage of approved expansions in dated graph slices | `si_mechanism_coverage_slices.py` | `SI_mechanism_coverage_slices.{png,pdf}` |
| Figure S22 | Clinical milestones of complement inhibitors | `si_complement_milestones.py` | `SI_complement_milestones.{png,pdf}` |
| Figure S23 | Evidence behind three rare-common case pairs | `si_bio_gap_numbers.py` (case scores) | `SI_explainability.png` (static), `SI_explainability_case_scores.tsv` |
| Figure S24 | How the training graph is built and scored | `si_graph_channels_schematic.py` | `SI_graph_channels_schematic.{png,pdf}` |
| Figure S25 | Test of compositionality | `si_compositionality_recompute.py` (tables), `si_compositionality.py` (figure) | `SI_compositionality.{png,pdf}` |

### Supplementary tables

| Item | Content | Script | Output |
|---|---|---|---|
| Table S1 | Disease-gene association sources | none (descriptive) | in the SI text |
| Table S2 | Disease-symptom association sources | none (descriptive) | in the SI text |
| Table S3 | Trials behind Figure 3a | `tables/SI/table_fig3a_trials.py` | `tables/SI/table_fig3a_trials.{md,tsv}` |
| Table S4 | Ontology-only baselines | `tables/SI/table_ontology_baseline.py` | `tables/SI/table_ontology_baseline.md` |
| Table S5 | Node-disjoint DisAtlas vs language-model systems by subgroup | `tables/SI/table_llm_subgroups.py` | `tables/SI/table_llm_subgroups.{md,tsv}` |
| Table S6 | Degree-matched and within-disease AUROC | `tables/SI/table_degree_controls.py` | `tables/SI/table_degree_controls.{md,tsv}` |
| Table S7 | Positive and negative control pairs | `figures/SI/si_cluster_rarecommon_numbers.py` | `figures/SI/SI_cluster_rarecommon_numbers.json` (key `controls`) |
| Table S8 | AUROC by disease-rarity pairing | `figures/SI/si_cluster_rarecommon_numbers.py` | `figures/SI/SI_cluster_rarecommon_numbers.json` |
| Tables S9, S10 | Event-time holdout: pooled advantages; subsets | `figures/SI/si_event_time_tables.py` | `figures/SI/SI_event_time_tables.md` |
| Table S11 | Advantage with current and slice-dated GO annotations | `figures/SI/si_go_release_table.py` | `figures/SI/SI_go_release_table.{md,tsv}` |
| Table S12 | FDA-approved pairs scored with dated models | `figures/SI/si_fda_prioritisation_text.py` | `figures/SI/SI_fda_prioritisation_text.md` |
| SI numbers | Numbers quoted in the SI text | `figures/SI/si_*_numbers.py` | `figures/SI/SI_*_numbers.json`; asserted by `tests/test_sibench_parity.py`, `tests/test_sibio_parity.py` (specs in `tests/parity/`) |

## Repository structure

```
run_pipeline.sh          one command for figures, tables, tests and notebooks (see Run modes)
requirements.txt         pinned Python packages
data/
  raw/channels/broad/    broad channel: disease-gene and disease-phenotype annotations and broad edges (post-threshold)
  raw/channels/specific/ specific channel: disease-module links, disease-target-mechanism evidence, mechanism classes
  raw/benchmarks/        IxIDN+/IxIDN- pairs, IxIDN-Matched pairs and curation verdicts, FDA controls, comparator outputs
  raw/ontology/          ontology (OBO) files read by the cold broad-channel rebuild
  derived/               small precomputed result tables read by the figure, table and test scripts
  zenodo/                (download, not in git) training graph, frozen 5-seed model, score matrices, large tables
  README.md              column dictionaries, licence and blinding
  ZENODO_SHA256SUMS.txt  checksums of every file in data/zenodo
src/disatlas/            library: graph build, model loading and scoring, retraining, blinding, paths, figure style
figures/main, figures/SI scripts for the main and Supplementary figures, and their outputs
tables/main, tables/SI   scripts for Table 1 and the Supplementary tables, and their outputs
notebooks/               four exploratory notebooks (data and channels, IxIDN benchmarks, SI robustness, SI biology)
tests/                   parity tests and the main-text number registry
scripts/                 notebook runner; scripts/maintainer/ shows how the shipped data were exported (not needed to reproduce)
```
`data/raw/ontology/hp.obo` is a field slice of the Human Phenotype Ontology read by the cold broad-channel rebuild:
20,413 term stanzas with only `id`, `alt_id`, `is_a`, `is_obsolete`, `replaced_by` and `consider` (no header,
release tag, dates, names or definitions). The Gene Ontology is not shipped: the GO-, HPA- and Pfam-based pair
distances ship precomputed in `data/zenodo/channels_raw/broad/pair_gene_distances.parquet`.

## Quick start

Requirements: Linux or macOS, Python 3.11, about 16 GB of RAM, and about 1 GB of disk for the Zenodo bundle (TBD once the channel inputs are added).

```bash
git clone <repository URL> DisAtlas && cd DisAtlas
python3.11 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Download the Zenodo record (DOI: `ZENODO_DOI_TBD`, <https://doi.org/ZENODO_DOI_TBD>) and unpack it so that
`data/zenodo/` contains `INDEX.tsv`, `graph/`, `model/`, `scores/`, `figures_main/`, `si_bio/`. Then check it:

```bash
cd data/zenodo && sha256sum -c ../ZENODO_SHA256SUMS.txt && cd ../..   # macOS: shasum -a 256 -c
./run_pipeline.sh
```

`run_pipeline.sh` uses `.venv/bin/python` (override with `PYTHON=/path/to/python`) and sets `PYTHONPATH=src`.
To run a single script by hand: `export PYTHONPATH=$PWD/src MPLBACKEND=Agg` and then, for example,
`.venv/bin/python figures/main/fig3_ixidn_benchmark.py`.

## Run modes

| Mode | Command | What it does | Agreement with the manuscript |
|---|---|---|---|
| WARM (default) | `./run_pipeline.sh` | Uses the processed data and frozen 5-seed model from `data/zenodo` (training graph, model, benchmark and all-pairs scores); rebuilds the training graph and IxIDN-Matched as checks | exact |
| COLD | `./run_pipeline.sh --cold` | Rebuilds the channels and the training graph from `data/raw` into `build/`, retrains the 5-seed ensemble and re-scores the benchmarks, compares `build/` with the frozen record, then makes all figures and tables and runs the tests | channels and graph exact; retrained scores close (tolerances below) |
| QC rerun (opt-in) | `./run_pipeline.sh --cold-qc` (either mode) | Stage `qc` always rebuilds IxIDN-Matched from the stored QC verdicts and checks it equals the shipped file; `--cold-qc` also reruns the QC with your own key (`ANTHROPIC_API_KEY`, `QC_MODEL`) into `build/qc/` using a prompt reconstructed from the Methods; drug names and trial phases from the licensed trial source are optional inputs (`--drugs`, `--phases`) | stored verdicts exact; rerun at a tolerance (`tests/test_ixidn_matched_qc.py`) |

Figures, tables, notebooks and the manuscript numbers always come from the frozen record (they need
`data/zenodo`), in both modes, because a retrained ensemble reproduces the published scores closely but not bit for bit.

Stages, in order:

| Mode | Stages |
|---|---|
| WARM | `verify graph qc figures si tests notebooks` |
| COLD | `verify channels graph qc retrain score figures si tests notebooks` |

- `verify`: checks `data/zenodo` against `data/ZENODO_SHA256SUMS.txt` (in COLD mode it is skipped with a message if `data/zenodo` is absent; comparisons with the frozen record are then skipped).
- `channels` (COLD): `python -m disatlas.channels.<broad|specific> --out build/channels/<track> --cold`.
  - specific: rebuilt from the post-filter evidence in `data/raw/channels/specific`.
  - broad: recomputes the phenotype distance `dist_HPO` from `data/raw/channels/broad/disease_phenotype.parquet`
    and `data/raw/ontology/hp.obo`, joins the shipped gene-based pair distances (`dist_GO`, `px_HPA`, `px_DOM`,
    `data/zenodo/channels_raw/broad/pair_gene_distances.parquet`), then applies the CRank/rank-mean fusion,
    class quotas, floor and repair. The result equals the shipped broad channel edge for edge (13,021 diseases,
    332,082 edges). The broad channel stage alone takes 8 min 45 s wall and 4.2 GB peak memory.
- `graph`: `python -m disatlas.graph.build` (WARM: `--check` against the shipped graph; COLD: from `build/channels`, or from `data/raw/channels` if the channel stage was not run).
- `qc`: `python -m disatlas.benchmarks.ixidn_matched --out build/ixidn_matched --check` (stored verdicts); with `--cold-qc` also `python -m disatlas.benchmarks.matched_qc --out build/qc --model $QC_MODEL`.
- `retrain` (COLD): one process per seed (4707-4711) in parallel, node2vec dim 32, walk length 20, window 2, p = 1, q = 0.5, Hadamard MLP decoder.
- `score` (COLD): `python -m disatlas.rescore --out build/scores --all-pairs`.
- `figures`, `si`: all scripts in `figures/main`, `tables/main`, then `figures/SI`, `tables/SI`.
- `tests`: `pytest tests` (in COLD mode first `tests/test_cold_exact.py` and `tests/test_cold_tolerance.py` on `build/`).
- `notebooks`: executes the four notebooks.

Start at a stage with `--from STAGE` and run one stage with `--only STAGE`, for example
`./run_pipeline.sh --only tests` or `./run_pipeline.sh --cold --from retrain`.

What COLD reproduces:

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

Measured on one cluster node (16 cores, 160 GB):

| Mode | Wall time | Peak memory |
|---|---|---|
| WARM | 14 min (tests 1 min); up to 46 min on a node shared with other jobs | 8.7 GB |
| COLD, `--from graph` | 37 min: graph 1 min; retrain 16 min (five seeds in parallel, about 2 GB each); score 7 min, mostly the all-pairs matrix; figures, tables, tests and notebooks 13 min; up to 2 h on a shared node | 10.4 GB |
| COLD, full run | 53 min on 16 CPUs: channels 10 min (broad 8 min 45 s), then as `--from graph`; measured on a fresh clone | 10.5 GB |

## Data and licensing

- `data/raw` holds minimal tables: only the columns the code reads, and only the rows that pass the published
  thresholds (no pre-filter data).
- Content from licence-restricted sources is blinded, not dropped: restricted targets and modules are replaced by
  stable salted tokens (`T_<12 hex>`, `M_<12 hex>`), so every overlap and graph computation is preserved exactly but
  the token cannot be mapped back. `disatlas.blinding.is_blinded(x)` tells you whether an identifier is a token.
- Results that depend on licensed trial and approval records are shipped as precomputed tables in `data/derived`.
- LLM-derived inputs (IxIDN-Matched curation verdicts, comparator system scores) are the authors' stored outputs,
  shipped in `data/raw/benchmarks` and used as-is; the pipeline calls no language model. The prompts are given in the
  manuscript and Supplementary Information, and you can run them with your own LLM on the blinded, shuffled benchmark
  files listed in [data/README.md, 'LLM-derived inputs'](data/README.md#llm-derived-inputs) and compare the results
  with the stored outputs.
- Column dictionaries, blinding details and the Zenodo file list: [data/README.md](data/README.md).
  Licences of the data sources: [LICENSES.md](LICENSES.md). <!-- LICENSES_MD_TBD -->

## Tests and parity

```bash
.venv/bin/python -m pytest -q tests           # or ./run_pipeline.sh --only tests (figures must be built first)
```

| Test | Checks |
|---|---|
| `tests/test_data_model.py` | blinding tokens, benchmark counts, training graph rebuilt edge for edge, ensemble reproduces the shipped IxIDN-Matched scores |
| `tests/test_main_parity.py` | every number of the main text, from `tests/main_text_numbers.tsv` and the figure/table value files |
| `tests/test_sibench_parity.py` | numbers quoted in the SI benchmark sections |
| `tests/test_sibio_parity.py` | numbers quoted in the SI biology sections (specs in `tests/parity/`) |
| `tests/test_cold_exact.py` | COLD: rebuilt channels and graph equal the shipped ones edge for edge |
| `tests/test_cold_tolerance.py` | COLD: retrained scores within the tolerances above |

A WARM run ends with `278 passed, 11 skipped` (the skipped tests are COLD-only checks) and `DisAtlas pipeline finished (warm)`.

## Troubleshooting

- `data/zenodo missing: download the Zenodo bundle`: unpack the Zenodo record into `data/zenodo` (see Quick start).
- Checksum failure in `verify`: re-download the listed file; do not edit `data/ZENODO_SHA256SUMS.txt`.
- Memory: the all-pairs score matrix (13,021 x 13,021) is expanded in memory by several scripts; WARM peaks at about
  8 GB. On smaller machines run stages one at a time with `--only`.
- Threads: the retrain stage sets `OMP_NUM_THREADS` to the number of cores divided by five; lower it on shared machines.
- `disatlas.channels.<track> is not available in this checkout`: resume with `./run_pipeline.sh --cold --from graph`,
  which builds the training graph from the shipped `data/raw/channels`.
- A figure script fails with a display error: set `MPLBACKEND=Agg` (the pipeline does this).

## How to cite

```bibtex
@article{disatlas,
  title   = {A Mechanistic Disease Representation Framework for Therapeutic Knowledge Transfer Across Rare and Common Conditions},
  author  = {AUTHORS_TBD},
  journal = {JOURNAL_TBD},
  year    = {YEAR_TBD},
  doi     = {DOI_TBD}
}
```

Data and model: Zenodo, DOI `ZENODO_DOI_TBD`.

## Licence

Code: see [LICENSE](LICENSE). <!-- CODE_LICENCE_TBD --> Data: see [LICENSES.md](LICENSES.md) and [data/README.md](data/README.md).
