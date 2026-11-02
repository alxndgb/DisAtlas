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
<!-- ONTOLOGY_TBD: data/raw/ontology is empty on main; filled by the broad-channel track -->

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

<!-- RUN_MODES_BRANCH: --cold, --from, --only, build/, disatlas.rescore and tests/test_cold_*.py are on the run-modes branch -->

| Mode | Command | What it does | Agreement with the manuscript |
|---|---|---|---|
| WARM (default) | `./run_pipeline.sh` | Uses the processed data and frozen 5-seed model from `data/zenodo`; rebuilds the training graph as a check | exact |
| COLD | `./run_pipeline.sh --cold` | Rebuilds the channels and the training graph from `data/raw`, retrains the 5-seed ensemble and re-scores the benchmarks into `build/`, then makes all figures and tables and runs the tests | channels and graph exact; retrained scores statistical (tolerances below) |

Stages, in order:

| Mode | Stages |
|---|---|
| WARM | `verify graph figures si tests notebooks` |
| COLD | `verify channels graph retrain score figures si tests notebooks` |

- `verify`: checks `data/zenodo` against `data/ZENODO_SHA256SUMS.txt` (in COLD mode it is skipped with a message if `data/zenodo` is absent).
- `channels` (COLD): `python -m disatlas.channels.<broad|specific> --out build/channels/<track> --cold`.
  The cold run rebuilds the specific channel from raw; broad channel source: see [data/README.md](data/README.md).
  <!-- BROAD_COLD_TBD -->
- `graph`: `python -m disatlas.graph.build` (WARM: `--check` against the shipped graph; COLD: from `build/channels`, or from `data/raw/channels` if the channel stage was not run).
- `retrain` (COLD): one process per seed (4707-4711), node2vec dim 32, walk length 20, window 2, p = 1, q = 0.5, Hadamard MLP decoder.
- `score` (COLD): `python -m disatlas.rescore --out build/scores --all-pairs`.
- `figures`, `si`: all scripts in `figures/main`, `tables/main`, then `figures/SI`, `tables/SI`.
- `tests`: `pytest tests` (in COLD mode also `tests/test_cold_exact.py` and `tests/test_cold_tolerance.py`).
- `notebooks`: executes the four notebooks.

Start at a stage with `--from STAGE` and run one stage with `--only STAGE`, for example
`./run_pipeline.sh --only tests` or `./run_pipeline.sh --cold --from retrain`.

Measured on a 16-core node (Linux):

| Mode | Wall time | Peak memory |
|---|---|---|
| WARM | 45 min | 8.3 GB |
| COLD, retrain stage | 15-17 min per seed (8 threads, seeds in parallel) | 1.4 GB per seed |
| COLD, full run | TBD | TBD |

COLD tolerances (`tests/test_cold_tolerance.py`). The retrain uses different random streams from the frozen run, so
scores are not bit-identical. The test asserts |AUROC(cold) - AUROC(frozen)| <= 0.03 on IxIDN-Matched (all pairs,
rare-rare, rare-common), IxIDN+/IxIDN- and the FDA development controls, <= 0.04 on the IxIDN-Matched test split, and
Spearman rho >= 0.90 between cold and frozen ensemble scores on the 4,364 IxIDN-Matched pairs and >= 0.85 on 10^6
fixed random disease pairs. Observed in one cold retrain (5-seed ensemble vs frozen): |AUROC difference| 0.003-0.022 across the six
benchmark sets above (largest on the IxIDN-Matched test split), Spearman rho 0.946 on IxIDN-Matched.

## Data and licensing

- `data/raw` holds minimal tables: only the columns the code reads, and only the rows that pass the published
  thresholds (no pre-filter data).
- Content from licence-restricted sources is blinded, not dropped: restricted targets and modules are replaced by
  stable salted tokens (`T_<12 hex>`, `M_<12 hex>`), so every overlap and graph computation is preserved exactly but
  the token cannot be mapped back. `disatlas.blinding.is_blinded(x)` tells you whether an identifier is a token.
- Results that depend on licensed trial and approval records are shipped as precomputed tables in `data/derived`.
- LLM-derived inputs (IxIDN-Matched curation verdicts, comparator system scores) are the authors' stored outputs,
  shipped in `data/raw/benchmarks` and used as-is; no language model is called and there is no option to re-run them.
  See [data/README.md](data/README.md#llm-derived-inputs), section 'LLM-derived inputs'.
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
