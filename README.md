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
node2vec ensemble embeds every disease; a pair is scored by the cosine similarity of its two disease vectors.

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

Every script runs from the repository root without arguments, except Figures S19 and S20, which are drawn by the Figure 7 script with `--si hits10` / `--si weighted` (`./run_pipeline.sh` calls them in order) and writes
its outputs next to itself. Most scripts also write a `*_values.tsv` / `*_values.json` / `*_numbers.json` file with
every plotted or quoted number; the tests read those files.

### Main text

| Item | Content | Script | Output |
|---|---|---|---|
| Figure 1 | Knowledge transfer between rare and common diseases (illustration) | none (static image) | `figures/main/fig1_knowledge_transfer.png` |
| Figure 2 | DisAtlas pipeline and its two evaluations | `figures/main/fig2_disatlas_pipeline.py` | `figures/main/disatlas_pipeline.{png,pdf}` |
| Figure 3 | The IxIDN-Matched benchmark | `figures/main/fig3_ixidn_benchmark.py` | `figures/main/ixidn_benchmark.{png,pdf}` |
| Figure 4 | Tracing a disease-disease similarity to its evidence | `figures/main/fig4_evidence_trace.py`; prose numbers `figures/main/fig4_numbers.py` | `figures/main/fig4_evidence_trace.{png,pdf}`, `figures/main/fig4_numbers.json` |
| Figure 5 | Disease map, clusters, and DisAtlas vs TxGNN by disease rarity (panel c) | `figures/main/fig5_disease_clustering_map.py` | `figures/main/disease_clustering_map.{png,pdf}` |
| Figure 6 | Cross-system comparison on IxIDN-Matched | `figures/main/fig6_method_comparison.py` | `figures/main/fig5_method_comparison.{png,pdf}` (same file name as the manuscript figure) |
| Figure 7 | Event-time test: rank of the case-study disease and conditional-logit odds ratios by year and method | `figures/main/fig7_prospective.py` | `figures/main/fig7_prospective.{png,pdf}` |
| Table 1 | Protocols of the cross-system comparison | `tables/main/table1_protocols.py` | `tables/main/table1_protocols.{md,tsv}` |
| Main-text numbers | All 568 registered number tokens of the main text, each with its kind and check | `tests/main_quantities.py`, `tests/test_main_parity.py` | registry `tests/main_text_numbers.tsv` |

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
| Figure S16 | Performance by disease rarity and anchor disease; TxGNN recovery of documented rare-common partners | `si_rare_common_transfer.py` | `SI_rare_common_transfer.{png,pdf}` |
| Figure S17 | Rare-common neighbour coverage | `si_rare_common_neighbour_coverage.py` | `SI_rare_common_neighbour_coverage.{png,pdf}` |
| Figure S18 | Event-time holdout: per-drug probability by method and decision year (matched analysis) | `si_event_time_auc_methods.py` | `SI_event_time_auc_methods.{png,pdf}`, `SI_event_time_auc_methods_values.json` |
| Figure S19 | Event-time holdout: top-10% share by decision year, with and without matching | `../main/fig7_prospective.py --si hits10` | `SI_event_time_hits10.{png,pdf}` |
| Figure S20 | Event-time holdout: top-10% share weighted for disease popularity (post hoc) | `../main/fig7_prospective.py --si weighted` | `SI_event_time_weighted.{png,pdf}` |
| Figure S21 | Later-trialled diseases within popularity strata | `si_popularity_strata.py` | `si_popularity_strata.{png,pdf}` |
| Figure S22 | Approved-repurposing benchmark (schematic) | `si_fda_schematic.py` | `SI_approved_repurposing_schematic.{png,pdf}` |
| Figure S23 | Evidence channels on FDA-approved repurposing pairs | `si_fda_channels.py` | `SI_fda_channels.{png,pdf}` |
| Figure S24 | Mechanism coverage of approved expansions in dated graph slices | `si_mechanism_coverage_slices.py` | `SI_mechanism_coverage_slices.{png,pdf}` |
| Figure S25 | Clinical milestones of complement inhibitors | `si_complement_milestones.py` | `SI_complement_milestones.{png,pdf}` |
| Figure S26 | Evidence behind three rare-common case pairs | `si_bio_gap_numbers.py` (case scores) | `SI_explainability.png` (static), `SI_explainability_case_scores.tsv` |
| Figure S27 | How the training graph is built and scored | `si_graph_channels_schematic.py` | `SI_graph_channels_schematic.{png,pdf}` |
| Figure S28 | Test of compositionality | `si_compositionality_recompute.py` (tables), `si_compositionality.py` (figure) | `SI_compositionality.{png,pdf}` |

### Supplementary tables


Evidence-trace and explainability items (final-manuscript numbering; rows above use the earlier numbering and are renumbered at integration):

| Item | Content | Script | Output |
|---|---|---|---|
| Figure S4 | Lowest common ancestor in the disease ontology (illustration) | none (static image) | `figures/SI/LCA.png` |
| Figure S6 | Embedding pair score on IxIDN-Matched pairs | `figures/SI/si_mechdist_embedding_pair_score.py` | `figures/SI/SI_mechdist_embedding_pair_score.png`, `_values.json` |
| Figure S27 | Explainability case studies (cosine scores; TxGNN panels static) | `figures/SI/si_explainability.py` | `figures/SI/SI_explainability.png`, `_values.json` |
| Table S22 | Figure 4 items changed by specificity weighting | `tables/SI/table_SI_trace_order.py` | `tables/SI/table_SI_trace_order.{md,tsv}` |
| Table S23 | Figure 4 displayed items with walk-mass share and specificity weight | `tables/SI/table_SI_trace_order.py` | `tables/SI/table_SI_trace_order_values.{md,tsv}` |

| Item | Content | Script | Output |
|---|---|---|---|
| Table S1 | Disease-gene association sources | none (descriptive) | in the SI text |
| Table S2 | Disease-symptom association sources | none (descriptive) | in the SI text |
| Table S3 | Trials behind Figure 3a | `tables/SI/table_fig3a_trials.py` | `tables/SI/table_fig3a_trials.{md,tsv}` |
| Table S4 | Ontology-only baselines | `tables/SI/table_ontology_baseline.py` | `tables/SI/table_ontology_baseline.md` |
| Table S5 | Node-disjoint DisAtlas vs language-model systems by subgroup | `tables/SI/table_llm_subgroups.py` | `tables/SI/table_llm_subgroups.{md,tsv}` |
| Table S6 | Degree-matched and within-disease AUROC | `tables/SI/table_degree_controls.py` | `tables/SI/table_degree_controls.{md,tsv}` |
| Table S7 | Positive and negative control pairs | `tables/SI/table_cluster_controls.py` | `tables/SI/table_cluster_controls.{md,tsv}` |
| Table S8 | Theme of each of the 55 disease clusters | `tables/SI/table_cluster_themes.py` | `tables/SI/table_cluster_themes.{md,tsv}`, `table_cluster_themes_numbers.json` |
| Table S9 | AUROC by disease-rarity pairing | `tables/SI/table_rare_common_strata.py` | `tables/SI/table_rare_common_strata.{md,tsv}` |
| Tables S10, S11 | Event-time holdout: pooled advantages; subsets of new trials | `figures/SI/si_event_time_tables.py` | `figures/SI/SI_event_time_tables.md` |
| Tables S12–S17 | Event-time: NMO rank by slice; yearly odds ratios; AUROC (pooled and yearly); early decision years; small comparison sets | `figures/SI/si_event_time_or_tables.py` | `figures/SI/SI_event_time_or_tables.md`, `SI_event_time_or_tables_numbers.json` |
| Table S18 | Advantage with current and slice-dated GO annotations | `figures/SI/si_go_release_table.py` | `figures/SI/SI_go_release_table.{md,tsv}` |
| Table S19 | FDA-approved pairs scored with dated models | `figures/SI/si_fda_prioritisation_text.py` | `figures/SI/SI_fda_prioritisation_text.md` |
| SI numbers | Numbers quoted in the SI text | `figures/SI/si_*_numbers.py` | `figures/SI/SI_*_numbers.json`; asserted by `tests/test_sibench_parity.py`, `tests/test_sibio_parity.py` (specs in `tests/parity/`) |

## Repository structure

```
run_pipeline.sh          one command for figures, tables, tests and notebooks (see Run modes)
requirements.txt         pinned Python packages
data/
  raw/                   inputs of the graph build: public disease-gene sources, Open Targets 26.09, HPO,
                         the 8,343 released knowledge-graph disease-gene edges, broad-channel gene tables and the
                         released specific-channel edges (channels/specific/); see data/raw/README.md
  raw/benchmarks/        IxIDN+/IxIDN- pairs, IxIDN-Matched pairs and curation verdicts, FDA controls, comparator outputs
  derived/               small precomputed result tables read by the figure, table and test scripts
  zenodo/                (download, not in git) broad graph, term-distance slices, frozen model and its outputs,
                         score matrices, large tables
  README.md              column dictionaries, licence and blinding
  ZENODO_SHA256SUMS.txt  checksums of data/zenodo (model outputs: ZENODO_MODEL_OUTPUTS_SHA256SUMS.txt)
src/disatlas/            library: graph build, model loading and scoring, retraining, blinding, paths, figure style
figures/main, figures/SI scripts for the main and Supplementary figures, and their outputs
tables/main, tables/SI   scripts for Table 1 and the Supplementary tables, and their outputs
notebooks/               four exploratory notebooks (data and graph, IxIDN benchmarks, SI robustness, SI biology)
tests/                   parity tests and the main-text number registry
scripts/                 parity replays, Habcd retrain template, notebook runner; scripts/maintainer/ shows how the
                         shipped data were exported (not needed to reproduce)
docs/                    model outputs, optional retrains, errata
release/                 data-minimisation and knowledge-graph exposure ledgers
```

### Training graph

The training graph (P_full: 532,958 edges, 13,042 diseases) is the union of two disease-disease channels:

- **broad** (369,405 edges, shipped as `data/zenodo/graph/graph_edges.parquet`): its disease-gene table of record
  has 24,450 tier-1 rows, 604 of them licence-restricted pairs that are not redistributed (see
  [Data and licensing](#data-and-licensing)). `python -m disatlas.graph.raw_build` rebuilds the broad channel from
  `data/raw`, whose disease-gene table is 21,066 rows from public sources (GenCC, Orphanet, Open Targets 26.09
  direct) plus the 8,343 released knowledge-graph disease-gene edges, and nothing else. Of its 23,840 tier-1 rows,
  16,975 are public and 6,865 come from the released edges; the rebuilt broad graph has 364,288 edges. The GO and HPO term-distance slices it reads ship in
  `data/zenodo/specific_channels/`.
- **specific** (189,903 edges): shipped as disease-disease edges with their weight only
  (`data/raw/channels/specific/specific_edges.parquet`). Rebuilding it from its inputs would need about 1.06 million
  rows of licence-restricted knowledge-graph evidence (literature, pathology, MetaBase, CTD and mechanism nodes),
  which is why only its edges are released; a rebuild from public sources only is FAR from the released model
  (test AUROC -0.066).

`python -m disatlas.graph.pfull_build` combines them (outer union, FDA-pair strip, walk weight `w_arm`,
gene-anchored step) and reproduces the recorded P_full graph exactly (edge set, row order and `w_arm`,
`tests/test_pfull_graph.py`). The released disease-gene edges alone are not sufficient: the training graph also
needs the specific-channel edges.

## Quick start

Requirements: Linux or macOS, Python 3.11, about 16 GB of RAM, and about 1.2 GB of disk for the unpacked Zenodo bundle.

```bash
git clone <repository URL> DisAtlas && cd DisAtlas
python3.11 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

`requirements.txt` pins the package versions used to render the figures; with them a fresh run regenerates every committed
figure byte for byte (other matplotlib versions can shift anti-aliasing by a few grey levels without changing any number).

Download the Zenodo record (DOI: `ZENODO_DOI_TBD`, <https://doi.org/ZENODO_DOI_TBD>) and unpack it so that
`data/zenodo/` contains `graph/`, `specific_channels/`, `model/`, `model_outputs/`, `scores/`, `figures_main/`, `si_bio/`. Then check it:

```bash
cd data/zenodo && sha256sum -c ../ZENODO_SHA256SUMS.txt && sha256sum -c ../ZENODO_MODEL_OUTPUTS_SHA256SUMS.txt && cd ../..   # macOS: shasum -a 256 -c
./run_pipeline.sh
```

`run_pipeline.sh` uses `.venv/bin/python` (override with `PYTHON=/path/to/python`) and sets `PYTHONPATH=src`.
To run a single script by hand: `export PYTHONPATH=$PWD/src MPLBACKEND=Agg` and then, for example,
`.venv/bin/python figures/main/fig3_ixidn_benchmark.py`.

## Run modes

| Mode | Command | What it does | Agreement with the manuscript |
|---|---|---|---|
| WARM (default) | `./run_pipeline.sh` | Uses the processed data and frozen model outputs from `data/zenodo`; rebuilds the P_full training graph and IxIDN-Matched as checks | exact |
| COLD | `./run_pipeline.sh --cold` | Rebuilds the P_full training graph into `build/` (with `DISATLAS_RAWGRAPH=1` also the broad graph from `data/raw`), retrains the 5-seed ensemble and re-scores the benchmarks, checks the graph against its recorded signature and the retrain with the pre-registered parity test, then makes all figures and tables and runs the tests | graph exact; retrain judged by the parity test |
| QC rerun (opt-in) | `./run_pipeline.sh --cold-qc` (either mode) | Stage `qc` always rebuilds IxIDN-Matched from the stored QC verdicts and checks it equals the shipped file; `--cold-qc` also reruns the QC with your own key (`ANTHROPIC_API_KEY`, `QC_MODEL`) into `build/qc/` using a prompt reconstructed from the Methods; drug names and trial phases from the licensed trial source are optional inputs (`--drugs`, `--phases`) | stored verdicts exact; rerun at a tolerance (`tests/test_ixidn_matched_qc.py`) |

Figures, tables, notebooks and the manuscript numbers always come from the frozen record (they need
`data/zenodo`), in both modes, because a retrained ensemble reproduces the published scores closely but not bit for bit.

Stages, in order:

| Mode | Stages |
|---|---|
| WARM | `verify graph qc figures si tests notebooks` |
| COLD | `verify graph qc retrain score figures si tests notebooks` |

- `verify`: checks `data/zenodo` against `data/ZENODO_SHA256SUMS.txt` (in COLD mode it is skipped with a message if `data/zenodo` is absent; comparisons with the frozen record are then skipped).
- `graph`: `python -m disatlas.graph.pfull_build --out build/graph` (broad graph from `data/zenodo/graph`; in COLD
  mode with `DISATLAS_RAWGRAPH=1`, first `python -m disatlas.graph.raw_build --out build/rawgraph`, about 1 h, and
  the P_full graph is built from that). Both builds check their output against the recorded signatures.
- `qc`: `python -m disatlas.benchmarks.ixidn_matched --out build/ixidn_matched --check` (stored verdicts); with `--cold-qc` also `python -m disatlas.benchmarks.matched_qc --out build/qc --model $QC_MODEL`.
- `retrain` (COLD): `python -m disatlas.retrain --graph build/graph/pfull_graph_edges.parquet`, one process per seed (4707-4711) in parallel; weighted node2vec on `w_arm`, dim 32, walk length 20, window 2, p = 1, q = 0.5, Hadamard MLP decoder.
- `score` (COLD): `python -m disatlas.rescore --out build/scores --all-pairs`.
- `figures`, `si`: all scripts in `figures/main`, `tables/main`, then `figures/SI`, `tables/SI`.
- `tests`: `pytest tests` (in COLD mode first `tests/test_cold_exact.py` and `tests/test_cold_tolerance.py` on `build/`).
- `notebooks`: executes the four notebooks.

Start at a stage with `--from STAGE` and run one stage with `--only STAGE`, for example
`./run_pipeline.sh --only tests` or `./run_pipeline.sh --cold --from retrain`.

What COLD reproduces:

- **Exact** (`tests/test_cold_exact.py`, `tests/test_pfull_graph.py`, `tests/test_raw_graph.py`): the P_full
  graph rebuilt from the released inputs equals its recorded signature (`tests/pfull_rebuild_signature.json`:
  532,929 edges with the shipped broad graph, 528,008 with `DISATLAS_RAWGRAPH=1`; edge-set, row-order and `w_arm`
  SHA-256); the shipped broad graph (369,405 edges) and the signature of the released model's P_full (532,958
  edges, 13,042 nodes) are checked unchanged.
- **Parity** (`tests/test_cold_tolerance.py`): the retrained ensemble is compared with the released model by the
  pre-registered rule (5-seed mean endpoint cosine on the 1,094 IxIDN-Matched test pairs; 95 % disease-cluster
  bootstrap CI of the AUROC difference, B = 2,000; PARITY if the lower bound is > -0.02 and coverage >= 95 %).
  The authors' 5-seed retrain on P_full gives test AUROC 0.8589 against 0.8605 for the released model,
  ΔAUROC -0.002 [-0.008, +0.005], coverage 100 %: PARITY. `python scripts/pfull_parity.py` replays this from
  the shipped cosines (`tests/test_pfull_parity.py`). Walk sampling and floating-point sums depend on the
  thread count and BLAS build, so a retrain is judged by this rule, not by identity.
- The warm test suite then runs as in a warm run, so the manuscript numbers stay exact.

Measured on fresh clones, one cluster node each (peak memory is the Slurm MaxRSS of the job):

| Mode | Wall time | CPUs (memory allocated) | Peak memory |
|---|---|---|---|
| WARM | WARM_TIME_TBD | 16 (64 GB) | WARM_MEM_TBD |

## Data and licensing

Of the 918 curated disease–gene pairs used for the graph of record, the 604 supported by OMIM alone are not redistributed for licence reasons, and 6 more with no public source are not released. The 308 released pairs (`data/raw/public_disease_gene/`) are all in the public Open Targets 26.09 direct evidence, below the thresholds used to build the shipped Open Targets files. The remaining released BIKG-derived data are shared with permission from the AstraZeneca data office. All reported results reproduce exactly from the shipped graph files and model outputs. Rebuilding the gene-anchored graph from the public inputs here changes 36,868 edges (6.9%) relative to the training graph of the released model; a 5-seed retrain on a rebuild that differs from this one by 314 edges (0.06%) is within retraining noise (IxIDN-Matched test ΔAUROC −0.0016, 95% CI −0.0074 to +0.0041). The note "rebuilt from shipped inputs" in the shipped `model_config.json` predates the exclusion of the OMIM-only pairs; the final graph reproduces from the shipped graph files.

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
| `tests/test_data_model.py` | blinding tokens, benchmark counts, broad-graph build statistics, ensemble reproduces the shipped IxIDN-Matched scores |
| `tests/test_pfull_graph.py` | P_full training graph rebuilt from the repository equals the recorded graph (edge set, row order, `w_arm`) |
| `tests/test_raw_graph.py` | broad graph rebuilt from `data/raw` equals the shipped one (`DISATLAS_RAWGRAPH=1`, about 1 h) |
| `tests/test_pfull_parity.py` | the recorded P_full retrain parity result replays from the shipped cosines |
| `tests/test_main_parity.py` | every number of the main text, from `tests/main_text_numbers.tsv` and the figure/table value files |
| `tests/test_sibench_parity.py` | numbers quoted in the SI benchmark sections |
| `tests/test_sibio_parity.py` | numbers quoted in the SI biology sections (specs in `tests/parity/`) |
| `tests/test_cold_exact.py` | COLD: the P_full graph in `build/graph` equals the recorded one |
| `tests/test_cold_tolerance.py` | COLD: the retrained ensemble passes the parity test |

A WARM run ends with `WARM_TESTS_TBD` and `DisAtlas pipeline finished (warm)`.

## Troubleshooting

- `data/zenodo missing: download the Zenodo bundle`: unpack the Zenodo record into `data/zenodo` (see Quick start).
- Checksum failure in `verify`: re-download the listed file; do not edit `data/ZENODO_SHA256SUMS.txt`.
- Memory: the all-pairs score matrix (13,021 x 13,021) is expanded in memory by several scripts; WARM peaks at about
  13 GB (measured 12.3-12.6 GiB in fresh-clone runs). On smaller machines run stages one at a time with `--only`.
- Threads: the retrain stage sets `OMP_NUM_THREADS` to the number of cores divided by five; lower it on shared machines.
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
