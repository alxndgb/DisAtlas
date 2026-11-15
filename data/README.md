# DisAtlas data

The small tables live in this repository under `data/raw` and `data/derived`. The broad graph, the
term-distance slices, the frozen model and its outputs and the large score files are in the Zenodo record and unpack into `data/zenodo`, which is
not tracked by git. A warm run (default) reads them; a cold run rebuilds them into `build/` (`README.md`). `data/ZENODO_SHA256SUMS.txt`
lists every Zenodo file with its checksum (`cd data/zenodo && sha256sum -c ../ZENODO_SHA256SUMS.txt`); `release/ZENODO_MANIFEST.tsv` adds sizes
and licences (`python scripts/maintainer/zenodo_manifest.py --check`). `data/ZENODO_MODEL_OUTPUTS_SHA256SUMS.txt` (repository-relative paths) is the subset for `model_outputs/`.

## Folder map

```
data/
  raw/              inputs of the training-graph build (column dictionaries and row counts: raw/README.md)
    gencc/, orphanet/, public_disease_gene/, opentargets/   public disease-gene edges (Open Targets 26.09 direct)
    bikg_habcd/       the 8,343 released knowledge-graph disease-gene edges
    hpo/              disease-phenotype edges
    broad_channels/   HPA cell-type profiles and Pfam domains of disease genes
    specific_channels/ GO annotations of disease genes (gene2go.parquet) and the GO/HPO term-distance parameters; inputs of the broad channel. The folder name is historical: it holds no specific-channel data (the specific channel is released as edges only, in channels/specific/).
    channels/specific/ the released specific-channel edges (u, v, w)
    mondo/, go/       MONDO and GO extracts read by the cluster-theme and MONDO-concordance rebuilds
    training_graph/   training graph of the v3.2 embedding (benchmark and FDA pairs stripped)
    benchmarks/
      ixidn/            IxIDN+ and IxIDN- pair lists (removed from the training graph)
      ixidn_matched/    IxIDN-Matched evaluation pairs (4,364) and the language-model curation verdicts
      fda_repurposing/  FDA-approved repurposing positive controls (holdout split withheld)
      comparators/      blinded comparator input and stored per-run comparator outputs on IxIDN-Matched
      disease_labels.tsv
  derived/          intermediate and result tables behind the figures, tables and quoted numbers (one folder per item)
  zenodo/           (download) graph/, ontology_term_term_distance/, model/, model_outputs/, figures_main/, si_bio/, event_time/
                    (warm run; a cold run rebuilds the graph and the model into build/)
```

## Licence and blinding

Open Targets identifiers are released unchanged. Targets and modules that come from licence-restricted
sources (curated pathway databases, the literature-mining knowledge graph and similar) cannot be
redistributed. They are replaced by stable salted tokens: `T_<12 hex>` for targets and `M_<12 hex>` for
modules. The same source identifier always maps to the same token, so overlaps between diseases and
all graph computations are preserved exactly, but the token cannot be mapped back without the private
salt. For these rows the source is given only at group level (`source_group`, `families_group`). No
publication identifiers or evidence sentences are included. Disease identifiers are public: MONDO, or Orphanet where no MONDO term exists, and for 149 diseases EFO,
MeSH or OMIM; every file uses the same identifier for a disease.
`disatlas.blinding.is_blinded(x)` tells you whether an identifier is a token.

Trial identifiers (ClinicalTrials.gov NCT ids) and trial dates are released as shown in the paper's tables and figures
(author decision 2026-11-15), for example Table S3, Table S14, Figure 7 and Figures S25-S26 and the derived tables they read.

## LLM-derived inputs

Four inputs come from language-model or agent runs made by the authors. The pipeline uses these stored
outputs as they are and makes no language-model calls by default. Where a prompt is printed in the
manuscript, users who want to can run it with their own LLM on the input file below and compare the result
with the stored output. For the IxIDN-Matched curation the opt-in `./run_pipeline.sh --cold-qc` reruns
the QC with the user's own key from a prompt reconstructed from the Methods (`disatlas.benchmarks.matched_qc`,
README "Run modes"); no other rerun code is shipped.

| table | what the model produced | prompt in the manuscript | input file | stored output (columns) |
|---|---|---|---|---|
| comparator scores | a pair score for each IxIDN-Matched pair, per run (Claude Sci, Claude Chat, Gemini Pro; Phylo Biomni agent seeds) | Methods, "Comparator systems" (sec:methods_comparators), Protocol 2 and Claude Sci; Phylo Biomni (Protocol 3) is described there without a prompt | `raw/benchmarks/comparators/ixidn_matched_blind_input.tsv` (blind_id, a, b, disease_a_name, disease_b_name; rows shuffled, no labels) | `raw/benchmarks/comparators/ixidn_matched_comparator_scores.tsv` (pair_id = blind_id, one column per `system\|run`); read via `derived/fig6_comparison/per_run_scores.tsv` |
| IxIDN-Matched curation verdicts | STRONG / PLAUSIBLE / WEAK / SUSPECT / NOT_IN_GRAPH for each of the 4,302 candidate positives | not printed; described in Methods, "IxIDN-Matched benchmark" (sec:methods_ixidn_matched) and SI, "How IxIDN-Matched is built" (sec:SI_ixidn_matched_build) | none shipped (the shared-drug evidence comes from licence-restricted trial data) | `raw/benchmarks/ixidn_matched/llm_curation_verdicts.tsv` (disease_a, disease_b, verdict, rejudged, first_pass_verdict, keep); read via `derived/fig3_ixidn/matched_candidate_verdicts.tsv` (disease_a, disease_b, verdict) |
| disease labels | rare / common and cancer flags for diseases whose rule signals were missing or conflicting (literature-tool agent) | not printed; described in Methods, "Disease set" (sec:methods_disease_set) | none shipped (the agent used literature lookups) | `raw/benchmarks/disease_labels.tsv` (disease_id, is_rare, is_cancer, is_malignant, disease_class) and `derived/fig5_clustering/disease_map_k55.tsv.gz` (disease_class, is_rare) |
| approved indications per label version | the diseases named in each version of a drug's public US label, used for first-label dates | not printed; SI, "Mechanism coverage of approved indication expansions" (sec:SI_mechanism_coverage) | none shipped (public DailyMed label versions) | `derived/si_mechanism_timelines/mechanism_timeline_approved_indications.tsv` (disease_id, first_label_date, first_label_setid, first_label_version and the milestone columns) |

## Column dictionaries

### raw (graph inputs)

Column dictionaries, row counts and the assembly rule are in [raw/README.md](raw/README.md). The disease-gene
table assembled from them (`disatlas.graph.raw_build.gene_edges()`) is 21,066 public rows plus the 8,343 released
knowledge-graph rows, and nothing else. Of its 23,840 tier-1 rows, 16,975 are public and 6,865 come from the
released edges. Licence-restricted disease-gene pairs (604 tier-1 rows of the table of record) are not redistributed ([raw/README.md](raw/README.md)).

### raw/channels/specific

| file | columns |
|---|---|
| specific_edges.parquet | u, v (u < v; MONDO, or for 49 diseases a native EFO/Orphanet/MeSH identifier), w: specific-channel term of the walk weight after mechanism-class balancing; 189,903 edges |

Only the disease-disease edges of the specific channel are released. A rebuild from its inputs needs about
1.06 million rows of licence-restricted knowledge-graph evidence (literature, pathology, MetaBase, CTD and
mechanism nodes), and a rebuild from public sources only loses 0.066 test AUROC.

### raw/benchmarks

| file | columns |
|---|---|
| ixidn/ixidn_plus_pairs.parquet, ixidn_minus_pairs.parquet | disease_a, disease_b |
| ixidn_matched/ixidn_matched_pairs.tsv | pair_id, disease_a, disease_b, y (1 positive / 0 negative), split (tune / test / cross: both diseases in the tune half, both in the test half, or one in each), pair_class (rare/common composition) |
| ixidn_matched/llm_curation_verdicts.tsv | disease_a, disease_b, verdict, rejudged, first_pass_verdict, keep (kept positives = IxIDN-Matched positives): the stored language-model QC of the 4,302 candidate positives; NOT_IN_GRAPH = no shared trial drug, assigned without a call; first_pass_verdict = earlier verdict of the 944 re-judged pairs (961 filled) |
| ixidn_matched_qc/matching_pairs.parquet | disease_a, disease_b, y, cell, pair_id: veto-passing IxIDN- pairs (y=0) and the QC candidate positives (y=1) in the frozen row order of the seeded (seed 0) cell match; cell = seven matching keys joined by '\|' (lower/higher disease degree decile, shared organ system, shared Orphanet area, name-similarity bin, semantic-similarity bin, gene coverage); pair_id = benchmark pair id or empty. `python -m disatlas.benchmarks.ixidn_matched` rebuilds ixidn_matched_pairs.tsv from it and the verdicts |
| fda_repurposing/fda_repurposing_pairs.tsv | disease_a, disease_b: all 2,000 control pairs, unlabelled, removed from the training graph (needed to rebuild it exactly) |
| fda_repurposing/fda_repurposing_controls.tsv | disease_a, disease_b, stratum (hard / easy), pair_class, split (dev / easy_unsplit / complement_six), fda_verified_pair, in_model; the 560 holdout pairs are sealed and not included |
| comparators/ixidn_matched_blind_input.tsv | blind_id (= pair_id), a, b, disease_a_name, disease_b_name: the blinded, row-shuffled file given to the comparator systems |
| comparators/ixidn_matched_comparator_scores.tsv | pair_id, then one column per `system\|run` |
| disease_labels.tsv | disease_id, is_rare, is_cancer, is_malignant, disease_class |

Column provenance of `comparators/ixidn_matched_comparator_scores.tsv`. Each column is the `DDA` column of the named stored per-run file, an exact match on blind_id (= pair_id), checked on 2026-11-14 and pinned by per-column checksum in `tests/test_comparator_scores.py`. The original files are not shipped. NA means the system returned no score; analyses rank such pairs last (Figure 6a: one value below the lowest scored pair; `data/derived/fig5_clustering/fig5c_txgnn_pair_scores.tsv` stores TxGNN's 213 NA pairs as that sentinel, -9.929e-10).

| column | original file | NA pairs |
|---|---|---|
| `Claude-Sci Opus 5.5\|Claude-Sci-Opus-5-5_RUN1` | `Claude-Sci-Opus-5-5_RUN1_ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Claude-Sci Opus 5.5\|Claude-Sci-Opus-5-5_RUN2` | `Claude-Sci-Opus-5-5_RUN2_ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Claude-Sci Opus 5.5\|Claude-Sci-Opus-5-5_RUN3` | `Claude-Sci-Opus-5-5_RUN3_ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Claude-Sci Opus 5.5\|Claude-Sci-Opus-5-5_RUN4` | `Claude-Sci-Opus-5-5_RUN4_ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Claude-Sci Opus 5.5\|Claude-Sci-Opus-5-5_RUN5` | `Claude-Sci-Opus-5-5_RUN5_ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Claude-Chat Opus 5.5\|Claude-Chat-Opus-5-5` | `Claude-Chat-Opus-5-5_RUN1_ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Claude-Chat Opus 5.5\|Claude-Chat-Opus-5-5_RUN2` | `Claude-Chat-Opus-5-5_RUN2_ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Claude-Chat Opus 5.5\|Claude-Chat-Opus-5-5_RUN3` | `Claude-Chat-Opus-5-5_RUN3_ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Claude-Chat Opus 5.5\|Claude-Chat-Opus-5-5_RUN4` | `Claude-Chat-Opus-5-5_RUN4_ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Claude-Chat Opus 5.5\|Claude-Chat-Opus-5-5_RUN5` | `Claude-Chat-Opus-5-5_RUN5_ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Gemini Pro 3.1\|Gemini-Pro-3-1-RUN1` | `Gemini-Pro-3-1-RUN1-ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Gemini Pro 3.1\|Gemini-Pro-3-1-RUN2` | `Gemini-Pro-3-1-RUN2-ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `Gemini Pro 3.1\|Gemini-Pro-3-1-RUN3` | `Gemini-Pro-3-1-RUN3-ixidn_expanded_BLIND_with_DDA.tsv` | 0 |
| `TxGNN published ckpt\|DistMult` | `txgnn_pair_scores.tsv.gz` (column `score`; NA where `mapped` is false (187) or the pair `collapsed` (26)) | 213 |
| `Phylo\|Phylo-Deepwalk-DDA-seed42` | `Phylo-Deepwalk-DDA-ixidn_expanded_BLIND_with_DDA_seed42.tsv` | 291 |
| `Phylo\|Phylo-Deepwalk-DDA-seed43` | `Phylo-Deepwalk-DDA-ixidn_expanded_BLIND_with_DDA_seed43.tsv` | 291 |
| `Phylo\|Phylo-Deepwalk-DDA-seed44` | `Phylo-Deepwalk-DDA-ixidn_expanded_BLIND_with_DDA_seed44.tsv` | 291 |
| `Phylo\|Phylo-Deepwalk-DDA-seed45` | `Phylo-Deepwalk-DDA-ixidn_expanded_BLIND_with_DDA_seed45.tsv` | 291 |
| `Phylo\|Phylo-Deepwalk-DDA-seed46` | `Phylo-Deepwalk-DDA-ixidn_expanded_BLIND_with_DDA_seed46.tsv` | 291 |

TxGNN NA pairs (213, all kept as unmappable). 187 pairs have `mapped` false: at least one of the two diseases has no TxGNN disease node. 26 pairs are `collapsed`: both diseases map to the same TxGNN disease node (idx_a = idx_b in the original file), so the published checkpoint can only return a self-score for them (cosine 1). The original file does carry a DistMult value for these 26, but it does not compare two diseases, so the merged table sets them to NA and every analysis ranks them with the other unmappable pairs (last). Methods and Results report the 213 together as the pairs the checkpoint cannot map (ranked last); the 4,151 scored pairs are the "can map" subset.


### zenodo/graph

| file | columns |
|---|---|
| graph_edges.parquet | u, v: the broad graph of record (369,405 edges, 13,042 diseases); `python -m disatlas.graph.raw_build` rebuilds it from `data/raw` without the licence-restricted pairs (364,288 edges; `DISATLAS_RAWGRAPH=1 pytest tests/test_raw_graph.py`) |
| candidate_pairs.parquet, graph_signature.json, graph_build_stats.json | candidate pairs of the broad build, the recorded signature and build counts |

The training graph of the v3.2 embedding is shipped as `raw/training_graph/training_graph_edges.parquet` (benchmark and FDA pairs
stripped; `raw/README.md`). The graph of record (532,958 edges, 13,042 diseases; signature in
`tests/pfull_graph_signature.json`) used the licence-restricted pairs. `python -m disatlas.graph.pfull_build --out build/graph`
builds it from this file, `raw/channels/specific/specific_edges.parquet` and the released disease-gene table
(532,929 edges, 13,042 diseases; columns u, v, in_broad, in_spec, w_arm) and checks it against
`tests/pfull_rebuild_signature.json`.

### zenodo/ontology_term_term_distance

GO BP (semi-structural), MF and CC and HPO term-term distance slices (int16, 89.9 MB in four files), read by the broad-graph build (`raw_build`). They are not read by the cluster-theme rebuild. Named `specific_channels/` until 2026-11-15; it holds no specific-channel data (the released specific-channel edges are `data/raw/channels/specific/`).

### zenodo/model_outputs

The outputs of the final (gene-anchored) model and of the P_full parity run: per-seed embeddings, IxIDN-Matched
cosine scores, cosine matrices and `pfull_reference/`. Contents and source roles:
`model_outputs/README.md`, `MANIFEST.tsv` and [docs/MODEL_OUTPUTS.md](../docs/MODEL_OUTPUTS.md).

### zenodo/event_time

Per-event inputs of the event-time rebuild (`src/disatlas/rebuild/event_time.py`, `tests/test_rebuild_event_time.py`, skipped when
absent), staged by `scripts/maintainer/stage_event_time.py`. These files carry only the first-trial (drug, disease, year)
events that the rebuild reads, without trial ids. Licence: CC BY-SA 3.0, inheriting the ChEMBL share-alike terms (author decision 2026-11-15; author to confirm).

| file | rows | columns |
|---|---|---|
| event_percentiles.parquet | 1,541,932 | pool, y (decision year), compound_id, event_v3 (disease), method, score, pct_p1, p1_merge_level, p1_small, pct_s1, pct_s2, sd_s2, s2_widen, s2_short, pct_s3, pct_s4, n_cand, pct_s5, n_events_drug_year, covered |
| event_top10_precision.parquet | 1,501,920 | pool, compound_id, method, pop_stratum, prec10 (per-drug top-10% precision), y |

Table S13 (`derived/si_event_time/table_pooled_top10_precision.tsv`) is rebuilt from `event_top10_precision.parquet`
(gene-anchored models, cosine readout, as for every other event-time number; `tests/test_rebuild_event_time.py`).

### zenodo/model/disatlas_v3.2_ensemble

`embeddings/node2vec_s<seed>.npz` (ids, X: unit-norm 32-d embeddings), `nodes.txt`, `model_config.json` (recipe and
SHA-256 of every file). Seeds 4707-4711. Load with `disatlas.model.load_ensemble("data/zenodo/model/disatlas_v3.2_ensemble")`;
`score_pairs(A, B)` returns the 5-seed mean cosine and `top_partners(d, k)` the highest-cosine partners.

This is the earlier v3.2 ensemble (13,021 diseases, trained on an earlier 498,155-edge graph that is not shipped;
`model_config.json` records its edge-set hash, and under `final_model_training_graph` points to the P_full build
output `build/graph/pfull_graph_edges.parquet` with its recorded signature). It is kept as the
reference layout of a retrain and for the earlier-graph comparison of `tests/test_pfull_graph.py`.
The final model of the manuscript, and the default of `load_ensemble()`, is the gene-anchored ensemble trained on P_full; its
embeddings and cosine scores are in `zenodo/model_outputs`, and the P_full retrain is compared with it. DisAtlas has no
decoder: every pair score is the 5-seed mean cosine of the two embeddings, and the clustering distance is 1 - cosine.

## Table S10 (cluster themes) and the disease-gene table

Table S10 gene-share, GO-process and HPO overlap folds are descriptive overlaps of each disease's gene, GO
biological-process and HPO annotations within a cluster against random pairs; they use the annotation tables that are
inputs of the broad channel (disease-gene, gene2go, disease-HPO), not channel edges or term-term distances. The authors'
disease-gene table holds 610 disease-gene pairs that are not redistributed (604 supported by OMIM alone, 6 with no public
source; see [raw/README.md](raw/README.md)). They are inputs of the broad channel and of the gene-anchored step of the
training graph; they are not edges of either channel. The authors' table also lacks 4,047 edges that the released table
(29,409 rows) carries (2,193 relaxed-tier and 1,854 tier-1 Open Targets edges), and the released HPO table has 202,778 rows
against 244,394 (no propagated rows, 1,722 fewer tier-1 rows).

Without the 610 pairs the gene-share contrast (anatomy- versus process-defined clusters, two-sided Mann-Whitney) is
p = 0.052 instead of the printed 0.032: omitting them alone, with every other input as in the authors' tables, gives
p = 0.052, and the other input differences do not move it across 0.05. With all released inputs the medians are 6.1 and
13.1 (printed 7.8 and 17.2), and the HPO and GO-process contrasts are p = 0.004 and p = 0.002 (printed 0.005 and 0.001),
unchanged in conclusion. The GO BP rank-1 term matches the stored table in 42 of 55 clusters (55 with the authors'
tables); this drop comes from the three disease-gene table differences together (each alone gives 50, 47 or 48), not from
the 610 pairs alone. No cluster's theme changes. The printed values are those of the authors' tables;
`python -m disatlas.rebuild.cluster_themes` rebuilds the released-input values.

## Retraining (cold run)

`./run_pipeline.sh --cold` builds the P_full training graph into `build/graph` (with `DISATLAS_RAWGRAPH=1` the
broad graph is first rebuilt from `data/raw`), retrains the 5-seed ensemble (`python -m disatlas.retrain --graph
build/graph/pfull_graph_edges.parquet`; seeds 4707-4711) and re-scores the benchmarks into `build/`. It is not
needed to reproduce any result: figures, tables and numbers come from the frozen outputs in `data/zenodo`. The
graph is reproduced exactly; the retrain is judged by the pre-registered parity test against the released model
(`tests/test_cold_tolerance.py`, `scripts/pfull_parity.py`; the authors' 5-seed retrain: test AUROC 0.8589 vs
0.8605, ΔAUROC -0.002 [-0.008, +0.005], PARITY). Smoke test:
`python -m disatlas.retrain --graph build/graph/pfull_graph_edges.parquet --out runs/smoke --seeds 4707 --epochs 1 --num-walks 2` (needs `requirements-retrain.txt`).

## Read by

<!-- READ_BY:BEGIN (generated by scripts/maintainer/data_index.py; do not edit) -->

| folder | files | read by |
|---|---|---|
| `derived/ablation/` | 10 | `figures/SI/si_ablation_31arm.py`, `figures/SI/si_mechanism_sensitivity.py`, `notebooks/03_si_robustness.ipynb` |
| `derived/control_pairs/` | 2 | `figures/main/fig5_disease_clustering_map.py`, `scripts/maintainer/stage_fig5.py` |
| `derived/degree_controls/` | 1 | `src/disatlas/rebuild/degree_controls.py`, `tables/SI/table_degree_controls.py` |
| `derived/degree_distributions/` | 5 | `figures/SI/si_mechdist_degree_by_channel.py`, `src/disatlas/rebuild/degree_distributions.py`, `src/disatlas/rebuild/powerlaw_csn.py` |
| `derived/fig1_open_targets/` | 5 | `figures/main/fig1_knowledge_transfer.py` |
| `derived/fig3_ixidn/` | 12 | `figures/SI/si_ixidn_matched_funnel.py`, `figures/main/fig3_ixidn_benchmark.py`, `notebooks/02_benchmarks_ixidn.ipynb`, `scripts/maintainer/stage_benchmarks.py`, `scripts/maintainer/stage_matched_qc.py`, `src/disatlas/benchmarks/ixidn_matched.py`, `src/disatlas/benchmarks/matched_qc.py`, `src/disatlas/rescore.py` |
| `derived/fig3a_trials/` | 1 | `tables/SI/table_fig3a_trials.py` |
| `derived/fig4_evidence_trace/` | 8 | `figures/SI/si_cluster_rarecommon_numbers.py`, `figures/main/fig4_evidence_trace.py`, `figures/main/fig4_numbers.py`, `tables/SI/table_SI_trace_order.py` |
| `derived/fig5_clustering/` | 3 | `figures/main/fig5_disease_clustering_map.py`, `figures/main/fig5c_export_txgnn_scores.py`, `notebooks/02_benchmarks_ixidn.ipynb`, `scripts/maintainer/stage_fig5.py` |
| `derived/fig6_comparison/` | 10 | `figures/main/fig3_ixidn_benchmark.py`, `figures/main/fig5_disease_clustering_map.py`, `figures/main/fig5c_export_txgnn_scores.py`, `figures/main/fig6_method_comparison.py`, `notebooks/02_benchmarks_ixidn.ipynb`, `scripts/maintainer/txgnn_partner_percentile.py`, `src/disatlas/rebuild/ccplus_strata.py`, `src/disatlas/rebuild/degree_controls.py`, `src/disatlas/rebuild/disease_bootstrap.py`, `src/disatlas/rebuild/fig6_inputs.py`, `src/disatlas/rebuild/mechdist_pair_score.py`, `src/disatlas/rebuild/placement_exclusion.py`, `src/disatlas/rebuild/si_baselines.py`, `src/disatlas/rebuild/si_heur_vs_ai.py`, `src/disatlas/rebuild/uncertainty.py`, `tables/main/table1_protocols.py` |
| `derived/fig7_prospective/` | 23 | `figures/SI/si_event_time_auc_methods.py`, `figures/SI/si_event_time_numbers.py`, `figures/SI/si_event_time_or_tables.py`, `figures/main/fig7_prospective.py`, `notebooks/04_si_biology.ipynb`, `scripts/maintainer/fig7_clogit_fit.R`, `scripts/maintainer/stage_fig7_or.py`, `src/disatlas/rebuild/event_time.py` |
| `derived/input_sources/` | 1 | `tables/SI/table_input_sources.py` |
| `derived/ixidn_funnel/` | 2 | `figures/SI/si_ixidn_matched_funnel.py`, `figures/main/fig3_ixidn_benchmark.py`, `notebooks/02_benchmarks_ixidn.ipynb` |
| `derived/ixidn_robustness/` | 4 | `figures/SI/si_node_disjoint.py`, `notebooks/03_si_robustness.ipynb`, `scripts/maintainer/stage_benchmarks.py`, `scripts/maintainer/stage_figD.py`, `src/disatlas/ontology_baseline.py`, `src/disatlas/rebuild/disease_bootstrap.py`, `src/disatlas/rebuild/si_heur_vs_ai.py`, `tables/SI/table_llm_subgroups.py` |
| `derived/mechanism_channels/` | 8 | `figures/SI/si_mechanism_sensitivity.py`, `notebooks/03_si_robustness.ipynb` |
| `derived/ontology_baseline/` | 5 | `scripts/maintainer/stage_ontology_baseline.py`, `src/disatlas/ontology_baseline.py`, `tables/SI/table_ontology_baseline.py` |
| `derived/pipeline_counts/` | 2 | `figures/SI/si_pipeline_detailed.py` |
| `derived/placement_confidence/` | 4 | `figures/SI/si_placement_confidence.py`, `scripts/maintainer/stage_figD.py`, `scripts/maintainer/stage_placement_confidence.py`, `src/disatlas/rebuild/placement_exclusion.py` |
| `derived/popularity_strata/` | 4 | `figures/SI/si_popularity_strata.py`, `scripts/maintainer/stage_figD.py` |
| `derived/si_baselines/` | 2 | `scripts/maintainer/stage_baseline_pair_scores.py`, `scripts/maintainer/stage_figD.py`, `src/disatlas/rebuild/si_baselines.py`, `src/disatlas/rebuild/si_heur_vs_ai.py`, `tables/SI/table_si_baselines.py` |
| `derived/si_cluster/` | 7 | `figures/SI/si_cluster_map_vs_model.py`, `figures/SI/si_cluster_rarecommon_numbers.py`, `figures/SI/si_cluster_validation.py`, `notebooks/04_si_biology.ipynb`, `scripts/maintainer/stage_si_cluster.py`, `scripts/slurm/k_selection_boot_array.sh`, `src/disatlas/rebuild/cluster_mondo_concordance.py`, `src/disatlas/rebuild/cluster_within_between.py`, `src/disatlas/rebuild/k_selection_boot.py`, `tables/SI/table_cluster_controls.py` |
| `derived/si_cluster_themes/` | 6 | `scripts/maintainer/stage_si_cluster.py`, `src/disatlas/rebuild/cluster_themes.py`, `tables/SI/table_cluster_themes.py` |
| `derived/si_event_time/` | 8 | `figures/SI/si_cluster_rarecommon_numbers.py`, `figures/SI/si_event_time_numbers.py`, `figures/SI/si_event_time_or_tables.py`, `figures/SI/si_event_time_tables.py`, `figures/SI/si_rare_common_neighbour_coverage.py`, `figures/main/fig7_prospective.py`, `notebooks/04_si_biology.ipynb`, `scripts/maintainer/stage_event_time.py`, `scripts/maintainer/stage_fig7_or.py`, `scripts/maintainer/stage_si_rare_common.py`, `scripts/pfull_parity.py`, `src/disatlas/rebuild/event_time.py`, `src/disatlas/rebuild/rc_coverage.py` |
| `derived/si_explainability/` | 2 | `figures/SI/si_explainability.py`, `src/disatlas/rebuild/case_scores.py` |
| `derived/si_fda/` | 7 | `figures/SI/si_bio_gap_numbers.py`, `figures/SI/si_fda_channels.py`, `figures/SI/si_fda_prioritisation_text.py`, `figures/SI/si_fda_schematic.py`, `figures/main/fig5_disease_clustering_map.py`, `notebooks/04_si_biology.ipynb`, `scripts/maintainer/stage_figD.py` |
| `derived/si_go_release/` | 3 | `figures/SI/si_go_release_table.py` |
| `derived/si_graph_channels/` | 2 | `figures/SI/si_graph_channels_schematic.py` |
| `derived/si_graph_channels/provenance_f/` | 3 | `figures/SI/si_graph_channels_schematic.py` |
| `derived/si_heur_vs_ai/` | 1 | `scripts/maintainer/stage_figD.py`, `src/disatlas/rebuild/si_heur_vs_ai.py`, `tables/SI/table_heuristic_vs_ai.py` |
| `derived/si_limitations/` | 2 | `figures/SI/si_bio_gap_numbers.py` |
| `derived/si_mechanism_timelines/` | 7 | `figures/SI/si_bio_gap_numbers.py`, `figures/SI/si_complement_milestones.py`, `figures/SI/si_mechanism_coverage_slices.py` |
| `derived/si_mechdist_pair_score/` | 4 | `figures/SI/si_mechdist_embedding_pair_score.py`, `src/disatlas/rebuild/mechdist_pair_score.py` |
| `derived/si_rare_common/` | 17 | `figures/SI/si_cluster_rarecommon_numbers.py`, `figures/SI/si_rare_common_neighbour_coverage.py`, `figures/SI/si_rare_common_transfer.py`, `notebooks/04_si_biology.ipynb`, `scripts/maintainer/stage_ccplus_pairs.py`, `scripts/maintainer/stage_partner_universe.py`, `scripts/maintainer/stage_rc_annotation_counts.py`, `scripts/maintainer/stage_si_rare_common.py`, `scripts/maintainer/txgnn_partner_percentile.py`, `src/disatlas/rebuild/ccplus_strata.py`, `src/disatlas/rebuild/rc_coverage.py`, `tables/SI/table_rare_common_strata.py` |
| `derived/threshold_sensitivity/` | 4 | `figures/SI/si_threshold_sensitivity.py`, `notebooks/03_si_robustness.ipynb` |
| `derived/training_graph/` | 1 | `figures/SI/si_mechdist_degree_by_channel.py`, `figures/SI/si_mechdist_embedding_pair_score.py`, `scripts/maintainer/stage_model.py`, `src/disatlas/rebuild/degree_controls.py`, `src/disatlas/rebuild/mechdist_pair_score.py`, `src/disatlas/rebuild/si_baselines.py` |
| `derived/uncertainty/` | 2 | `figures/SI/si_uncertainty.py`, `notebooks/03_si_robustness.ipynb`, `scripts/maintainer/stage_figD.py`, `src/disatlas/rebuild/uncertainty.py` |
| `raw/benchmarks/comparators/` | 2 | `scripts/maintainer/stage_benchmarks.py`, `src/disatlas/rebuild/fig6_inputs.py` |
| `raw/benchmarks/` | 1 | `scripts/maintainer/stage_benchmarks.py`, `scripts/maintainer/stage_matched_qc.py`, `scripts/maintainer/stage_partner_universe.py`, `src/disatlas/benchmarks/ixidn_matched.py`, `src/disatlas/benchmarks/matched_qc.py`, `src/disatlas/graph/pfull_build.py`, `src/disatlas/rebuild/fig6_inputs.py`, `src/disatlas/rescore.py`, `src/disatlas/retrain.py` |
| `raw/benchmarks/fda_repurposing/` | 2 | `scripts/maintainer/stage_benchmarks.py`, `src/disatlas/graph/pfull_build.py`, `src/disatlas/rescore.py` |
| `raw/benchmarks/ixidn/` | 2 | `src/disatlas/rescore.py`, `src/disatlas/retrain.py` |
| `raw/benchmarks/ixidn_matched/` | 2 | `figures/main/fig3_ixidn_benchmark.py`, `notebooks/02_benchmarks_ixidn.ipynb`, `run_pipeline.sh`, `scripts/maintainer/stage_benchmarks.py`, `scripts/maintainer/stage_matched_qc.py`, `src/disatlas/benchmarks/ixidn_matched.py`, `src/disatlas/benchmarks/matched_qc.py`, `src/disatlas/rescore.py` |
| `raw/benchmarks/ixidn_matched_qc/` | 1 | `scripts/maintainer/stage_matched_qc.py`, `src/disatlas/benchmarks/ixidn_matched.py` |
| `raw/bikg_habcd/` | 2 | `notebooks/01_data_and_channels.ipynb`, `src/disatlas/graph/raw_build.py` |
| `raw/broad_channels/` | 2 | `src/disatlas/graph/raw_build.py` |
| `raw/channels/specific/` | 2 | `figures/SI/si_mechdist_degree_by_channel.py`, `notebooks/01_data_and_channels.ipynb`, `src/disatlas/graph/pfull_build.py`, `src/disatlas/rebuild/degree_distributions.py`, `src/disatlas/rebuild/powerlaw_csn.py` |
| `raw/gencc/` | 1 | `src/disatlas/graph/raw_build.py` |
| `raw/go/` | 1 | `figures/SI/si_rare_common_transfer.py`, `scripts/maintainer/stage_cluster_theme_inputs.py`, `src/disatlas/channels/broad.py`, `src/disatlas/rebuild/cluster_themes.py` |
| `raw/hpo/` | 2 | `figures/SI/si_rare_common_transfer.py`, `src/disatlas/channels/broad.py`, `src/disatlas/graph/raw_build.py`, `src/disatlas/rebuild/cluster_themes.py` |
| `raw/mondo/` | 4 | `scripts/maintainer/stage_cluster_theme_inputs.py`, `scripts/maintainer/stage_mondo_body_system.py`, `scripts/maintainer/stage_si_cluster.py`, `src/disatlas/rebuild/cluster_mondo_concordance.py`, `src/disatlas/rebuild/cluster_themes.py`, `tables/SI/table_cluster_themes.py` |
| `raw/opentargets/` | 2 | `src/disatlas/graph/raw_build.py` |
| `raw/orphanet/` | 1 | `scripts/maintainer/stage_fig7_or.py`, `src/disatlas/fig1_open_targets.py`, `src/disatlas/graph/raw_build.py` |
| `raw/public_disease_gene/` | 1 | `src/disatlas/graph/raw_build.py` |
| `raw/specific_channels/` | 2 | `scripts/maintainer/stage_cluster_theme_inputs.py`, `src/disatlas/graph/raw_build.py`, `src/disatlas/rebuild/cluster_themes.py` |
| `raw/training_graph/` | 1 | `figures/SI/si_mechdist_degree_by_channel.py`, `figures/SI/si_mechdist_embedding_pair_score.py`, `scripts/maintainer/stage_model.py`, `src/disatlas/rebuild/degree_controls.py`, `src/disatlas/rebuild/mechdist_pair_score.py`, `src/disatlas/rebuild/si_baselines.py` |

Per-file list: [DATA_INDEX.tsv](DATA_INDEX.tsv).
<!-- READ_BY:END -->
