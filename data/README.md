# DisAtlas data

The small tables live in this repository under `data/raw` and `data/derived`. The training graph, the
frozen model and the large score files are in the Zenodo record and unpack into `data/zenodo`, which is
not tracked by git. `data/zenodo/INDEX.tsv` lists every Zenodo file, and `data/ZENODO_SHA256SUMS.txt`
gives their checksums (`cd data/zenodo && sha256sum -c ../ZENODO_SHA256SUMS.txt`).

## Folder map

```
data/
  raw/
    channels/
      broad/        broad (shared gene / shared phenotype) channel, post-threshold
      specific/     specific (shared mechanism module / mechanism node) channel, post-threshold
    benchmarks/
      ixidn/            IxIDN+ and IxIDN- pair lists (removed from the training graph)
      ixidn_matched/    IxIDN-Matched evaluation pairs (4,364) and the language-model curation verdicts
      fda_repurposing/  FDA-approved repurposing positive controls (holdout split withheld)
      ordo/             ontology distances for IxIDN-Matched pairs
      comparators/      stored per-run outputs of the comparator systems on IxIDN-Matched
      disease_labels.tsv
  derived/
    event_time/     event-time (prospective) evaluation result tables, precomputed
  zenodo/           (download) graph/, model/, scores/
```

## Licence and blinding

Open Targets identifiers are released unchanged. Targets and modules that come from licence-restricted
sources (curated pathway databases, the literature-mining knowledge graph and similar) cannot be
redistributed. They are replaced by stable salted tokens: `T_<12 hex>` for targets and `M_<12 hex>` for
modules. The same source identifier always maps to the same token, so overlaps between diseases and
all graph computations are preserved exactly, but the token cannot be mapped back without the private
salt. For these rows the source is given only at group level (`source_group`, `families_group`). No
publication identifiers or evidence sentences are included. Disease identifiers are MONDO identifiers (EFO where no MONDO term exists).
`disatlas.blinding.is_blinded(x)` tells you whether an identifier is a token.

## Column dictionaries

### raw/channels/broad

| file | column | meaning |
|---|---|---|
| broad_edges.parquet | u, v | disease pair joined in the broad channel (u < v) |
| broad_nodes.parquet | disease_id | the 13,021 diseases of the network |
| disease_gene.parquet | disease_id, target | disease-gene annotation (target: Ensembl/symbol, or a `T_` token) |
| | evidence_tier, tier_weight | evidence tier (1 = strongest) and the weight used to build the channel |
| | source_group | public source(s) joined by `|`; `literature` / `tier_rescue` for restricted rows |
| disease_phenotype.parquet | disease_id, hpo_id, evidence_tier, tier_weight | disease-phenotype (HPO) annotation |

### raw/channels/specific

| file | column | meaning |
|---|---|---|
| disease_module_links.parquet | disease_id, module_id | disease linked to a mechanism module (REACTOME, named and mechanism-node modules are raw, restricted ones `M_` tokens) |
| | module_order | row order within the build, needed to reproduce floating-point sums exactly |
| | scope | module scope used by the projection |
| | n_fam, docs_max | number of evidence families and maximum supporting-document count; link weight e = n_fam + log1p(docs_max) |
| | families_group | evidence families at group level |
| disease_target_mechanism.parquet | disease_id, mechanism_node, target, target_type, source_group, evidence_count | disease-target evidence placed on mechanism nodes |
| mechanism_node_classes.tsv | node, mechanism_class | mechanism class of each mechanism node (used for class balancing) |
| gene_set_counts.parquet | disease_id, module_id, module_order | disease x gene set: REACTOME leaf pathways with 5-500 background genes and at least one disease gene (all tested pathways, enriched or not), and kept NAMED gene sets |
| | n_genes, docs_max, families | disease genes in the set, maximum supporting-document count over them, evidence families present (no gene identifiers) |
| reactome_background.parquet | kind, id, n_genes | enrichment background: `all` = background genes N, `pathway` = pathway size K, `disease` = disease genes n |
| direct_module_evidence.parquet | disease_id, module_id, module_order, family, docs | kept literature, pathology-text (`M_` tokens) and mechanism-node modules per evidence family; docs = dated-document count (no document ids, dates or text) |

`python -m disatlas.channels.specific --out DIR` rebuilds disease_module_links.parquet from the last three tables (hypergeometric REACTOME enrichment, BH q < 0.05 over the pathways with k >= 1; keep n_fam >= 2 or docs >= 2); tests/test_specific_channel.py asserts it equals the shipped file.

### raw/benchmarks

| file | columns |
|---|---|
| ixidn/ixidn_plus_pairs.parquet, ixidn_minus_pairs.parquet | disease_a, disease_b |
| ixidn_matched/ixidn_matched_pairs.tsv | pair_id, disease_a, disease_b, y (1 positive / 0 negative), split (tune / test), pair_class (rare/common composition) |
| ixidn_matched/llm_curation_verdicts.tsv | disease_a, disease_b, verdict, rejudged, first_pass_verdict, keep (kept positives = IxIDN-Matched positives) |
| fda_repurposing/fda_repurposing_pairs.tsv | disease_a, disease_b: all 2,000 control pairs, unlabelled, removed from the training graph (needed to rebuild it exactly) |
| fda_repurposing/fda_repurposing_controls.tsv | disease_a, disease_b, stratum (hard / easy), pair_class, split (dev / easy_unsplit / complement_six), fda_verified_pair, in_model; the 560 holdout pairs are sealed and not included |
| ordo/ixidn_matched_ontology_distance.tsv | pair_id, ordo_dist (-1 = undefined), mondo_dist |
| ordo/ordo_distance_summary.tsv | set, n, q1, median, q3, mean (random-pair reference distribution) |
| comparators/ixidn_matched_comparator_scores.tsv | pair_id, then one column per `system|run` |
| disease_labels.tsv | disease_id, is_rare, is_cancer, is_malignant, disease_class |

### derived/event_time

Result tables of the event-time evaluation (E1: percentile rank of the later indication among candidates;
E2: precision measures), per pool, scope (pooled or year), analysis and comparator, with bootstrap
confidence intervals. The trial and approval events behind them come from a licensed source and are not
redistributed, so these tables are provided precomputed.

### zenodo/graph

| file | columns |
|---|---|
| union_graph_edges.parquet | u, v, in_broad, in_spec, w_spec (specific-channel weight), f (class-balancing factor), w_arm (walk weight); 498,155 edges, 13,021 diseases |
| training_graph_edges.parquet | same columns, after removing IxIDN+/IxIDN- pairs; this is the graph the ensemble was trained on |
| mechanism_class_balancing.tsv, graph_build_stats.json | class factors and build counts |

Rebuild and check: `python -m disatlas.graph.build --data data --out <dir> --check`.

### zenodo/model/disatlas_v3.2_ensemble

`embeddings/node2vec_s<seed>.npz` (ids, X: unit-norm 32-d embeddings), `decoders/decoder_s<seed>.safetensors`
(Hadamard MLP decoder, config in the file metadata), `nodes.txt`, `model_config.json` (recipe and
SHA-256 of every file). Seeds 4707-4711. Load with `disatlas.model.load_ensemble()`; `score_pairs(A, B)`
returns the ensemble probability and `top_partners(d, k)` the highest-scoring partners.

### zenodo/scores

| file | content |
|---|---|
| disatlas_all_pairs_ensemble.npz | ids and p_upper: strict upper triangle of the 13,021 x 13,021 ensemble probability matrix (`disatlas.pair_matrix.square`) |
| ixidn_matched_disatlas_full_scores.tsv | pair_id, p_s4707..p_s4711, p_ensemble (benchmark diseases in training, benchmark pairs removed) |
| ixidn_matched_disatlas_nodedisjoint_scores.tsv | same, benchmark diseases held out of training (five folds) |

## Retraining (optional)

`python -m disatlas.retrain --out runs/retrain` re-runs the published recipe (node2vec dim 32, walk length
20, window 2, p = 1, q = 0.5, seeds 4707-4711, Hadamard MLP decoder) on the shipped training graph. It is
not needed to reproduce any result. Its random streams differ from the frozen run, so the retrained
scores match the shipped ones statistically, not exactly. Smoke test:
`python -m disatlas.retrain --out runs/smoke --seeds 4707 --max-nodes 600 --epochs 1 --num-walks 2`.
