# data/raw — inputs of the MechDist graph build (Habcd release)

The broad channel is rebuilt from the files below by `python -m disatlas.graph.raw_build --out build/rawgraph`
(test: `tests/test_raw_graph.py`); the training graph is then built from it and the released specific-channel
edges by `python -m disatlas.graph.pfull_build --out build/graph` (test: `tests/test_pfull_graph.py`). All tables are post-threshold: only rows that enter the build are shipped,
and only the columns the build reads. Row and column reductions are listed in `release/MINIMISATION_rawgraph.tsv`.
Disease identifiers are MONDO (with EFO/Orphanet/MeSH/OMIM where no MONDO term exists); genes are Ensembl
(`ENSEMBL:ENSG…`).

## Disease-gene edges (one folder per source)

| folder | file | rows | content |
|---|---|---|---|
| `gencc/` | `disease_gene.parquet` | 7,319 | tier-1 curated edges supported by GenCC |
| `orphanet/` | `disease_gene.parquet` | 8,053 | tier-1 curated edges present in public Orphanet en_product6 |
| `public_disease_gene/` | `disease_gene.parquet` | 308 | further tier-1 curated edges, each present in the public Open Targets 26.09 direct evidence below the thresholds of the `opentargets/` files; `source`: Open Targets |
| `opentargets/` | `disease_gene_direct.parquet` | 3,694 | Open Targets 26.09 direct associations: overall score without literature (Europe PMC) ≥ 0.6, or, when a single data source supports the pair, that source's score ≥ 0.7; at most 100 per disease; only pairs not already curated (tier 1) |
| `opentargets/` | `disease_gene_relaxed.parquet` | 4,091 | Open Targets 26.09 direct, tier-2/3 rescue ladder relaxed to score ≥ 0.3, k = 1 (setting R3); `evidence_tier` 2–4 |
| `bikg_habcd/` | `habcd_bikg_gene_edges.parquet` | 8,343 | the released knowledge-graph edges (`disease_id`, `gene_symbol`, `evidence_tier`): the approved release with its `layer` column removed, rows and values unchanged; checksum in `.sha256` |

The gencc and orphanet files overlap (an edge with both supports is in both). Assembly: union of the
tier-1 files, plus the relaxed rows (lowest tier kept per pair), then the released edges for pairs not yet covered. `tier_weight` is derived from `evidence_tier`
(1: 1.0, 2: 0.6, 3: 0.35, 4: 0.15).

The assembled disease-gene table (`raw_build.gene_edges()`, 29,409 rows) is 21,066 rows from the public sources
above plus the 8,343 released knowledge-graph rows, and nothing else. Of its 23,840 tier-1 rows, 16,975 are public
and 6,865 come from the released edges.

Of the 918 curated disease–gene pairs used for the graph of record, the 604 supported by OMIM alone are not redistributed for licence reasons, and 6 more with no public source are not released. The 308 released pairs (`data/raw/public_disease_gene/`) are all in the public Open Targets 26.09 direct evidence, below the thresholds used to build the shipped Open Targets files. The remaining released BIKG-derived data are shared with permission from the AstraZeneca data office. All reported results reproduce exactly from the shipped graph files and model outputs. Rebuilding the gene-anchored graph from the public inputs here changes 36,868 edges (6.9%) relative to the training graph of the released model; a 5-seed retrain on a rebuild that differs from this one by 314 edges (0.06%) is within retraining noise (IxIDN-Matched test ΔAUROC −0.0016, 95% CI −0.0074 to +0.0041). The note "rebuilt from shipped inputs" in the shipped `model_config.json` predates the exclusion of the OMIM-only pairs; the final graph reproduces from the shipped graph files.

## Disease-phenotype edges (`hpo/`)

`disease_phenotype.parquet` (202,778 rows): disease-phenotype edges supported by public HPO annotations
(phenotype.hpoa: aspect P, not negated, frequency ≥ 1 %, or an OMIM-derived HPOA annotation). `disease_phenotype_propagated.parquet` (39,471 rows, tier 9): annotations propagated to a
disease from a MONDO child disease whose own annotation is kept.

## Ontology term distances and gene GO annotations (`specific_channels/`; historical name) — inputs of the broad channel

`gene2go.parquet` (gene, go_id, namespace): GO annotations of genes with a disease edge, restricted to terms in
the term-distance slices. `term_distance_params.json`: d_max of the full term-distance matrices (the slices
cover only the terms in use, so d_max is fixed here), the obsolete/alt-id remap and the namespace weights.
The term-distance slices (int16,
GO BP semi-structural, MF, CC and HPO; 6,431 / 2,134 / 1,162 / 10,646 terms) are large and ship in the Zenodo
bundle under `data/zenodo/ontology_term_term_distance/` (checksums in `release/ZENODO_SHA256SUMS_rawgraph.txt`).

## Broad channels (`broad_channels/`) — gene-proximity distances

`gene_hpa_celltype_profile.parquet`: Human Protein Atlas single-cell-type profiles (centred, unit norm) of the
5,871 disease genes with a profile. `gene_pfam.parquet`: UniProt-reviewed Pfam domains of disease genes.

## Specific channel (`channels/specific/`) — released edges

`specific_edges.parquet` (189,903 rows; `u`, `v`, `w`): the disease-disease edges of the specific channel (shared
mechanism modules) and their specific-channel weight, as used by the P_full training graph; checksum in `.sha256`. Endpoints are MONDO identifiers except 49
native identifiers without a MONDO term, kept unmapped. Rebuilding this channel from its inputs would need about
1.06 million rows of licence-restricted knowledge-graph evidence (literature, pathology, MetaBase, CTD and
mechanism nodes), which is why only its disease-disease edges are released; a rebuild from public sources only
loses 0.066 test AUROC. The released edges on their own do not reproduce the manuscript model: the training graph
also needs the broad channel rebuilt above.

## Ontology extracts (`mondo/`, `go/`) — inputs of the cluster-theme and MONDO-concordance rebuilds

`mondo/`: MONDO term table, per-disease ancestors and body-system classes (MONDO release 2026-02-06), staged by
`scripts/maintainer/stage_mondo_body_system.py` and `scripts/maintainer/stage_cluster_theme_inputs.py`. `go/go_bp_terms.tsv.gz`:
Gene Ontology terms annotated in `specific_channels/gene2go.parquet` and their ancestors (go-basic.obo dated 2026-07-26).
Columns, counts and readers: `mondo/README.md`, `go/README.md`.

## Training graph of the v3.2 embedding (`training_graph/`)

`training_graph_edges.parquet` (488,103 rows, 13,021 diseases; `u`, `v`, `in_broad`, `in_spec`, `w_arm`): the graph the
v3.2 five-seed embedding was trained on, with the IxIDN-Matched benchmark pairs and the FDA pairs stripped. It differs
from the graph that `disatlas.graph.pfull_build` rebuilds from the inputs above (flagged for the authors);
`data/derived/training_graph/training_graph_stats.json` describes it before the strip (498,155 edges). Read by the
rebuilds that need training-graph degrees, neighbours or graph baselines (`degree_controls`, `mechdist_pair_score`,
`si_baselines`). Source, checksum and licence status: `training_graph/SOURCE.md`.
