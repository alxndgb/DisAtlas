# Ontologies

Pinned ontology files, byte-identical to the releases used. Checksums: `MD5SUMS.txt` (`zenodo:` = file is in `data/zenodo/00_raw/ontologies/`).

| File | Location | Content |
|---|---|---|
| `hp.obo` | git | Human Phenotype Ontology (OBO 1.2), 11.2 MB |
| `mondo.obo` | zenodo | MONDO (OBO 1.2), 52.0 MB; disease hierarchy for the IxIDN+ parent/child co-investigation check and the MONDO ontology baseline |
| `mondo_propagation.obo` | zenodo | MONDO (OBO 1.2, later release), 53.1 MB; used only for descendant-to-ancestor symptom propagation |
| `go-basic.obo` | zenodo | Gene Ontology go-basic (OBO 1.2), 32.2 MB; under 50 MB but routed to Zenodo as a borderline file |
| `ancestors/go_bp_semi_ancestors.parquet` | git | GO BP ancestor closure (is_a/part_of weight 1, regulates weight 2): `node, ancestor, dist`; 480,463 rows |
| `ancestors/go_mf_ancestors.parquet` | git | GO MF ancestor closure: `node, ancestor, dist`; 73,818 rows |
| `ancestors/go_cc_ancestors.parquet` | git | GO CC ancestor closure: `node, ancestor, dist`; 46,358 rows |
| `ancestors/hpo_ancestors.parquet` | git | HPO is_a ancestor closure: `node, ancestor, dist`; 222,576 rows |
| `ancestors/go_bp_semi_nodes.csv` | git | GO BP term list (matrix order): `go_id`; 24,129 rows |
| `ancestors/go_mf_nodes.csv` | git | GO MF term list: `go_id`; 10,041 rows |
| `ancestors/go_cc_nodes.csv` | git | GO CC term list: `go_id`; 4,075 rows |
| `ancestors/hpo_nodes.csv` | git | HPO term list: `hpo_id`; 19,836 rows |
| `ancestors/go_bp_semi_term_remap.csv` | git | obsolete/alt GO BP IDs to canonical: `source_id, canonical_id, reason`; 3,494 rows |
| `ancestors/go_mf_term_remap.csv` | git | as above, GO MF; 2,371 rows |
| `ancestors/go_cc_term_remap.csv` | git | as above, GO CC; 589 rows |
| `ancestors/hpo_term_remap.csv` | git | as above, HPO; 4,020 rows |
| `ordo/` | git | Minimal ORDO (Orphanet) extract for the appendix ORDO/HPO dependence, ontology-only baseline and negative-set realism analyses; see `ordo/README.md` |

The `ancestors/` tables are produced from `go-basic.obo` / `hp.obo` by `code/00_ontology/` and are shipped as a parity reference.
