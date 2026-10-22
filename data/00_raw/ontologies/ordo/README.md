# ORDO (Orphanet classification) — minimal inputs for the max-LCA distance

Post-filter extract of the Orphanet Rare Disease Ontology needed to recompute the ORDO distances reported
in the appendix ("Dependence between ORDO and HPO", "Ontology-only classifier and concordance checks",
"Negative-set realism"). Orphanet codes are public integer IDs (`ORPHA:<code>`). Used by `code/00_ontology/max_lca.py`.

| File | Content |
|---|---|
| `ordo_nodes.csv` | `orpha_id`: the 7,836 active disorder-level ORDO classes (distance axis) |
| `ordo_ancestors.parquet` | `node, ancestor, dist`: longest-path ancestor-or-self closure of each disorder in the classification graph (`part_of` + `subClassOf` edges whose parent is not a typing metaclass; pseudo-root `-1` one hop above the 30 top categories); 129,631 rows |
| `mondo_lca_axis.csv` | `mondo_id`: MONDO diseases that receive a MONDO max-LCA distance |
| `random_pair_diseases.csv` | `orpha_code`: ordered disease list from which the random reference pairs of Figure SI_negrealism (a) are drawn |
| `mechdist_hpo_annotation_counts.csv` | `disease_id, n_hpo_annotations`: HPO annotation count per disease of the MechDist graph (4,710 diseases, 125,998 annotations); non-public disease IDs are opaque (`DIS_B…`) |

Distance: `dist(x, y) = min over common ancestors c of longest(c→x) + longest(c→y)`.
MONDO diseases map to ORDO via `data/00_raw/diseases/orpha_mondo_crosswalk.parquet` (first matching row).
