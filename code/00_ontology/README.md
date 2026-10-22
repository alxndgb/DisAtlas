# code/00_ontology

Ontology term distances used by MechDist (Methods, mechanism distance) and the ontology baselines
(Appendix, IxIDN-Matched ontology baseline). Reads only the pinned OBO files in
`data/00_raw/ontologies/` and `data/zenodo/00_raw/ontologies/`.

## API (`ontology.py`)

```python
import sys; sys.path.insert(0, "code/00_ontology")
import ontology as O

nodes, anc, remap = O.load_ancestors("go_bp_semi")   # shipped closure; also go_mf, go_cc, hpo
M = O.distance_matrix(nodes, anc)                    # int16 n x n, row/col order = nodes
d_max = int(M.max())                                 # normaliser (go_bp_semi 23, go_mf 18, go_cc 18, hpo 25)
canon = remap.get(term, term)                        # obsolete / alt id -> active id
```

- `build_go()` -> `{key: (nodes, ancestors, remap_df)}` and `build_hpo()` -> `(nodes, ancestors, remap_df)`
  rebuild the `ancestors/` tables from the OBO files (GO: is_a / part_of weight 1, regulates family
  weight 2, per namespace; HPO: is_a). Output equals the shipped tables row for row after sorting.
- `distance_matrix(nodes, anc)`: min-ancestral-path, `dist(x, y) = min_c d(x, c) + d(y, c)`.
- `lca_pair(ax, ay)`: distance and all tied lowest common ancestors for one pair.
- `read_obo`, `term_remap`, `parent_graph`, `ancestor_closure`, `depth_from_root`: building blocks.

## Disease-level distances (`disease_distance.py`)

```python
import disease_distance as DD
pairs = DD.pairs_sharing_entity(disease_gene, "gene_id")        # sorted (a, b), a < b
go  = DD.go_disease_distances(pairs, disease_gene, gene_go)     # dist_BP, dist_MF, dist_CC, dist_GO
hpo = DD.hpo_disease_distances(DD.pairs_sharing_entity(disease_symptom, "hpo_id"), disease_symptom)
```

- GO: gene-gene distance = min term LCA distance / d_max; disease-disease = symmetric best-match
  average (BMA) over genes; namespaces summed with weights BP 1.0, MF 0.25, CC 0.1. Terms annotated to
  more than 2% of disease genes are masked at distance time (a gene annotated only to masked terms
  keeps them). GO terms are remapped (alt_id / replaced_by / consider) before lookup.
- HPO: BMA of symptom LCA distances / d_max, no gene layer. As in the production network, symptom
  IDs are looked up on the matrix axis directly; obsolete or alternate HPO IDs are dropped.

## HPO semantic similarity (`semsim.py`)

`hpo_scores(ids, disease_symptom, anc, remap, M, idx, d_max)` -> `{hpo_resnik_bma, hpo_lin_bma,
hpo_simgic, mechdist_hpo}` (n x n float32). Corpus IC over all annotated diseases; Resnik / Lin with
symmetric BMA; SimGIC = IC-weighted Jaccard of closures; `mechdist_hpo` = negated HPO BMA distance.

## Scripts

| Script | Output (`data/00_rebuilt_ontology/`) | Manuscript section |
|---|---|---|
| `build_ancestors.py` | `ancestors/{key}_ancestors.parquet`, `_nodes.csv`, `_term_remap.csv`, `ancestor_build_summary.csv` | Methods, ontological distance; Figure SI_lca |
| `build_disease_distances.py` | `disease_go_distance.parquet`, `disease_hpo_distance.parquet` | Methods, MechDist edges E1 and E4 |
| `build_hpo_baseline.py` | `ixidn_matched_hpo_scores.csv`, `ixidn_matched_hpo_baseline.csv` | Appendix SI_ixidn_matched_ontology_baseline (Table SI_ontology_baseline, HPO rows) |

All three are deterministic and call no external service. Parity tests: `tests/test_ontology_parity.py`.
