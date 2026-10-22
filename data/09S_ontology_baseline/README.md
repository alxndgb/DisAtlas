# 09S — ontology-only baseline and ORDO/HPO concordance (stored outputs)

Stored production tables behind the appendix sections "Dependence between ORDO and HPO" and
"Ontology-only classifier and concordance checks" (Table SI_ontology_baseline). Pair IDs (`blind_id`)
are those of `data/08S_eval/ixidn4364_scores.tsv`. Each table is rebuilt byte-identically by
`code/00_ontology/` (outputs in `data/00_rebuilt_ontology/`; checked by `tests/test_ontology_parity.py`).

| File | Content | Manuscript number(s) | Rebuilt by |
|---|---|---|---|
| `ixidn4364_pairs_scored.csv` | per-pair HPO / MONDO SimGIC similarity | MONDO SimGIC label-only control; HPO rows of Table SI_ontology_baseline | `build_hpo_baseline.py` |
| `ixidn_matched_ordo_mondo_distance.tsv` | `blind_id, ordo_dist, mondo_dist`: ORDO / MONDO max-LCA distance (longest-path legs), -1 = not in the ontology | input to all checks below and to Negative-set realism | `max_lca.py` (via `build_ontology_dependency.py`) |
| `ontology_only_baselines.csv` | ontology-only classifiers (HPO, ORDO, MONDO): AUROC on scorable pairs and on all pairs (unscorable ranked last), DisAtlas on the same pairs, 95% disease-clustered bootstrap CIs | HPO 35% coverage, 0.69-0.72, all pairs 0.54, DisAtlas +0.03-0.06 (n.s.) | `build_ontology_dependency.py` |
| `ordo_hpo_correlations.csv` | Spearman correlation between ORDO / MONDO distance, HPO similarity and the DisAtlas score | Spearman 0.51-0.62 | same |
| `ordo_hpo_least_concordant.csv` | AUROC on the least-concordant half / third of ORDO-vs-HPO pairs; pairs without HPO annotation | DisAtlas 0.70-0.76, HPO 0.67-0.73, ORDO at chance | same |
| `outside_orphanet_control.csv` | DisAtlas AUROC inside / outside the Orphanet classification | 62% of pairs, 0.786 (0.810 test split) | same |
| `hpo_annotation_orphanet_share.csv` | share of MechDist HPO annotations on Orphanet-keyed diseases | 99.1% | same |

The ontologies themselves are public (MONDO `mondo.obo` in `data/zenodo/00_raw/ontologies/`; the ORDO
classification as the minimal ancestor table in `data/00_raw/ontologies/ordo/`).
