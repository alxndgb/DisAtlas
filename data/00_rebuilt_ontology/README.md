# data/00_rebuilt_ontology

Outputs of `code/00_ontology/` rebuilt from `data/00_raw/` and `data/zenodo/00_raw/ontologies/`.
Disease IDs are the strings used in `data/00_raw/associations/`; pair IDs (`blind_id`) are those of
`data/08S_eval/ixidn4364_scores.tsv`.

| File | Content | Manuscript section |
|---|---|---|
| `ancestors/{go_bp_semi,go_mf,go_cc,hpo}_ancestors.parquet` | `node, ancestor, dist` ancestor-or-self closure (GO weighted by relation; HPO is_a) | Methods, ontological distance; Figure SI_lca |
| `ancestors/*_nodes.csv`, `*_term_remap.csv` | matrix axis; obsolete / alternate ID -> active ID | same |
| `ancestors/ancestor_build_summary.csv` | terms, closure rows, d_max (23 / 18 / 18 / 25) | same |
| `disease_go_distance.parquet` | 114,961 pairs sharing a gene: `dist_BP, dist_MF, dist_CC, dist_GO` | Methods, MechDist edge E1 |
| `disease_hpo_distance.parquet` | 13,207,984 pairs sharing a symptom: `dist_HPO` (about 120 MB; Zenodo-sized) | Methods, MechDist edge E4 |
| `ixidn_matched_hpo_scores.csv` | per-pair HPO Resnik-BMA, Lin-BMA, SimGIC, direct MechDist HPO similarity; empty = unscorable | Appendix SI_ixidn_matched_ontology_baseline |
| `ixidn_matched_hpo_baseline.csv` | AUROC with 95% disease-clustered bootstrap CIs, scorable subset and all pairs, all / test-test scopes | Table SI_ontology_baseline (HPO rows) |
| `ixidn_matched_ordo_mondo_distance.tsv` | per-pair ORDO and MONDO max-LCA distance (longest-path legs); -1 = not in the ontology | Appendix, Dependence between ORDO and HPO; Negative-set realism |
| `ontology_dependency_baselines.csv` | ontology-only classifiers (HPO, ORDO, MONDO): AUROC on scorable pairs and on all pairs, with DisAtlas on the same pairs and 95% disease-clustered bootstrap CIs | Table SI_ontology_baseline |
| `ontology_dependency_correlations.csv` | Spearman correlation between ORDO / MONDO distance, HPO similarity and the DisAtlas score | Dependence between ORDO and HPO (Spearman 0.51-0.62) |
| `ontology_dependency_discordant.csv` | least-concordant ORDO-vs-HPO subsets; pairs without HPO annotation | Dependence between ORDO and HPO; Table SI_ontology_baseline |
| `ontology_dependency_orphanet.csv` | pairs inside / outside the Orphanet classification | Outside-Orphanet control (0.786 / 0.810) |
| `ontology_dependency_annotation_source.csv` | share of MechDist HPO annotations on Orphanet-keyed diseases (99.1%) | Dependence between ORDO and HPO |
| `negative_realism_ordo_distance.tsv.gz` | ORDO distance of random ORDO disease pairs and IxIDN-Matched positives / negatives | Negative-set realism, Figure SI_negrealism (a) |
| `negative_realism_hard_negative_auroc.tsv` | DisAtlas AUROC against progressively harder negative subsets, and within-stratum baselines | Negative-set realism, Figure SI_negrealism (b) |
| `negative_realism_summary.tsv` | medians and shared organ-system / ORDO-area shares quoted in the text | Negative-set realism |

All values equal the production tables exactly (see `tests/test_ontology_parity.py`).
