# Approved-repurposing (FDA) benchmark

Small precomputed tables that cannot be regenerated inside this repository (upstream model runs or benchmark construction). Each line gives the content and the scripts that read the file.

- `benchmark_pairs_split.tsv`: benchmark disease pairs with control class and split (counts in the FDA schematic). Read by: si_bio_gap_numbers.py, si_fda_schematic.py.
- `broad_only_fda_dev_pair_scores_ens5.tsv`: per-pair percentiles for the full, specific-off and broad-off arms (five-seed ensembles). Read by: si_fda_channels.py.
- `broad_only_fda_dev_results.tsv`: recovery at the 25/5/1% cut-offs per arm with bootstrap intervals. Read by: si_fda_channels.py.
- `channel_arm_universe.tsv`: number of diseases in each channel-ablation arm's node set. Read by: si_bio_gap_numbers.py.
- `endpoints.tsv`: prioritisation endpoints per method. Read by: si_fda_channels.py, si_fda_prioritisation_text.py.
- `flow.tsv`: benchmark construction flow counts. Read by: si_fda_prioritisation_text.py.
- `paired_differences.tsv`: paired method differences with bootstrap intervals. Read by: si_fda_prioritisation_text.py.

Licence: a pair is shipped with real MONDO ids and labels when its shared drug is on the US label for both diseases (openFDA) or both diseases are Open Targets/ChEMBL launched indications of a common drug. Diseases that occur in any pair resting on licensed launch evidence are replaced by stable blinded tokens (`T_<12hex>`, src/disatlas/blinding.py, the same token in both tables) and their labels are left empty; counts, scores, classes and splits are unchanged.
