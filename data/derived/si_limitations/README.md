# Limitations: single-gene pairs and the rare-shared-gene link ablation

Small precomputed tables that cannot be regenerated inside this repository (upstream model runs or benchmark construction). Each line gives the content and the scripts that read the file.

- `gene_anchor_ablation_summary.tsv`: five-seed retrain adding a specific-layer link for pairs sharing a rare high-confidence gene: top-100 share, target pairs, FDA development pairs lost, IxIDN-Matched AUROC change (paired bootstrap). Read by: si_bio_gap_numbers.py.
- `single_gene_pairs_census.tsv.gz`: census of 14,805 classic single-gene disease pairs with rank and flag columns. Read by: si_bio_gap_numbers.py.
