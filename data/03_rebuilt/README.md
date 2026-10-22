# MechDist rebuild outputs (Methods, mechanism distance; Results, MechDist network)

Written by `code/03_network/` from `data/00_raw/`.

| File | Content |
|---|---|
| `mechdist_network_edges.parquet` | disease_a, disease_b, channel_class (which of GO/HPO/px the pair has), crank (fused rank), origin (`topN`, `floor`), weight |
| `mechdist_network_nodes.parquet` | disease_id, data_support |
| `S3_thresholds.json` | distance cut-offs and per-channel distribution summary |
| `S5_px_budget.json` | per-channel-class edge quotas used in S6 |
| `S7_repair_report.json` | connectivity repair log |
| `work/` | S1-S5 intermediate tables (GO, HPO, candidate pairs, CRank frame); regenerated, not tracked |
| `frozen/disease_sort_rank.parquet` | stored input: each disease's rank in the original identifier order, which fixes the order of tied pairs (see `code/03_network/README.md`) |
| `reference/` | production network edge/node lists and a SHA-256 of the production S5 frame, used by `tests/test_mechdist_parity.py` |
