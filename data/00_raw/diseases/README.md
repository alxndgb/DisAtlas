# diseases

- `diseases.parquet` — one row per disease ID string used in any table under `data/00_raw/` (13,680 rows); columns disease_id, canonical_id (the `ORPHA:` row for an `ORPHANET:` spelling of the same Orphanet disease, otherwise disease_id), label (null for blinded `DIS_B…` IDs), is_disease_indication (trial-network diseases only), in_dd_network, in_mechdist_network, in_mechdist_universe.
- `orpha_mondo_crosswalk.parquet` — Orphanet-to-MONDO crosswalk; columns orpha_id, mondo_id, rank (0 = preferred); 8,496 rows (8,484 Orphanet codes).
