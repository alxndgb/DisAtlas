# Stage 1-4: GO/HPO disease distances and stage records

Stored tables (already thresholded; nothing here is downloaded). Disease IDs from
non-public vocabularies are blinded (`DIS_B...`), matching `data/00_raw/`.

- `disease_go_distance.parquet.gzip` -- pairwise GO disease distance for disease pairs sharing at least one
  gene; parity reference for `code/00_ontology/disease_distance.py` (`tests/test_ontology_parity.py`).
  The HPO counterpart (>50 MB) is `../zenodo/network/disease_hpo_distance.parquet.gzip`.
- `S0_inputs_manifest.json`, `S3_thresholds.json`, `S5_px_budget.json`, `S7_repair_report.json` --
  per-stage threshold records from the original pipeline run (CUT_GO = CUT_HPO = 0.2), kept for audit.
  `code/03_network/` regenerates S3/S5/S7 under `data/03_rebuilt/`.

The disease-gene and disease-symptom association tables and the disease labels are in
`data/00_raw/associations/` and `data/00_raw/llm_annotations/`.
