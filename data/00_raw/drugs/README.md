# drugs

- `drugs.parquet` — every drug node of the trial network; columns drug_id, name, inchikey (pipe-joined; name/inchikey given for ChEMBL IDs only); 38,389 rows.
- `drug_moa_decision.parquet` — per-drug mechanism-of-action filter decision used by IxIDN; columns drug_id, category, non_moa, is_biotherapeutic, disease_degree, decision (keep / flag_review / exclude); 38,389 rows.
