# associations

- `disease_gene.parquet` — disease-gene associations used by MechDist; columns disease_id, gene_id (Ensembl gene ID, `ENSEMBL:ENSG…`), evidence_tier, tier_weight, cap_restored, num_sources, source_codes (pipe-joined opaque source codes); 25,972 rows.
- `disease_symptom.parquet` — disease-phenotype associations used by MechDist; columns disease_id, hpo_id, evidence_tier, tier_weight, cap_restored, num_sources; 244,394 rows.
- `disease_gene_tier1_sources.parquet` — tier-1 disease-gene pairs, one row per supporting source; columns disease_id, gene_id, source_code; 32,494 rows (18,287 pairs).
- `disease_symptom_tier1_sources.parquet` — tier-1 disease-phenotype pairs, one row per supporting source; columns disease_id, hpo_id, source_code; 271,340 rows (167,798 pairs).
