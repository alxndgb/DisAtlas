# genes

- `gene_go.parquet` — gene-GO annotations used for GO distances; columns gene_id (`ENSEMBL:ENSG…`), go_id, namespace; 59,312 rows.
- `gene_pfam.parquet` — gene-Pfam domain memberships; columns gene_id (`ENSEMBL:ENSG…`), pfam_id; 33,621 rows.
- `gene_hpa.parquet` — gene cell-type expression profiles (centred, unit-norm log1p nCPM); columns gene_id (`ENSEMBL:ENSG…`) + 154 cell-type values; 17,294 rows.
- `gene_id_map.parquet` — gene identifier map joining the gene tables; columns gene_id (Ensembl ENSG), ncbi_gene_id, symbol; 19,558 rows.
