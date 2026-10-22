# Data licences and attribution

The tables in `data/00_raw/` are derived from the public resources below. Identifiers that come from non-public sources are replaced by opaque codes (`DIS_B…`, `DRUG_B…`, `TRIAL_B…`, `SRC_…`) and carry no third-party content.

| Resource | Used in | Licence |
|---|---|---|
| Mondo Disease Ontology | `ontologies/`, `diseases/`, MONDO disease IDs and labels | CC BY 4.0 |
| Human Phenotype Ontology | `ontologies/hp.obo`, `associations/disease_symptom*` | HPO licence (free use with attribution) |
| Gene Ontology | `ontologies/go-basic.obo`, `genes/gene_go.parquet` | CC BY 4.0 |
| Orphanet | Orphanet disease IDs, `diseases/orpha_mondo_crosswalk.parquet` | CC BY 4.0 |
| ChEMBL | ChEMBL drug IDs, names and InChIKeys in `drugs/drugs.parquet` | CC BY-SA 3.0 |
| NCBI Gene | `genes/gene_id_map.parquet` | Public domain (NCBI) |
| Ensembl | Ensembl gene IDs | No restrictions |
| Pfam / InterPro | `genes/gene_pfam.parquet` | CC0 |
| Human Protein Atlas | `genes/gene_hpa.parquet` | CC BY-SA 3.0 |
| ClinicalTrials.gov | NCT trial IDs in `trials/` | Public domain (U.S. NLM) |
| EU Clinical Trials Register / CTIS | EudraCT and CTIS trial IDs in `trials/` | EMA reuse terms (attribution) |

Derived tables (distances, labels, benchmark pairs, model outputs) are released under the repository licence.
