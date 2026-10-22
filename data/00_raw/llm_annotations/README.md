# Frozen LLM annotations

Disease IDs are MONDO or ORPHANET IDs as in the source data, or `DIS_B…` codes for diseases with neither. These match the IDs in every other table. Prompts, models and the parser are in `code/02_llm_annotations/`.

| File | Rows | Key | Columns |
|---|---|---|---|
| `ixidn_judge_verdicts.parquet` | 4,302 | disease_a, disease_b | verdict (STRONG/PLAUSIBLE/WEAK/SUSPECT/NOT_IN_GRAPH), verdict_source, prior_verdict, keep |
| `disease_labels.parquet` | 13,021 | disease_id | is_rare, disease_class, is_cancer, is_malignant, rare_source, cancer_source (rule/curated/agent), asked_rare, asked_cancer |
| `disease_label_responses.parquet` | 1,130 | disease_id | response_is_rare, response_is_cancer, response_is_malignant, response_confidence |

md5: see `md5sums.txt`.
