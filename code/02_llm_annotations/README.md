# 02 LLM annotations

The model outputs used in the paper are shipped as frozen tables in `data/00_raw/llm_annotations/`. No rebuild script calls an API.

| Use (manuscript section) | Model | Frozen table |
|---|---|---|
| Keep/drop check of the 4,302 candidate IxIDN-Matched positives (Methods, IxIDN-Matched construction) | Claude Sonnet 5 | `ixidn_judge_verdicts.parquet` |
| Rare/common, neoplasm and malignant labels for the 13,021 diseases (Methods, disease labels) | Claude Opus 5.5 | `disease_labels.parquet`, `disease_label_responses.parquet` |

- `prompts/`: system prompts, user-message templates and the label tool schema, exactly as sent.
- `models.json`: model IDs and request settings.
- `parse_responses.py`: rebuilds the shipped labels from the saved responses and checks they match. Run `python code/02_llm_annotations/parse_responses.py`.
- `call_llm.py`: optional re-run. The output is not deterministic and is not used by the rebuild.

## Labels

A disease's rare/common label comes from the curated sheet when there is one. Otherwise it comes from the rule. Diseases with missing or conflicting rule signals (`asked_rare`, `asked_cancer`) were sent to the model. A `rare`/`common` (or `yes`/`no`) answer replaced the rule value. An `uncertain` answer kept the rule value. `rare_source` and `cancer_source` record which of these set each label.

## Judge

A pair is kept when its verdict is `STRONG` or `PLAUSIBLE`. `NOT_IN_GRAPH` pairs were not sent to the model. `prior_verdict` holds the verdict from the earlier pass for the 944 re-judged pairs.
