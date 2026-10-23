# Comparator benchmarks (Figure 6)

Comparisons of DisAtlas with language-model systems (`llm/`) and the published TxGNN checkpoint (`txgnn/`) on
IxIDN-Matched (4,364 disease pairs, 50% positive). Every number is rebuilt from the stored outputs used in the paper,
with no network access, model calls or API keys:

    python benchmarks/verify_benchmarks.py

This runs `llm/score_llm.py` and `txgnn/score_txgnn.py` and asserts each manuscript value at reported precision.
The answer key and full-DisAtlas scores come from `data/08S_eval/ixidn4364_scores.tsv`. `_common.py` holds the shared
metrics: AUROC; Figure 6b ranking metrics with a 500-draw pair bootstrap (seed 4707); Figure 6c run agreement;
Figure 6d top-10% overlap; disease-clustered Poisson bootstrap.

Re-querying the models (`llm/score_llm.py --rerun-llm`) and recomputing TxGNN scores (`txgnn/fetch_txgnn.sh`,
`txgnn/run_txgnn_inference.py`) are optional, and neither is needed to reproduce the paper.

| Manuscript item | Output |
|---|---|
| Results, cross-system comparison (pooled AUROC); Figure 6a | `llm/results/protocol2_auroc_by_subset.tsv` (all pairs), `txgnn/results/ranking_metrics.tsv` |
| Figure 6b ranking metrics with 95% CI | `llm/results/ranking_metrics.tsv`, `txgnn/results/ranking_metrics.tsv` |
| Figure 6c run-to-run agreement | `llm/results/run_agreement.tsv`, `run_agreement_summary.tsv` |
| Figure 6d top-10% overlap | `llm/results/top10_jaccard.tsv` |
| Supplementary "LLM subgroups" (test split, oncology, non-oncology) | `llm/results/protocol2_auroc_by_subset.tsv`, `protocol2_oncology_contrast.tsv` |
| Explainability contrast, TxGNN AUROC and coverage | `txgnn/results/ranking_metrics.tsv`, `coverage.tsv` |

Phylo, which appears in Figure 6, is scored with the other baselines and is not part of this folder.
