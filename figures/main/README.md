# Main-text figures

One script per computed main-text figure. Each script reads only repo-relative inputs, calls no model or
external service, and writes `<stem>.png` (300 dpi) and `<stem>.pdf` next to itself. Rebuilt PNGs match the manuscript
figures: Figure 4 is bit-identical in every pixel; Figure 5 differs only by antialiasing (≤ 2 of 255 intensity levels in
< 0.03% of pixels); Figures 2 and 6 are bit-identical when rendered with matplotlib 3.11.2 (the version the manuscript
figures were drawn with), while other matplotlib versions shift their text layout slightly. Figures 1 and 7 are the
published images.

```
python figures/main/fig2.py   # Figure 2  DisAtlas pipeline overview
python figures/main/fig4.py   # Figure 4  UPNA-DDAL architecture
python figures/main/fig5.py   # Figure 5  disease map and rarity-stratified transfer
python figures/main/fig6.py   # Figure 6  cross-system comparison on IxIDN-Matched
```

Requires numpy, pandas, scipy, matplotlib and pyarrow.

| Figure | Manuscript section | Script | Output | Inputs |
|---|---|---|---|---|
| 1 | Introduction (knowledge transfer between common and rare diseases) | static | `fig1_knowledge_transfer.png` | composed outside Python; shipped as the published image |
| 2 | Methods overview (DisAtlas pipeline) | `fig2.py` | `fig2_pipeline` | schematic; benchmark counts asserted against `data/08S_eval/ixidn4364_scores.tsv` and network size against `data/03_network/mechdist_network_nodes.parquet.gzip`; coverage percentages are the values stated in the text |
| 4 | Methods, UPNA-DDAL | `fig4.py` | `fig4_upna_ddal_architecture` | schematic; panel b is illustrative (only the aHUS-PNH edge is a reported result) |
| 5 | Results, disease map and rarity-stratified transfer | `fig5.py` | `fig5_disease_map`, `fig5c_stratified_auroc.tsv` | `data/fig5_display_map.tsv.gz`, `data/fig5_control_pairs.tsv`, `data/ixidn_matched_go_resnik.tsv.gz`, `data/08S_eval/ixidn4364_scores.tsv`, `data/00_rebuilt_ontology/ixidn_matched_hpo_scores.csv`, `data/00_raw/llm_annotations/disease_labels.parquet` |
| 6 | Results, cross-system comparison | `fig6.py` | `fig6_method_comparison`, `fig6b_ranking_metrics.tsv`, `fig6c_run_agreement.tsv`, `fig6d_top10_jaccard.tsv` | `data/08S_eval/ixidn4364_scores.tsv`, `benchmarks/llm/data/llm_run_scores.tsv.gz`, `benchmarks/llm/data/disatlas_unseen_per_seed.tsv.gz`, `benchmarks/txgnn/data/txgnn_pair_scores.tsv.gz`, `tables/inputs/comparator_run_scores.tsv` |
| 7 | Results, explainability case studies | static | `fig7_explainability.png` | composite of the case-study table and TxGNN explainer views; shipped as the published image |

Shared modules: `disatlas_style.py` (house fonts, sizes, panel letters; strips dates from file metadata) and
`fig_metrics.py` (weighted AUROC equal to scikit-learn `roc_auc_score`, disease-level node bootstrap).

## Checks built into the scripts

- `fig5.py` asserts the panel c values (AUROC, 95% CI, HPO/GO baselines, n+/n−) against the manuscript and
  re-derives each pair's rarity stratum from the stored disease labels. Panel c CI: disease-level bootstrap,
  1,000 draws, seed 4707, pair weight = product of its two diseases' multiplicities. Baseline pairs without a
  score are ranked last.
- `fig6.py` asserts n = 4,364 (2,182 positive) and the node-disjoint AUROC of 0.742 (0.754 on the test split).
  Panel b CIs: 500-draw pair bootstrap, seed 4707. Multi-run methods use the mean over runs; unscored pairs are
  ranked last; ties at a top-k cut-off are shared fractionally.

## Notes

- LLM and TxGNN columns are the stored predictions used in the manuscript. Re-querying the models is a separate,
  opt-in path in `benchmarks/llm/` (`--rerun-llm`, own API key, non-deterministic) and is not used by any figure.
- Disease identifiers that are not MONDO / Orphanet are blinded codes (`DIS_B...`), consistent with `data/00_raw/`.

## Matched shortcut baselines and the rarity contrast

`rarity_vs_matched_shortcuts.py` (Results: matched shortcut baselines and rarity stratification, Figure 5c;
appendix shortcut audit and the rare/common summary table) writes

- `matched_shortcut_baselines.tsv`: single-score AUROC of the four non-degree matched shortcut baselines
  (MONDO SimGIC, name TF-IDF, gene count log-sum, gene coverage minimum) overall and per rarity stratum.
  All pairs: 0.500-0.542. MONDO SimGIC drops the 7 pairs without a MONDO term (4,357 pairs, 0.5418).
- `rarity_vs_matched_shortcuts.tsv`: DisAtlas (full) AUROC, the combined matched-shortcut AUROC (out-of-fold
  logistic model on the four features) and DisAtlas minus combined with a 1,000-draw paired disease-cluster
  bootstrap CI.

Run with `python figures/main/rarity_vs_matched_shortcuts.py` (needs scikit-learn).
