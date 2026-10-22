# 01_build_ixidn — IxIDN+, IxIDN-, IxIDN-Matched

Rebuilds the clinical-trial co-investigation benchmark from the raw tables in `data/00_raw/` and
writes everything to `data/01_rebuilt/` (staged data is never overwritten).

```
PYTHON=python bash code/01_build_ixidn/run_all.sh     # about 2.5 minutes on a laptop
python -m pytest tests/test_ixidn_parity.py
```

| Step | Script | Output (`data/01_rebuilt/`) | Manuscript |
|---|---|---|---|
| 1 | `01_build_dd_network.py` | `dd_nodes.parquet`, `dd_survivors.parquet` | Methods, disease–drug trial network |
| 2 | `02_build_ixidn_pos.py` | `ixidn_pos_edges.parquet` (43,112 edges, 1,379 diseases), `disease_coinvestigated_drugs.parquet` | Methods, IxIDN+ |
| 3 | `03_build_ixidn_neg.py` | `ixidn_neg_edges.parquet` (335,557 pairs), `ixidn_eval_pairs.parquet`, `ixidn_train_strip_pairs.parquet` | Methods, IxIDN- and E_IxIDN removal (13,021 diseases, 325,837 training edges) |
| 4 | `04_match_ixidn.py` | `ixidn_matched.tsv` (4,364 pairs) | Methods, IxIDN-Matched |
| 5 | `05_shortcut_baselines.py` | `ixidn_matched_shortcut_baselines.tsv`, `ixidn_matched_shortcut_stack.tsv` | Results and appendix, matched shortcut baselines; rare/common table contrast |

`tests/test_ixidn_parity.py` joins the stored model scores (`data/08S_eval/ixidn4364_scores.tsv`)
to `ixidn_matched.tsv` and checks AUROC 0.765 (all pairs) and 0.782 (held-out test split).

**LLM judge.** Step 4 uses the stored judge verdicts
(`data/00_raw/llm_annotations/ixidn_judge_verdicts.parquet`); positives with verdict STRONG or
PLAUSIBLE are kept. `python 04_match_ixidn.py --rerun-llm` re-judges the candidates through
`code/02_llm_annotations/call_llm.py judge`; only this switch needs an API key (your own), and
because LLM output is non-deterministic the resulting benchmark (`ixidn_matched_rerun_llm.tsv`)
can differ from the published one.

**Shortcut baselines.** Four non-degree scores (MONDO SimGIC, name TF-IDF, causal-gene count, gene
coverage); AUROC 0.500–0.542. SimGIC is undefined for 7 pairs (one disease has no SimGIC value);
those pairs are dropped for its single-score AUROC (0.5418, n = 4,357) and median-filled inside the
stack. The supplementary common–common rows need `data/00_raw/frozen/ccplus_pair_rows.tsv.gz`.

All paths are repository-relative (`ixidn_lib.py`); seeds are fixed, so the default run is deterministic.
