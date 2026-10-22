# Exploratory notebooks

Short, executed notebooks for browsing the shipped data. Each notebook reads only files under
`data/` (the repository root is found by walking up to the folder that contains `data/00_raw/`),
prints small summaries, and `assert`s every manuscript number it recomputes, so a notebook that
runs to the end has reproduced those numbers.

| Notebook | Content | Manuscript sections |
|---|---|---|
| `01_raw_data_tour.ipynb` | Trial triples, diseases, rare/common and cancer labels, drugs, tier-1 disease–gene and disease–symptom sources, propagation pass, gene features (HPA, Pfam, GO), pinned ontologies, LLM curation verdicts | Methods: data sources, disease set, rare/common labels, IxIDN-Matched curation; Supplementary input-source tables |
| `02_ixidn_construction.ipynb` | Co-investigation graph (IxIDN+), distant-pair negatives and veto (IxIDN-), matching universe, LLM-curated positives, matching cells, graph distances, MD5-parity test split | Methods: IxIDN-Matched benchmark, held-out test split; Supplementary: how IxIDN-Matched is built (funnel and shortest-path figure) |
| `03_mechdist_network.ipynb` | MechDist topology (nodes, edges, mean degree, clustering), edges by channel class, candidate pairs, benchmark strip and training graph | Methods: MechDist network, embedding; Supplementary: MechDist network degree distribution |
| `04_headline_evaluation.ipynb` | Headline AUROC on all 4,364 pairs and the 1,094-pair test split, seed SD, rarity strata, pair counts by pool and oncology label | Results: IxIDN-Matched evaluation; Supplementary: seed and sampling uncertainty, performance by rarity stratum (summary table) |

Data read by the notebooks: `data/00_raw/` (raw inputs), the staged benchmark tables in
`data/02S_build_IxIDN/`, `data/02S_build_IxIDN_neg/` and `data/06S_final_benchmark/`, the staged
network in `data/03_network/`, and the frozen ensemble scores and metric tables in `data/08S_eval/`.
Disease-cluster bootstrap intervals are read from those frozen tables, not resampled. Files from
the Zenodo bundle (`data/zenodo/`) are used only if present (OBO checksums in notebook 01, the
candidate-pair count in notebook 03).

To re-execute (about 20 s in total; no Jupyter server needed):

    python notebooks/execute_notebooks.py

The notebooks also open and run in Jupyter or VS Code with the repository environment.
