# Optional retrain: the final training graph (P_full)

The manuscript figures and numbers come from the frozen model outputs in `data/zenodo/model_outputs`, and
`./run_pipeline.sh` reproduces them without retraining. This page covers the optional retrain of the final training graph.

## The final training graph: Habcd broad graph + released specific edges

The training graph of the released model (P_full, 532,958 edges, 13,042 diseases) is the Habcd broad graph
(`data/zenodo/graph/graph_edges.parquet`) plus the released specific-channel edges
(`data/raw/channels/specific/specific_edges.parquet`). Its gene-anchored step used licence-restricted disease-gene pairs,
which are not redistributed, so the build below gives the release rebuild of P_full (532,929 edges; 528,008 with
the broad graph rebuilt by `python -m disatlas.graph.raw_build`):

    python -m disatlas.graph.pfull_build --out build/graph       # exact; checked against tests/pfull_rebuild_signature.json
    python -m disatlas.retrain --graph build/graph/pfull_graph_edges.parquet --out runs/pfull --seeds 4707 4708 4709 4710 4711
    python scripts/pfull_parity.py                               # replay of the authors' result from the shipped cosines

`./run_pipeline.sh --cold` runs the first two steps and the parity test (`tests/test_cold_tolerance.py`). The
authors' 5-seed retrain on P_full against the released model, on the 1,094 IxIDN-Matched test pairs: test AUROC
0.8589 vs 0.8605, ΔAUROC -0.002 [-0.008, +0.005], coverage 100 %: PARITY. The specific channel cannot be rebuilt
from its inputs (licence-restricted literature, pathology, MetaBase, CTD and mechanism-node evidence, about
1.06 million rows); a rebuild from public sources only is FAR (ΔAUROC -0.066). The 8,343 released disease-gene
edges are therefore not sufficient on their own: the retrain needs the released specific edges as well.

Pre-registered parity rule (`scripts/pfull_parity.py`):

| | |
|---|---|
| endpoint | AUROC on the held-out IxIDN-Matched test pairs (1,094) |
| reference | the released final-model embeddings (`data/zenodo/model_outputs/embeddings/gene_anchored`) |
| CI | disease-cluster bootstrap (clusters = first endpoint), 2,000 resamples, 95 % percentile |
| PARITY | CI lower bound of the AUROC difference > -0.02 and coverage >= 95 % |
