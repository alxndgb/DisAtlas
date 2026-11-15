# training_graph_edges.parquet (new raw input, flagged for author confirmation)

Training graph of the DisAtlas v3.2 embedding (graph the node2vec model was trained on), benchmark (IxIDN-Matched) and FDA
pairs already stripped; 488,103 undirected edges, 13,021 diseases. Columns: u, v (MONDO ids, u < v), in_broad, in_spec
(membership of the broad / specific channel), w_arm (edge weight used for training). No dates, no trial data.
Source: authors' training-graph file (training_graph_edges.parquet of the v3.2 run), columns u, v, in_broad, in_spec, w_arm kept;
edge weights verified identical to the source. sha256 of this file: 5201c6b9efc63962a88a3cfc9b444d6543cf5da7a5ce13dd0c755f50da6705e6.
Licence: CC BY 4.0 for the derived salted edges only (author decision 2026-11-15; author to confirm; data-owner clearance still pending). Derived from the disease-disease evidence graph; no third-party text/records.
Read by: src/disatlas/rebuild/degree_controls.py (degree of each disease for the degree-matched subsets).
Note: this is the graph behind data/derived/training_graph/training_graph_stats.json (stats there are of the unstripped graph).
