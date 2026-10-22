# MechDist network (edge and node lists)

- `mechdist_network_edges.parquet.gzip`, `mechdist_network_nodes.parquet.gzip` — edge and node lists
  of the MechDist network (13,021 diseases, 332,082 edges), so the network can be inspected without
  the full graph pickle (`../zenodo/network/mechdist_graph_v2.1.0.pickle`).
- The fused candidate-pair tables (`../zenodo/network/dis_dist_all_df*.parquet.gzip`) are the input to
  `code/03_network/build_mechdist_network.py`, which rebuilds the network (rank-mean fusion,
  class-budgeted top-N selection and connectivity repair) under `data/03_rebuilt/`.
