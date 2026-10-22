# MechDist network (Methods, mechanism distance; Results, MechDist network)

Rebuilds the MechDist disease network from the raw tables in `data/00_raw/` and the ontology
closures read through `code/00_ontology/ontology.py`. Run from the repo root:

```
python code/03_network/mechdist_distances.py      # S1-S5, about 12 min on a laptop, ~10 GB RAM
python code/03_network/build_mechdist_network.py  # S6-S7, about 2 min
python -m pytest -q tests/test_mechdist_parity.py
```

| Stage | Script | What it does |
|---|---|---|
| S1 | `mechdist_distances.py` | GO distance between diseases sharing at least one gene: per-namespace best-match average of gene-gene minimum term distance / d_max, weighted BP 1.0, MF 0.25, CC 0.1; GO terms annotated to more than 2% of disease genes are masked at distance time (a gene annotated only to masked terms keeps its terms) |
| S2 | `mechdist_distances.py` | HPO distance between diseases sharing at least one phenotype: best-match average of term distance / d_max |
| S3 | `mechdist_distances.py` | cut-off 0.2 on GO and on HPO distance |
| S4 | `mechdist_distances.py` | union of pairs passing either cut-off; proximal channels px_HPA ((1 - cosine)/2 of cell-type expression profiles) and px_DOM (1 - Jaccard of Pfam domain sets) |
| S5 | `mechdist_distances.py`, `crank.py` | CRank rank aggregation (Zitnik et al. 2018; prior 0.05, 100 bins); the GO+HPO CRank sets the per-channel-class edge budget |
| S6 | `build_mechdist_network.py` | rank-mean fusion of the four channels (missing channel ranks last, exact ties keep CRank order), per-channel-class top-200,000, per-disease floor of 3 neighbours |
| S7 | `build_mechdist_network.py` | connectivity repair (up to 2 bridges per island from ranks 200,000-1,500,000, then the full pool); largest connected component |

Output: 13,021 diseases and 332,082 edges (`data/03_rebuilt/`), edge-for-edge identical to the
production network.

Tie order. CRank output has many exact ties, so the stable sort in S5/S6 depends on the row
order of the candidate pairs. Production sorted pairs by the original disease-ID strings, and blinded
IDs sort differently: in testing, blinded-string order changed 16 of 332,082 edges.
`data/03_rebuilt/frozen/disease_sort_rank.parquet` stores each disease's position in the original
string order, and S4 orders pairs by it. With that order every S1-S5 value and the edge set match
production exactly.
