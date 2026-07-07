![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

![Logo](images/Alexion_AZ_Logo.png)

# Dis2Vec

## Project Description

Dis2Vec (Disease to Vector) is a representation learning framework that generates biologically grounded disease embeddings from human genetic and phenotypic data. This repository provides the code and implementation for the Dis2Vec framework introduced in the manuscript "Representation Learning of Human Diseases for Indication Expansion and Investment Decisions."

### Key Features

Dis2Vec addresses a fundamental challenge in translational medicine: computationally modeling complex human diseases to accelerate therapeutic development. The framework introduces a novel strategy that partitions the disease landscape into rare and non-rare categories, enabling systematic knowledge repurposing both within and between these groups.

**Core Capabilities:**

1. **Disease-Disease Association Learning (DDAL)**: Predicts therapeutic transferability using real-world drug repurposing investment decisions from clinical trials as validation benchmarks.

2. **Unsupervised Disease Clustering**: Reveals shared biological mechanisms across diseases, organizing the disease landscape into mechanistically coherent clusters that integrate rare and non-rare disorders.

3. **Cross-Evidence Generalization**: Provides biologically grounded representations validated against independent data modalities not used during training.

**Why Dis2Vec?**

- **Mechanism-Centered**: Unlike traditional knowledge graphs that rely on superficial associations, Dis2Vec learns from the hierarchical structures of Gene Ontology (GO) and Human Phenotype Ontology (HPO) to capture higher-order biological relationships.
- **Addresses Data Scarcity**: Designed to operate effectively in rare disease contexts where data is limited, making it suitable for therapeutic development across the full disease spectrum.
- **Foundation for Agentic AI**: Establishes a biologically grounded disease-representation and benchmarking layer that enables trustworthy agentic biomedical AI systems in rare-disease indication expansion.

This repository contains the implementation code, cached DDAL predictions, and notebooks for querying high-confidence disease-disease links for indication expansion and rare-to-non-rare transfer hypotheses.

## Software Requirements

Python 3.11 is recommended.

Install dependencies from `requirements.txt`:

```bash
python -m pip install -r requirements.txt
```

## How to use/Run the code

Run the query notebook to inspect cached DDAL predictions:

1. Open `notebooks/query_cached_UPNA-DDAL_preds.ipynb`.
2. Execute cells top to bottom.
3. Update disease pair IDs in the query cells to inspect additional hypotheses.

The notebook reads:

- `data/dis2vec_all_preds.parquet.gzip` (cached DDAL predictions)
- `data/dis2vec_clust.xlsx` (disease clustering metadata)

# Data

## Large Data Files

Due to file size constraints, large data files are hosted on Zenodo and are available at:

**[https://zenodo.org/records/21251215](https://zenodo.org/records/21251215)**

### Available Data Files on Zenodo (905.8 MB total)

| Filename | Size | Description |
|----------|------|-------------|
| **dis2vec_all_preds.parquet.gzip** | 244.6 MB | Compressed parquet file with all disease-disease association predictions from UPNA-DDAL model |
| **dis2vec_clust.xlsx** | 3.8 MB | Excel file with disease clustering results, metadata, UMAP coordinates, and cluster labels |
| **go_bp_dist_matrix.parquet.gzip** | 461.8 MB | Compressed distance matrix computed from Gene Ontology (GO) Biological Processes hierarchy using Lowest Common Ancestor (LCA) minimization |
| **hpo_dist_matrix.parquet.gzip** | 183.6 MB | Compressed distance matrix computed from Human Phenotype Ontology (HPO) hierarchy using LCA minimization |
| **clustering_rare_and_common.ipynb** | 1.4 MB | Jupyter notebook demonstrating unsupervised clustering and UMAP visualization of rare and non-rare diseases |
| **query_UPNA-DDAL_preds_cached.ipynb** | 35.6 KB | Jupyter notebook for querying cached DDAL predictions |
| **Dis2Vec_map_dis_vec_clust_UMAP.html** | 5.3 MB | Interactive HTML visualization of disease clustering via UMAP (colored by cluster label) |
| **Dis2Vec_map_rare_UMAP.html** | 5.3 MB | Interactive HTML visualization of disease clustering via UMAP (colored by rare vs. non-rare status) |
| **requirements.txt** | 583 Bytes | Python dependencies file |

**License:** CC BY 4.0

## Local data/ Directory

The `data/` directory in this repository should contain the key artifacts downloaded from Zenodo for use by the notebooks:

1. `dis2vec_all_preds.parquet.gzip`: cached disease-disease association predictions
2. `dis2vec_clust.xlsx`: clustering metadata used for disease ID lookup and interpretation

Download these files from the Zenodo record above and place them in the `data/` directory before running the notebooks.

# Notebooks

## notebooks/

`query_cached_UPNA-DDAL_preds.ipynb` provides examples for querying associations such as:

1. aHUS <-> PNH
2. Rheumatoid Arthritis <-> Antisynthetase Syndrome
3. Rheumatoid Arthritis <-> Polymyositis
4. Rheumatoid Arthritis <-> Dermatomyositis
5. Rheumatoid Arthritis <-> Mixed Connective Tissue Disease
6. Rheumatoid Arthritis <-> Overlap Myositis

# Research Context

This repository aligns with the manuscript:

Representation Learning of Human Diseases for Indication Expansion and Investment Decisions

Dis2Vec is positioned as a biologically grounded representation layer for:

1. Mechanistic disease similarity learning from genetic and phenotypic evidence.
2. Transfer-learning evaluation for DDAL.
3. Cross-disease knowledge repurposing between rare and non-rare disorders.

# Citation

If you use this repository, please cite the associated Dis2Vec manuscript when available and include a link to this repository.

```
@misc{dis2vec,
	title = {Representation Learning of Human Diseases for Indication Expansion and Investment Decisions},
	author = {Ravandi, Cyrus Babak and Mowrey, William R. and Chatterjee, Ayan and Khanshan, Fariba and Haddadi, Parham and Ughetto, Michaël and Barrett, Ian and Ding, Wei and Del Angel, Guillermo and Diethe, Tom and Mobarec, Juan Carlos and Lambden, Simon and Cheng, Seng and Eliassi-Rad, Tina and Risa, Gabriel and Ricchiuto, Piero},
	year = {2026},
	doi = {10.1101/2024.11.19.624381},
	url = {https://www.biorxiv.org/content/10.1101/2024.11.19.624381},
	journal = {In Prep}
}

```


