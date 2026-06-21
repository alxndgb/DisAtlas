![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

![Logo](images/Alexion_AZ_Logo.png)

# Dis2Vec

## Project Description

Dis2Vec (Disease to Vector) is a disease representation learning framework designed to generate biologically grounded disease embeddings from human genetic and phenotypic evidence.

The framework is used for:

1. Disease-Disease Association Learning (DDAL), including transfer-learning-style evaluation against independent positive and negative controls.
2. Unsupervised disease clustering to reveal shared mechanisms across rare and non-rare diseases.

This repository contains cached DDAL predictions and notebooks for querying high-confidence disease-disease links (for example, indication expansion and rare-to-non-rare transfer hypotheses).

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

## data/

Key artifacts used by the notebooks:

1. `dis2vec_all_preds.parquet.gzip`: cached disease-disease association predictions.
2. `dis2vec_clust.xlsx`: clustering metadata used for disease ID lookup and interpretation.

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


