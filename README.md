# DisAtlas

Code and data to reproduce every result in the revised manuscript and its Supplementary Information:

> **A Mechanistic Disease Representation Framework for Therapeutic Knowledge Transfer Across Rare and Common Conditions**
> Ravandi, Mowrey, Chatterjee, Khanshan, Haddadi, Ughetto, Barrett, Ding, Del Angel, Diethe, Mobarec, Lambden, Cheng, Eliassi-Rad, Risa, Ricchiuto.

DisAtlas builds a mechanism-distance network between diseases (MechDist) from gene and phenotype evidence, embeds it
with Node2Vec, and trains a pair decoder (UPNA-DDAL) to rank disease pairs for therapeutic knowledge transfer. It is
evaluated on IxIDN+ / IxIDN- (disease pairs that share, or never share, an investigational drug) and on the
degree-matched IxIDN-Matched benchmark (4,364 pairs).

The default run makes no internet requests and no language-model calls, and uses the frozen model. Every stage checks
its output against the value printed in the manuscript and stops if any value differs.

## Quick start

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# download the Zenodo bundle and unpack it into data/zenodo/ (see below)
./run_pipeline.sh
```

| Option | Effect |
|---|---|
| *(none)* | rebuild all numbers, tables and figures from shipped inputs, then run the parity tests |
| `--full` | also rebuild the 13.2M-pair HPO disease-distance table (slow; checked by `tests/test_ontology_parity.py`) |
| `--rerun-llm` | re-judge IxIDN-Matched candidates with your own `ANTHROPIC_API_KEY`; written next to the frozen verdicts, which stay the reference because model output is not deterministic |
| `--retrain` | retrain Node2Vec and the decoder for the five seeds into `data/zenodo/retrained/`; the frozen model is never overwritten |
| `--from STAGE` | start at `ontology`, `ixidn`, `llm`, `network`, `model`, `figures`, `tables` or `tests` |

Set `PYTHON` to choose the interpreter. The pipeline was run with Python 3.11 and the pinned versions in
`requirements.txt`.

## Pipeline

| Stage | Code | Reproduces |
|---|---|---|
| 1 ontology | `code/00_ontology/` | MONDO/HPO/GO ancestors, disease distances, ORDO mapping, ontology baselines |
| 2 benchmarks | `code/01_build_ixidn/` | IxIDN+, IxIDN-, IxIDN-Matched (4,364 pairs) and the matched shortcut baselines |
| 3 LLM annotations | `code/02_llm_annotations/` | disease labels and judge verdicts, parsed from stored model responses |
| 4 network | `code/03_network/` | MechDist distances and network |
| 5 model | `code/04_train_model/` | benchmark strip, Node2Vec, UPNA-DDAL decoder, five-seed ensemble evaluation (AUROC 0.765 all pairs; 0.782 held-out test split) |
| 6 figures | `figures/main/`, `figures/SI/` | main-text and Supplementary figures |
| 7 tables | `tables/`, `benchmarks/`, `notebooks/` | all tables; comparisons with language models and TxGNN (Figure 6); executed exploration notebooks |
| 8 tests | `tests/` | parity of the rebuilt benchmarks, network and ontology tables with the frozen production files |

Each folder has a README mapping its scripts to manuscript sections, figures and tables.

## Data

- `data/00_raw/` holds minimal inputs as parquet files: diseases, genes, associations, drugs, pseudonymised trial
  triples, ontologies, stored LLM responses and frozen benchmark inputs. Licensed sources are pseudonymised or reduced
  to the fields the analysis needs; each subfolder README lists its sources and columns.
- `data/01_*` to `data/08S_*` hold intermediate and evaluation tables that are small enough for git.
- `data/zenodo/` (not in git) holds the large files: the OBO ontologies (`00_raw/ontologies/`), the MechDist network
  and distance tables (`network/`), and the frozen five-seed model with its embeddings, decoders and evaluation logits
  (`models/`). Checksums are in `data/zenodo/models/manifest.json` (model package and graph) and `data/ZENODO_SHA256SUMS.txt` (every file in the bundle; check with `cd data/zenodo && sha256sum -c ../ZENODO_SHA256SUMS.txt`). The Zenodo record is listed below.

## Citation

```bibtex
@article{ravandi2026disatlas,
  title  = {A Mechanistic Disease Representation Framework for Therapeutic Knowledge Transfer Across Rare and Common Conditions},
  author = {Ravandi, Cyrus Babak and Mowrey, William R. and Chatterjee, Ayan and Khanshan, Fariba and
            Haddadi, Parham and Ughetto, Micha\"{e}l and Barrett, Ian and Ding, Wei and Del Angel, Guillermo and
            Diethe, Tom and Mobarec, Juan Carlos and Lambden, Simon and Cheng, Seng and Eliassi-Rad, Tina and
            Risa, Gabriel and Ricchiuto, Piero},
  year   = {2026},
  note   = {Code: \url{https://github.com/alxndgb/DisAtlas}}
}
```

- Code: <https://github.com/alxndgb/DisAtlas>
- Data: Zenodo record (link added on deposit)

## Contact

- Cyrus Babak Ravandi: Cyrus.Ravandi@Alexion.com
- Gabriel Risa: Gabriel.Risa@Alexion.com

## Acknowledgments

The authors acknowledge AstraZeneca's Scientific Computing Platform (SCP) for the computational infrastructure used
in this research.

## License

The code is provided for research and educational purposes. Data files on Zenodo are licensed under CC BY 4.0;
third-party source data keep their original licences (see `data/LICENSES.md` and the `data/00_raw/` READMEs).
