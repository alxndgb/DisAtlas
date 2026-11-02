#!/usr/bin/env bash
# DisAtlas v3.2: regenerate every manuscript figure, table and number from the shipped data.
#   ./run_pipeline.sh              verify Zenodo bundle, rebuild graph check, figures, tables, parity tests, notebooks
#   ./run_pipeline.sh --retrain    additionally retrain the 5-seed ensemble into retrain_out/ (slow; scores match
#                                  the shipped model statistically, not bit for bit)
#   ./run_pipeline.sh --use-own-llm-key   see README: reruns the stored LLM prompts with your own key (not needed
#                                  for parity; stored outputs in data/raw/benchmarks are used by default)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-.venv/bin/python}
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}" MPLBACKEND=Agg
RETRAIN=0; for a in "$@"; do case $a in --retrain) RETRAIN=1;; --use-own-llm-key) echo "see README (LLM rerun)"; exit 0;; esac; done
stage(){ echo; echo "== $*"; }

stage "1/6 Zenodo bundle checksums"
test -d data/zenodo/model || { echo "data/zenodo missing: download the Zenodo bundle (see data/README.md)"; exit 1; }
SUMS="$PWD/data/ZENODO_SHA256SUMS.txt"; (cd data/zenodo && sha256sum --quiet -c "$SUMS")

stage "2/6 training graph rebuilt from broad and specific channels"
$PY -m disatlas.graph.build --check

stage "3/6 main figures and Table 1"
for f in figures/main/fig*.py tables/main/*.py; do echo "-- $f"; $PY "$f"; done

stage "4/6 SI figures and tables"
for f in figures/SI/si_*.py tables/SI/table_*.py; do echo "-- $f"; $PY "$f"; done

stage "5/6 parity tests"
$PY -m pytest -q tests

stage "6/6 exploratory notebooks"
for nb in notebooks/*.ipynb; do echo "-- $nb"; $PY scripts/run_notebook.py "$nb"; done

if [ $RETRAIN = 1 ]; then stage "retrain (optional)"; $PY -m disatlas.retrain --out retrain_out; fi
echo; echo "DisAtlas pipeline finished"
