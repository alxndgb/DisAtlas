#!/usr/bin/env bash
# DisAtlas v3.2: regenerate every manuscript figure, table and number.
#   ./run_pipeline.sh                 WARM (default): verify the Zenodo bundle and use its processed data (training
#                                     graph, frozen 5-seed model, scores); figures, tables, parity tests, notebooks
#   ./run_pipeline.sh --cold          COLD: rebuild channels and training graph from data/raw, retrain the 5-seed
#                                     ensemble and re-score the benchmarks into build/, compare build/ with the
#                                     frozen record (exact graph, statistical scores); figures, tables and the
#                                     manuscript numbers always come from the frozen record
#   --cold-qc      opt-in, either mode: in stage qc, additionally rerun the IxIDN-Matched language-model QC
#                  (needs ANTHROPIC_API_KEY and QC_MODEL) and rebuild IxIDN-Matched from the new verdicts into
#                  build/qc/; without it, stage qc rebuilds IxIDN-Matched from the stored verdicts only
#   --from STAGE   start at STAGE     --only STAGE   run STAGE alone
#   stages  warm: verify graph qc figures si tests notebooks
#           cold: verify channels graph qc retrain score figures si tests notebooks
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-.venv/bin/python}
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}" MPLBACKEND=Agg
MODE=warm FROM= ONLY= COLDQC=0
while [ $# -gt 0 ]; do case $1 in
  --cold) MODE=cold;; --cold-qc) COLDQC=1;; --from) FROM=$2; shift;; --only) ONLY=$2; shift;;
  *) echo "unknown option $1"; exit 2;; esac; shift; done
export DISATLAS_DATA_MODE=warm   # figures, tables, notebooks: frozen record
if [ $MODE = cold ]; then STAGES="verify channels graph qc retrain score figures si tests notebooks"
else STAGES="verify graph qc figures si tests notebooks"; fi
for s in $FROM $ONLY; do case " $STAGES " in *" $s "*) ;; *) echo "unknown $MODE stage '$s' ($STAGES)"; exit 2;; esac; done
SEEDS="4707 4708 4709 4710 4711"
B=build                           # cold-mode outputs (gitignored)
STARTED=1; [ -n "$FROM" ] && STARTED=0
run(){  # run stage $1 unless excluded by --from / --only
  [ "$1" = "$FROM" ] && STARTED=1
  [ $STARTED = 1 ] || return 0
  [ -z "$ONLY" ] || [ "$ONLY" = "$1" ] || return 0
  echo; echo "== $MODE: $1"; "stage_$1"; }

stage_verify(){
  if [ ! -d data/zenodo/model ]; then
    [ $MODE = cold ] && { echo "data/zenodo absent: cold build continues; comparisons with the frozen record are skipped"; return 0; }
    echo "data/zenodo missing: download the Zenodo bundle (see data/README.md)"; exit 1; fi
  SUMS="$PWD/data/ZENODO_SHA256SUMS.txt"; (cd data/zenodo && sha256sum --quiet -c "$SUMS"); }

stage_channels(){
  for t in broad specific; do
    $PY -c "import disatlas.channels.$t" 2>/dev/null || {
      echo "disatlas.channels.$t is not available in this checkout; the channel rebuild code ships separately."
      echo "Resume with --cold --from graph to build the training graph from the shipped data/raw/channels."; exit 1; }
    $PY -m disatlas.channels.$t --out $B/channels/$t --cold
  done; }

stage_graph(){
  if [ $MODE = warm ]; then $PY -m disatlas.graph.build --check; return; fi
  CH=$B/channels; CHECK=
  [ -d $CH/broad ] && [ -d $CH/specific ] || { CH=data/raw/channels; echo "no $B/channels: using shipped $CH"; }
  [ -d data/zenodo/graph ] && CHECK=--check
  $PY -m disatlas.graph.build --channels $CH --out $B/graph $CHECK; }

stage_qc(){  # IxIDN-Matched from the stored QC verdicts (exact); --cold-qc also reruns the language-model QC
  $PY -m disatlas.benchmarks.ixidn_matched --out $B/ixidn_matched --check
  [ $COLDQC = 1 ] || return 0
  $PY -m disatlas.benchmarks.matched_qc --out $B/qc --model "${QC_MODEL:?set QC_MODEL to a model identifier}"
  $PY -m disatlas.benchmarks.ixidn_matched --out $B/qc --verdicts $B/qc/llm_curation_verdicts.tsv --check; }

stage_retrain(){  # one process per seed, in parallel
  M=$B/model/disatlas_v3.2_ensemble; mkdir -p $M; T=$(( $(nproc) / 5 )); [ $T -ge 1 ] || T=1
  for s in $SEEDS; do OMP_NUM_THREADS=$T DISATLAS_DATA_MODE=cold $PY -m disatlas.retrain --out $M --seeds $s > $M/retrain_s$s.log 2>&1 & done
  wait; cat $M/retrain_s*.log; DISATLAS_DATA_MODE=cold $PY -m disatlas.retrain --out $M --config-only; }

stage_score(){ DISATLAS_DATA_MODE=cold $PY -m disatlas.rescore --out $B/scores --all-pairs; }

stage_figures(){ for f in figures/main/fig*.py tables/main/*.py; do echo "-- $f"; $PY "$f"; done; }
stage_si(){ for f in figures/SI/si_*.py tables/SI/table_*.py; do echo "-- $f"; $PY "$f"; done; }

stage_tests(){
  if [ $MODE = warm ]; then $PY -m pytest -q tests; return; fi
  [ -d data/zenodo/model ] || { echo "data/zenodo absent: no frozen record to compare with"; return 0; }
  DISATLAS_DATA_MODE=cold $PY -m pytest -q tests/test_cold_exact.py tests/test_cold_tolerance.py  # build/ vs record
  $PY -m pytest -q tests; }                            # manuscript numbers from the frozen record (exact)

stage_notebooks(){ for nb in notebooks/*.ipynb; do echo "-- $nb"; $PY scripts/run_notebook.py "$nb"; done; }

for s in $STAGES; do run $s; done
echo; echo "DisAtlas pipeline finished ($MODE)"
