#!/usr/bin/env bash
# Rebuild every number, table and figure in the revised manuscript and its appendix from the
# shipped inputs (data/00_raw/, the OBO files and the Zenodo bundle in data/zenodo/).
#
#   ./run_pipeline.sh                 default: no internet, no LLM calls, frozen model
#   ./run_pipeline.sh --full          also rebuild the 13.2M-pair HPO disease-distance table
#   ./run_pipeline.sh --rerun-llm     re-judge IxIDN-Matched candidates with YOUR ANTHROPIC_API_KEY
#                                     (writes ixidn_matched_rerun_llm.tsv; model output is not
#                                     deterministic, so the frozen verdicts stay the reference)
#   ./run_pipeline.sh --retrain       retrain Node2Vec + decoder for the 5 seeds into
#                                     data/zenodo/retrained/ (the frozen model is never overwritten)
#   ./run_pipeline.sh --from STAGE    start at a stage: ontology ixidn llm network model figures tables tests
#
# Set PYTHON to choose the interpreter (default: python).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-python}"
FULL=0; RERUN_LLM=0; RETRAIN=0; FROM=ontology
while [ $# -gt 0 ]; do
  case "$1" in
    --full) FULL=1 ;;
    --rerun-llm) RERUN_LLM=1 ;;
    --retrain) RETRAIN=1 ;;
    --from) FROM="$2"; shift ;;
    -h|--help) sed -n 2,15p "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done
STAGES=(ontology ixidn llm network model figures tables tests)
ENABLED=""; on=0
for s in "${STAGES[@]}"; do [ "$s" = "$FROM" ] && on=1; [ "$on" = 1 ] && ENABLED="$ENABLED $s "; done
[ "$on" = 1 ] || { echo "unknown stage: $FROM" >&2; exit 2; }
want() { case "$ENABLED" in *" $1 "*) return 0 ;; esac; return 1; }

step() { echo; echo "==> $*"; }
in_dir() { local d="$1"; shift; ( cd "$ROOT/$d" && "$PY" "$@" ); }
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg

if want ontology; then
  step "1/8 ontology ancestors, disease distances, ORDO/HPO checks"
  in_dir code/00_ontology build_ancestors.py
  if [ "$FULL" = 1 ]; then in_dir code/00_ontology build_disease_distances.py
  else in_dir code/00_ontology build_disease_distances.py --skip-hpo; fi
  in_dir code/00_ontology build_ontology_dependency.py
  in_dir code/00_ontology build_hpo_baseline.py
  in_dir code/00_ontology build_negative_realism.py
fi

if want ixidn; then
  step "2/8 IxIDN+, IxIDN-, IxIDN-Matched and shortcut baselines"
  if [ "$RERUN_LLM" = 1 ]; then
    : "${ANTHROPIC_API_KEY:?--rerun-llm needs your own ANTHROPIC_API_KEY}"
    PYTHON="$PY" bash "$ROOT/code/01_build_ixidn/run_all.sh" --rerun-llm
  else
    PYTHON="$PY" bash "$ROOT/code/01_build_ixidn/run_all.sh"
  fi
fi

if want llm; then
  step "3/8 LLM annotations: parse stored responses and check the shipped labels"
  in_dir code/02_llm_annotations parse_responses.py
fi

if want network; then
  step "4/8 MechDist distances (S1-S5) and network (S6-S7)"
  in_dir code/03_network mechdist_distances.py
  in_dir code/03_network build_mechdist_network.py
fi

if want model; then
  step "5/8 UPNA-DDAL model"
  if [ "$RETRAIN" = 1 ]; then
    OUT="$ROOT/data/zenodo/retrained"; mkdir -p "$OUT"
    in_dir code/04_train_model 01_strip_benchmark.py --output "$OUT/G_train.pickle"
    for s in 4707 11 202 3407 8191; do
      in_dir code/04_train_model 02_train_node2vec.py --seed "$s" --graph "$OUT/G_train.pickle" --output-dir "$OUT"
      emb=$(ls "$OUT"/*"$s"*.npz | head -1)
      in_dir code/04_train_model 03_train_decoder.py --seed "$s" --embedding "$emb" --graph "$OUT/G_train.pickle" \
        --output "$OUT/decoder_s$s.pt"
    done
  fi
  in_dir code/04_train_model 04_evaluate_ensemble.py
fi

if want figures; then
  step "6/8 main and appendix figures"
  for f in "$ROOT"/figures/main/fig[0-9]*.py "$ROOT"/figures/main/rarity_vs_matched_shortcuts.py "$ROOT"/figures/SI/si_*.py; do
    echo "  - ${f#$ROOT/}"; in_dir "$(dirname "${f#$ROOT/}")" "$(basename "$f")"
  done
fi

if want tables; then
  step "7/8 tables, benchmark values and notebooks (each checks against the manuscript)"
  ( cd "$ROOT" && "$PY" tables/build_all.py && "$PY" benchmarks/verify_benchmarks.py && "$PY" notebooks/execute_notebooks.py )
fi

if want tests; then
  step "8/8 parity tests"
  if [ "$FULL" = 1 ]; then export ONTOLOGY_PARITY_FULL=1; fi
  ( cd "$ROOT" && "$PY" -m pytest -q -p no:cacheprovider -rs tests )
fi
echo; echo "Done: all stages reproduced the manuscript values."
