#!/usr/bin/env bash
# Rebuild IxIDN+, IxIDN-, IxIDN-Matched and the matched shortcut baselines into data/01_rebuilt/.
# Uses the stored LLM judge verdicts; pass --rerun-llm to step 4 to re-judge with your own key.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python}"
"$PY" 01_build_dd_network.py
"$PY" 02_build_ixidn_pos.py
"$PY" 03_build_ixidn_neg.py
"$PY" 04_match_ixidn.py "$@"
"$PY" 05_shortcut_baselines.py
