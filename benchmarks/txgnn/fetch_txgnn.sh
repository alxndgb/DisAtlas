#!/bin/bash
# Optional. Fetches the public inputs needed to re-run TxGNN inference (run_txgnn_inference.py).
# Not needed to reproduce the paper's numbers: score_txgnn.py reads the saved scores in data/.
#
#   bash benchmarks/txgnn/fetch_txgnn.sh
#
# Downloads into benchmarks/txgnn/external/ (gitignored):
#   TxGNN/                      TxGNN source at the commit used in the paper
#   primekg/{kg,node,edges}.csv PrimeKG release files from Harvard Dataverse
#   primekg/full_graph_42/      TxGNN full_graph split, seed 42 (built locally by TxGNN's TxData)
# The published checkpoint (model.pt, node_emb.pkl, drug_id.pkl, config.pkl) is distributed by the TxGNN
# authors (see the TxGNN README, "pre-trained model"); place it in external/checkpoint/. checksums.md5 lists
# the md5 of every file used in the paper so you can confirm you have the same checkpoint and split.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; X="$HERE/external"; mkdir -p "$X/primekg" "$X/checkpoint"
TXGNN_COMMIT=0e0203c3939f7d12da27aaf10e55daaff5424023
[ -d "$X/TxGNN" ] || git clone https://github.com/mims-harvard/TxGNN.git "$X/TxGNN"
git -C "$X/TxGNN" checkout -q "$TXGNN_COMMIT"
cd "$X/primekg"
for p in "7144484 kg.csv" "7144482 node.csv" "7144483 edges.csv"; do
    set -- $p; [ -s "$2" ] || curl -sfL --retry 3 -o "$2" "https://dataverse.harvard.edu/api/access/datafile/$1"
done
# full_graph split (seed 42), as in TxGNN: TxData(data_folder).prepare_split(split='full_graph', seed=42)
PYTHONPATH="$X/TxGNN" python -c "from txgnn import TxData; TxData('$X/primekg').prepare_split(split='full_graph', seed=42, no_kg=False)"
cd "$HERE" && md5sum -c checksums.md5 || echo "WARNING: some files differ from the ones used in the paper"
