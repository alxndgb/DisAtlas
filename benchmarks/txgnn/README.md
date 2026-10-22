# TxGNN published checkpoint (Protocol 4)

    python benchmarks/txgnn/score_txgnn.py        # default: saved scores only

`data/txgnn_pair_scores.tsv.gz` contains the scores of the published TxGNN checkpoint, used unchanged (not retrained),
on the 4,364 IxIDN-Matched pairs. Each pair has TxGNN disease indices, `mapped`, `collapsed` (both diseases on one
TxGNN node), `score` (TxGNN disease-disease DistMult) and `cos` (embedding cosine, not used in the paper). The scorer
uses `score`. Unmapped and collapsed pairs (213) are ranked last, which gives AUROC 0.566.

Optional recomputation, not needed for the paper's numbers:
1. `bash benchmarks/txgnn/fetch_txgnn.sh` clones TxGNN at the pinned commit, downloads PrimeKG from Harvard
   Dataverse and builds the `full_graph` split (seed 42). Place the published checkpoint from the TxGNN authors in
   `external/checkpoint/`. `checksums.md5` lists the md5 of the checkpoint and split files used in the paper.
2. `python benchmarks/txgnn/run_txgnn_inference.py` writes `rerun/txgnn_pair_scores.tsv.gz` and compares it with
   the saved scores.
