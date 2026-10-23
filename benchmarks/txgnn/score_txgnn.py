"""Score the published TxGNN checkpoint on IxIDN-Matched from the saved pair scores used in the paper.

Manuscript: Table "protocols" (Protocol 4); Results, "Explainability contrast with a knowledge graph
foundation model" (AUROC 0.566); Figure 6a, 6b, 6d.

    python benchmarks/txgnn/score_txgnn.py

Reads only benchmarks/txgnn/data/txgnn_pair_scores.tsv.gz and data/08S_eval/ixidn4364_scores.tsv.
Pairs the checkpoint cannot map, and pairs whose two diseases collapse onto one TxGNN node, have no
score and are ranked last (213 of 4,364). No network access, no GPU.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from _common import load_key, read_tsv, ranked_last, pair_bootstrap_metrics, write_tsv  # noqa: E402

NAME = "TxGNN (published checkpoint)"


def main():
    k = load_key()
    t = pd.read_csv(os.path.join(HERE, "data", "txgnn_pair_scores.tsv.gz"), sep="\t")
    d = k[["blind_id", "a", "b", "y"]].merge(t, on=["blind_id", "a", "b"], how="left")
    assert len(d) == 4364 and d.mapped.notna().all()
    scored = d.mapped & ~d.collapsed
    raw = d.score.where(scored).values
    cov = pd.DataFrame([dict(label=int(l), total=int((d.y == l).sum()), mapped=int(d.mapped[d.y == l].sum()),
                             collapsed=int(d.collapsed[d.y == l].sum()),
                             ranked_last=int((~scored)[d.y == l].sum())) for l in (0, 1)])
    write_tsv(cov, os.path.join(HERE, "results", "coverage.tsv"))
    y = d.y.values.astype(int)
    MT = pair_bootstrap_metrics(y, {NAME: ranked_last(raw)})
    write_tsv(MT, os.path.join(HERE, "results", "ranking_metrics.tsv"))
    print(cov.to_string(index=False))
    print(MT.round(3).to_string(index=False))
    return d, raw, cov, MT


if __name__ == "__main__":
    main()
