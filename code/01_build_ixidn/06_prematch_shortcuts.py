"""Step 6 — shortcut signal in the evaluation universe before matching.

Appendix section on IxIDN-Matched shortcut baselines: before matching, among the 27,415 positive and
155,375 negative candidate pairs of the evaluation universe, the sum of the two diseases' trial-graph
degrees separates the classes with AUROC 0.947 and MONDO similarity with 0.758 (pairs without a
SimGIC value dropped). Matching removes both (step 5).

Input:  data/00_raw/frozen/ixidn_matched_universe_features.parquet
Output: data/01_rebuilt/ixidn_universe_prematch_auroc.tsv   score, n, n_pos, n_neg, auroc
"""
import os
import sys

import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ixidn_lib as L  # noqa: E402

REPORTED = {"trial_graph_degree_sum": 0.947, "label_mondo_simgic": 0.758}


def main():
    F = pd.read_parquet(L.MATCH_FEATURES, columns=["label", "deg_a", "deg_b", "mondo_simgic"])
    assert (int((F.label == 1).sum()), int((F.label == 0).sum())) == (27415, 155375)
    scores = {"trial_graph_degree_sum": F.deg_a + F.deg_b, "label_mondo_simgic": F.mondo_simgic}
    rows = []
    for name, s in scores.items():
        m = s.notna()
        y = F.label[m]
        rows.append(dict(score=name, n=int(m.sum()), n_pos=int(y.sum()), n_neg=int((1 - y).sum()),
                         auroc=roc_auc_score(y, s[m])))
    T = pd.DataFrame(rows)
    for r in T.itertuples():
        assert round(r.auroc, 3) == REPORTED[r.score], (r.score, r.auroc)
    T.to_csv(L.out_path("ixidn_universe_prematch_auroc.tsv"), sep="\t", index=False, float_format="%.6f")
    print(T.to_string(index=False))


if __name__ == "__main__":
    main()
