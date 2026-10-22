"""Score the shipped IxIDN+/IxIDN- benchmark pairs with the frozen 5-seed UPNA-DDAL ensemble.

Thin wrapper around the frozen model's own `frozen_model_loader.FrozenMechDistModel`
(data/zenodo/models/), which verifies every file's sha256 against manifest.json before scoring.
"""
import argparse
import sys

import numpy as np
import pandas as pd


def auroc(y_true, y_score):
    order = np.argsort(-y_score)
    y = np.asarray(y_true)[order]
    n_pos, n_neg = y.sum(), len(y) - y.sum()
    tp, fp = np.cumsum(y), np.cumsum(1 - y)
    return float(np.trapz(tp / n_pos, fp / n_neg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", default="../../data/zenodo/models")
    ap.add_argument("--ixidn-pos", default="../../data/02S_build_IxIDN/ixidn_pos_edges.parquet.gz")
    ap.add_argument("--ixidn-neg", default="../../data/02S_build_IxIDN_neg/ixidn_neg_pinned_edges.parquet.gz")
    args = ap.parse_args()

    sys.path.insert(0, args.model_dir)
    from frozen_model_loader import FrozenMechDistModel  # noqa: E402

    pos = pd.read_parquet(args.ixidn_pos)[["disease_1", "disease_2"]]
    neg = pd.read_parquet(args.ixidn_neg)[["disease_1", "disease_2"]]
    A = pd.concat([pos["disease_1"], neg["disease_1"]]).tolist()
    B = pd.concat([pos["disease_2"], neg["disease_2"]]).tolist()
    y = np.array([1] * len(pos) + [0] * len(neg))

    m = FrozenMechDistModel(args.model_dir)
    p = m.score_pairs(A, B)
    ok = ~np.isnan(p)
    score = auroc(y[ok], p[ok])
    print(f"scored {ok.sum()}/{len(y)} shipped IxIDN+/- pairs (unembedded ids dropped)")
    print(f"ensemble raw AUROC: {score:.4f}")
    print("Reference (test_all_auc_ens, FROZEN_MECHDIST_MODEL.json): 0.7992 -- not directly "
          "comparable: that number is on the held-out TEST split only, this script scores all "
          "shipped pairs with no train/test partition.")


if __name__ == "__main__":
    main()
