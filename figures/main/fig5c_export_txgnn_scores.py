"""Figure 5 panel c: TxGNN pair scores on the 4,364 IxIDN-Matched pairs, with the protocol of the cross-system figure.

Source: data/derived/fig6_comparison/per_run_scores.tsv (the per-run score table read by fig6_method_comparison.py).
Protocol (fig6_method_comparison.py): a system's pair score is its mean over runs (NaN runs skipped); pairs the system
could not score are ranked last, score = min over all 4,364 pairs - 1e-9.
Writes data/derived/fig5_clustering/fig5c_txgnn_pair_scores.tsv (blind_id, score; read by
fig5_disease_clustering_map.py) and fig5c_comparator_summary.tsv next to this script (or in OUTDIR).
    python figures/main/fig5c_export_txgnn_scores.py [OUTDIR]
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import rankdata
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DV = os.path.join(ROOT, "data", "derived")
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
SYSTEM = "TxGNN published ckpt"   # per-run column prefix of the TxGNN published checkpoint

R = pd.read_csv(os.path.join(DV, "fig6_comparison", "per_run_scores.tsv"), sep="\t")
F = pd.read_csv(os.path.join(DV, "fig5_clustering", "ixidn_matched_pair_scores.tsv"), sep="\t")
assert len(R) == 4364 == len(F) and R.blind_id.is_unique and set(R.blind_id) == set(F.blind_id)
M = F[["blind_id", "y"]].merge(R.drop(columns=["a", "b"]), on="blind_id", how="left", validate="one_to_one", suffixes=("", "_r"))
assert (M.y.values == M.y_r.values).all(), "label mismatch between the Figure 5 and Figure 6 pair tables"


def auc(y, s):
    r = rankdata(s); n1 = y.sum(); n0 = len(y) - n1; return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


y = M.y.values.astype(int)
runs = [c for c in M.columns if c.startswith(SYSTEM + "|")]; assert runs, SYSTEM
v = M[runs].mean(axis=1).values
sc = np.where(np.isnan(v), np.nanmin(v) - 1e-9, v)
pd.DataFrame({"blind_id": M.blind_id, "score": sc}).to_csv(
    os.path.join(DV, "fig5_clustering", "fig5c_txgnn_pair_scores.tsv"), sep="\t", index=False, float_format="%.10g")
S = pd.DataFrame([dict(key="txgnn", system=SYSTEM, n_runs=len(runs), n_pairs=len(v), n_unscored=int(np.isnan(v).sum()),
                       auroc_all=auc(y, sc)),
                  dict(key="disatlas_fig6", system="DisAtlas", n_runs=1, n_pairs=len(M), n_unscored=int(M["DisAtlas"].isna().sum()),
                       auroc_all=auc(y, M["DisAtlas"].values))])
S.to_csv(os.path.join(OUT, "fig5c_comparator_summary.tsv"), sep="\t", index=False, float_format="%.4f"); print(S.to_string())
