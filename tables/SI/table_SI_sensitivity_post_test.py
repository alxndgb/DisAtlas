"""Appendix, "Sensitivity to settings fixed after test-split results" -- Table tbl:SI_sensitivity_post_test.

Inputs: tables/inputs/sensitivity_post_test_scores.tsv -- per-seed sigmoid probabilities on the
1,094 test-split IxIDN-Matched pairs for five seeds of each arm:
  reference         DisAtlas as shipped (identical to p_s* in data/08S_eval/ixidn4364_scores.tsv)
  retrained         same settings, embeddings and decoder retrained (run-to-run reference)
  window1/window5   node2vec context window 1 / 5 (embeddings and decoder retrained)
  degree_negatives  decoder retrained with degree-proportional negatives on the shipped embeddings
Ensemble = mean of the five per-seed probabilities. Deltas are differences in per-seed means:
decoder-only rows vs DisAtlas, window rows vs the retrained reference. SD uses ddof=1.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tablib import DATA, INPUTS, f3, signed3, write_and_check  # noqa: E402

SEEDS = [4707, 11, 202, 3407, 8191]
s = pd.read_csv(INPUTS / "sensitivity_post_test_scores.tsv", sep="\t")
assert len(s) == 1094 and s.y.sum() == 491

# the reference arm is the frozen headline model
head = pd.read_csv(DATA / "08S_eval" / "ixidn4364_scores.tsv", sep="\t").set_index("blind_id").loc[s.blind_id]
for k in SEEDS:
    assert np.allclose(head[f"p_s{k}"].values, s[f"reference_s{k}"].values)


def arm(name):
    M = s[[f"{name}_s{k}" for k in SEEDS]]
    auc = [roc_auc_score(s.y, M[c]) for c in M]
    ap = [average_precision_score(s.y, M[c]) for c in M]
    e = M.mean(axis=1)
    return dict(ens_auc=roc_auc_score(s.y, e), ens_ap=average_precision_score(s.y, e),
                auc_mean=np.mean(auc), auc_sd=np.std(auc, ddof=1), ap_mean=np.mean(ap))


A = {a: arm(a) for a in ["reference", "retrained", "window1", "window5", "degree_negatives"]}
ROWS = [("DisAtlas (window 2, uniform decoder negatives, unweighted loss)", "reference", "reference"),
        ("Same settings, retrained (run-to-run reference)", "retrained", "reference"),
        ("Node2vec window 1", "window1", "retrained"),
        ("Node2vec window 5", "window5", "retrained"),
        ("Degree-proportional decoder negatives", "degree_negatives", "reference")]
out = []
for label, a, ref in ROWS:
    r, b = A[a], A[ref]
    out.append([label, f3(r["ens_auc"]), f3(r["ens_ap"]), f"{f3(r['auc_mean'])} ± {f3(r['auc_sd'])}",
                signed3(r["auc_mean"] - b["auc_mean"]), signed3(r["ap_mean"] - b["ap_mean"])])
df = pd.DataFrame(out, columns=["Setting", "Ensemble AUROC", "Ensemble AUPRC",
                                "Per-seed AUROC (mean ± SD)", "ΔAUROC", "ΔAUPRC"])

EXPECTED = [
    ["DisAtlas (window 2, uniform decoder negatives, unweighted loss)", "0.782", "0.772", "0.776 ± 0.003", "+0.000", "+0.000"],
    ["Same settings, retrained (run-to-run reference)", "0.784", "0.772", "0.778 ± 0.007", "+0.002", "+0.001"],
    ["Node2vec window 1", "0.793", "0.778", "0.781 ± 0.007", "+0.004", "+0.002"],
    ["Node2vec window 5", "0.769", "0.761", "0.765 ± 0.005", "-0.012", "-0.009"],
    ["Degree-proportional decoder negatives", "0.784", "0.773", "0.778 ± 0.005", "+0.002", "+0.002"],
]

if __name__ == "__main__":
    write_and_check(df, __file__, "tbl_SI_sensitivity_post_test", EXPECTED,
                    "Sensitivity of test-split performance (1,094 pairs, 491 positive) to the two "
                    "settings fixed after test-split results had been seen.")
