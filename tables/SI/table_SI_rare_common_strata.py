"""Appendix, "Summary table" (rare/common) -- Table tbl:SI_rare_common_strata.

Point estimates (pairs, AUROC) are recomputed from the frozen DisAtlas ensemble scores:
  data/08S_eval/ixidn4364_scores.tsv           (IxIDN-Matched, column rare_common_pair)
  tables/inputs/ixidn_matched_ccplus_scores.tsv (supplementary common-common benchmark: IxIDN-Matched
                                                 pairs plus 309 added matched negative common-common pairs)
Disease-cluster bootstrap intervals and the comparison with the matched-shortcut stack are read from
the stored bootstrap outputs (not re-drawn here; the bootstrap is stochastic):
  tables/inputs/rare_common_stratified_bootstrap.tsv  95% CI of the IxIDN-Matched strata
  tables/inputs/rare_common_ccplus_bootstrap.tsv      CC+ CI and DisAtlas - matched-shortcut-stack deltas
A delta is "not significant" when its 95% interval includes zero.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tablib import DATA, INPUTS, f3, signed3, thousands, write_and_check  # noqa: E402

sc = pd.read_csv(DATA / "08S_eval" / "ixidn4364_scores.tsv", sep="\t")
cc = pd.read_csv(INPUTS / "ixidn_matched_ccplus_scores.tsv", sep="\t")
bs = pd.read_csv(INPUTS / "rare_common_stratified_bootstrap.tsv", sep="\t")
bs = bs[bs.scope == "All 4364 pairs"].set_index("stratum")
cp = pd.read_csv(INPUTS / "rare_common_ccplus_bootstrap.tsv", sep="\t")
cp = cp[cp.scope == "All pairs"].set_index(["benchmark", "stratum"])


def ci(lo, hi):
    return f"{f3(lo)}–{f3(hi)}"


def vs_shortcut(row):
    if row.delta_vs_stack_lo > 0 or row.delta_vs_stack_hi < 0:
        return (f"{signed3(row.delta_vs_stack)} ({signed3(row.delta_vs_stack_lo)} to "
                f"{signed3(row.delta_vs_stack_hi)})")
    return "not significant"


out = []
for stratum, label in [("rare-rare", "Rare–rare"), ("rare-common", "Rare–common"),
                       ("common-common", "Common–common (IxIDN-Matched)")]:
    g = sc[sc.rare_common_pair == stratum]
    auc = roc_auc_score(g.y, g.p)
    assert np.isclose(auc, bs.at[stratum, "auroc"]) and len(g) == bs.at[stratum, "n_pairs"]
    pairs = thousands(len(g))
    if stratum == "common-common":
        pairs += f" ({int(g.y.sum())} / {int((g.y == 0).sum())})"
        vs = "not interpretable (too few negatives)"
    else:
        vs = vs_shortcut(cp.loc[("IxIDN-Matched", stratum)])
    out.append([label, pairs, f3(auc), ci(bs.at[stratum, "ci_lo"], bs.at[stratum, "ci_hi"]), vs])

g = cc[cc.stratum == "common-common"]
n_added = int((g.source != "IxIDN-Matched").sum())
assert n_added == int((g.label[g.source != "IxIDN-Matched"] == 0).sum())  # added pairs are all negatives
auc = roc_auc_score(g.label, g.disatlas_ens5)
r = cp.loc[("IxIDN-Matched-CC+", "common-common")]
assert np.isclose(auc, r.auroc) and n_added == r.n_added_neg
out.append(["Common–common (supplementary, matched negatives)", f"{n_added} added negatives", f3(auc),
            ci(r.ci_lo, r.ci_hi), vs_shortcut(r)])
df = pd.DataFrame(out, columns=["Stratum", "Pairs", "AUROC", "95% CI", "vs. matched shortcuts"])

EXPECTED = [
    ["Rare–rare", "2,452", "0.722", "0.672–0.772", "+0.121 (+0.059 to +0.179)"],
    ["Rare–common", "1,241", "0.758", "0.668–0.835", "not significant"],
    ["Common–common (IxIDN-Matched)", "671 (645 / 26)", "0.804", "0.655–0.933", "not interpretable (too few negatives)"],
    ["Common–common (supplementary, matched negatives)", "309 added negatives", "0.796", "0.699–0.882", "not significant"],
]

if __name__ == "__main__":
    write_and_check(df, __file__, "tbl_SI_rare_common_strata", EXPECTED,
                    "Stratified AUROC of DisAtlas on IxIDN-Matched by disease-rarity pairing, and on the "
                    "supplementary common-common benchmark with 309 added matched negative pairs. CIs are 95% "
                    "disease-cluster bootstrap intervals.")
