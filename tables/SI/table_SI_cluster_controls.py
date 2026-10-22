"""Appendix, "The full control list" -- Table tbl:SI_cluster_controls (and the control-free
co-clustering ratios quoted in the same section).

Inputs:
  tables/inputs/cluster_labels_by_K.tsv   seed-0 k-medoids cluster labels of every MechDist-network
                                          disease at K = 10..80 (step 5); K = 50 is the shipped clustering
  tables/inputs/cluster_control_pairs.csv the 17 control pairs with their ensemble pairwise score (p),
                                          independence flag, and the percentile of p among
                                          co-investigated IxIDN-Matched test-split pairs (pct_ixpos_tt)
  data/08S_eval/ixidn4364_scores.tsv      IxIDN-Matched pairs and pools
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tablib import DATA, INPUTS, to_markdown, write_and_check  # noqa: E402

K_SHIP = "K50"
L = pd.read_csv(INPUTS / "cluster_labels_by_K.tsv", sep="\t").set_index("disease_id")
C = pd.read_csv(INPUTS / "cluster_control_pairs.csv")
C["scorable"] = C.p.notna() & C.a.isin(L.index) & C.b.isin(L.index)


def same(a, b, K):
    return bool(L.at[a, K] == L.at[b, K])


def yesno(x):
    return "Yes" if x else "No"


def above_below(pct):
    return "Above" if pct > 0.5 else "Below"


lal = C[C.tag == "LALD_Wolman-CESD"].iloc[0]
nf1 = C[C.tag == "NF1mut-NF2"].iloc[0]
new = C[C.origin == "new"]
ind = new[new.independent.astype(str).eq("True") & new.scorable]
ind_pos, ind_neg = ind[ind.role == "pos"], ind[ind.role == "neg"]
pos_share = sum(same(r.a, r.b, K_SHIP) for r in ind_pos.itertuples())
neg_sep = sum(not same(r.a, r.b, K_SHIP) for r in ind_neg.itertuples())
auc = roc_auc_score((ind.role == "pos").astype(int), ind.p)

rows = [
    ["LAL-D / CESD", "Known shared mechanism (illustrative; known before the clustering was finalised)",
     yesno(same(lal.a, lal.b, K_SHIP)), above_below(lal.pct_ixpos_tt)],
    ["NF1 (NF1-mutant) / NF2", "Related but mechanistically distinct (illustrative; known before the clustering was finalised)",
     yesno(same(nf1.a, nf1.b, K_SHIP)), above_below(nf1.pct_ixpos_tt)],
    [f"15 additional controls ({(new.role == 'pos').sum()} positive, {(new.role == 'neg').sum()} negative; "
     f"{len(ind_pos)} positive and {len(ind_neg)} negative independent of every selection step)",
     "See main text for identities",
     f"{pos_share}/{len(ind_pos)} positive pairs share a cluster; {neg_sep}/{len(ind_neg)} negative pairs "
     f"separated (AUROC {auc:.2f} overall)", "n/a"],
]
assert len(new) == 15
df = pd.DataFrame(rows, columns=["Control pair", "Relationship", "Shares cluster at K=50?",
                                 "Pairwise score vs. median co-investigated pair"])

# control-free check (Appendix text): co-cluster rate of co-investigated vs matched-negative test-split pairs
sc = pd.read_csv(DATA / "08S_eval" / "ixidn4364_scores.tsv", sep="\t")
t = sc[sc.pool == "test"]
ratio = {}
for K in L.columns:
    s = L.loc[t.a, K].values == L.loc[t.b, K].values
    ratio[int(K[1:])] = s[t.y.values == 1].mean() / s[t.y.values == 0].mean()
nf1_sep_other = sum(not same(nf1.a, nf1.b, K) for K in L.columns if K != K_SHIP)
text = pd.DataFrame([
    ["co-cluster ratio, test-split IxIDN+ vs IxIDN-, K = 50", f"{ratio[50]:.1f}", "8.2"],
    ["minimum co-cluster ratio over K >= 20", f"{min(v for k, v in ratio.items() if k >= 20):.1f}", "5.3"],
    ["independent controls: positive / negative / excluded", f"{len(ind_pos)}/{len(ind_neg)}/{len(new) - len(ind)}", "6/5/4"],
    ["NF1-mutant/NF2: other K values at which separated", f"{nf1_sep_other} of {len(L.columns) - 1}", "13 of the 15 (manuscript)"],
], columns=["quantity", "rebuilt", "manuscript"])

EXPECTED = [
    ["LAL-D / CESD", "Known shared mechanism (illustrative; known before the clustering was finalised)", "Yes", "Above"],
    ["NF1 (NF1-mutant) / NF2", "Related but mechanistically distinct (illustrative; known before the clustering was finalised)", "Yes", "Below"],
    ["15 additional controls (8 positive, 7 negative; 6 positive and 5 negative independent of every selection step)",
     "See main text for identities",
     "4/6 positive pairs share a cluster; 3/5 negative pairs separated (AUROC 0.83 overall)", "n/a"],
]

if __name__ == "__main__":
    write_and_check(df, __file__, "tbl_SI_cluster_controls", EXPECTED, "Summary of cluster-validation controls.")
    here = Path(__file__).resolve().parent
    text.to_csv(here / "tbl_SI_cluster_controls_text_values.csv", index=False)
    (here / "tbl_SI_cluster_controls_text_values.md").write_text(to_markdown(text))
    pd.DataFrame({"K": list(ratio), "cocluster_ratio_test": list(ratio.values())}).to_csv(
        here / "tbl_SI_cluster_controls_ratio_by_K.csv", index=False)
    print(text.to_string(index=False))
