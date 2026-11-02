"""Supplementary Table SI_ontology_baseline: ontology-only baselines and concordance checks on IxIDN-Matched.

Summarises the per-baseline results into the ranges quoted in the table:
  - four HPO similarity measures (Resnik-BMA, Lin-BMA, SimGIC, direct MechDist HPO distance), scored on the pairs with
    phenotype annotations on both diseases, against DisAtlas on the same pairs (all pairs and held-out test-test pairs);
  - the same HPO classifiers on all pairs, with unscorable pairs ranked last;
  - the least-concordant subsets where ORDO and HPO distances disagree most;
  - the outside-Orphanet control (at least one disease without an ORDO mapping).

Inputs: data/derived/ontology_baseline/ontology_{baselines,correlations,discordant,orphanet}.csv,
        data/derived/ixidn_robustness/ixidn4364_{per_pair_scores.tsv.gz,ontology_distance.tsv}
Output: table_ontology_baseline.md and table_ontology_baseline_values.json next to this script
"""
import json
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DV = HERE.parents[1] / "data/derived/ontology_baseline"
B = pd.read_csv(DV / "ontology_baselines.csv")
C = pd.read_csv(DV / "ontology_correlations.csv")
D = pd.read_csv(DV / "ontology_discordant.csv")
O = pd.read_csv(DV / "ontology_orphanet.csv").set_index("subset")

HPO = ["hpo_resnik_bma", "hpo_lin_bma", "hpo_simgic", "mechdist_hpo"]
r3 = lambda x: str(Decimal(str(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))  # inputs are stored at 4 dp
rng = lambda s, k=2: f"{s.min():.{k}f}-{s.max():.{k}f}"

h_all = B[(B.scope == "all") & B.column.isin(HPO)]
h_tst = B[(B.scope == "test-test") & B.column.isin(HPO)]
assert h_all.model_auroc_same_pairs.nunique() == 1 and h_tst.model_auroc_same_pairs.nunique() == 1
cov = h_all.coverage_pct.iloc[0]
v = {
    "hpo_coverage_pct": round(float(cov)),
    "model_hpo_pairs": round(float(h_all.model_auroc_same_pairs.iloc[0]), 2),
    "delta_hpo_pairs": "+" + rng(h_all.delta_model_minus),
    "n_ci_excl_zero": int((h_all.d_ci_lo > 0).sum()),
    "model_hpo_pairs_test": round(float(h_tst.model_auroc_same_pairs.iloc[0]), 2),
    "delta_hpo_pairs_test": "+" + rng(h_tst.delta_model_minus),
    "n_ci_excl_zero_test": int((h_tst.d_ci_lo > 0).sum()),
    "hpo_auroc": rng(h_all.auroc_scorable),
    "hpo_all_pairs": f"{h_all.auroc_all_unscorable_last.mean():.2f}",
    "hpo_all_pairs_range": rng(h_all.auroc_all_unscorable_last),
    "hpo_undefined_pct": round(100 - float(cov)),
}
disc = D[D.subset.str.startswith("least concordant")]
v["discordant_model"] = rng(disc.model)
v["discordant_ordo"] = rng(disc[disc.ontology == "neg_ordo_dist"].ontology_auroc)
v["discordant_hpo"] = rng(disc.hpo_auroc)
out = O.loc["at least one disease outside ORDO"]; out_t = O.loc["outside ORDO, held-out test-test"]
v["outside_orphanet_pct"] = round(100 * out.n / O.loc["all pairs"].n)
# Point AUROCs for the outside-Orphanet control are recomputed from the per-pair scores: the csv stores them at 4 dp,
# and rounding 0.9095 again would print 0.910 where the full-precision value 0.90949 gives 0.909.
from sklearn.metrics import roc_auc_score
_P = pd.read_csv(HERE.parents[1] / "data/derived/ixidn_robustness/ixidn4364_per_pair_scores.tsv.gz", sep="\t").merge(
    pd.read_csv(HERE.parents[1] / "data/derived/ixidn_robustness/ixidn4364_ontology_distance.tsv", sep="\t"), on="blind_id")
_out = ~(_P.ordo_dist >= 0); _tt = (_P.split_a == "test") & (_P.split_b == "test")
assert (int(_out.sum()), int((_out & _tt).sum())) == (out.n, out_t.n)
v["outside_orphanet_model"] = r3(roc_auc_score(_P.y[_out], _P.p[_out]))
v["outside_orphanet_model_test"] = r3(roc_auc_score(_P.y[_out & _tt], _P.p[_out & _tt]))
oh = C[(C.x == "neg_ordo_dist") & C.y.isin(["mechdist_hpo", "hpo_simgic"]) & (C.subset == "pooled")].spearman_rho
v["spearman_ordo_hpo"] = rng(oh)
assert float(v["hpo_all_pairs_range"].split("-")[0]) >= 0.535 and float(v["hpo_all_pairs_range"].split("-")[1]) < 0.545
assert disc[disc.ontology == "neg_ordo_dist"].ontology_auroc.between(0.40, 0.56).all()  # reported as "chance"

rows = [
    ("HPO-only classifier, pairs with phenotype annotations on both diseases "
     f"({v['hpo_coverage_pct']}% of pairs)",
     f"{v['model_hpo_pairs']:.2f} ({v['delta_hpo_pairs']}; 95% CI excludes zero for {v['n_ci_excl_zero']} of 4 HPO measures); "
     f"{v['model_hpo_pairs_test']:.2f} test-split ({v['delta_hpo_pairs_test']}; {v['n_ci_excl_zero_test']} of 4)",
     f"{v['hpo_auroc']} (HPO)"),
    ("HPO-only classifier, all pairs", "--", f"{v['hpo_all_pairs']} (HPO, undefined for {v['hpo_undefined_pct']}% of pairs)"),
    ("Least-concordant subset (ORDO vs. HPO disagree most)", v["discordant_model"], f"chance (ORDO); {v['discordant_hpo']} (HPO)"),
    (f"Outside-Orphanet control ({v['outside_orphanet_pct']}% of pairs, >=1 disease outside Orphanet)",
     f"{v['outside_orphanet_model']} ({v['outside_orphanet_model_test']} test-split)", "not defined (ORDO)"),
]
md = ["| Check | DisAtlas AUROC | Ontology-only baseline |", "|:---|:---|:---|"] + ["| " + " | ".join(r) + " |" for r in rows]
md.append(f"\nORDO and HPO distance are correlated at Spearman {v['spearman_ordo_hpo']}.")
(HERE / "table_ontology_baseline.md").write_text("\n".join(md) + "\n")
json.dump(v, open(HERE / "table_ontology_baseline_values.json", "w"), indent=1)
print("\n".join(md))
