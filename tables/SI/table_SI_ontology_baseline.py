"""Appendix, "Ontology-only classifier and concordance checks" -- Table tbl:SI_ontology_baseline
(HPO rows).

Built here: the two HPO-only-classifier rows.
  * DisAtlas AUROC on the HPO-scorable pairs (both diseases carry HPO annotations in
    data/00_raw/associations/disease_symptom.parquet) is recomputed from data/08S_eval/ixidn4364_scores.tsv.
  * HPO similarity AUROCs (Resnik-BMA, Lin-BMA, SimGIC, MechDist HPO distance) and the
    DisAtlas-minus-HPO deltas with their bootstrap CIs are read from tables/inputs/hpo_baselines.csv
    (stored baseline output; no HPO ontology file is shipped to recompute the similarities).
Not built: the least-concordant-subset row, the outside-Orphanet row and the Spearman sentence of the
caption all need ORDO-derived quantities, which are out of scope for this repository (tables/README.md).
"""
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import roc_auc_score

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tablib import DATA, INPUTS, RAW, f3, write_and_check  # noqa: E402

sc = pd.read_csv(DATA / "08S_eval" / "ixidn4364_scores.tsv", sep="\t")
hpo_dis = set(pd.read_parquet(RAW / "associations" / "disease_symptom.parquet", columns=["disease_id"]).disease_id)
m = sc.a.isin(hpo_dis) & sc.b.isin(hpo_dis)
auc_model = roc_auc_score(sc.y[m], sc.p[m])

b = pd.read_csv(INPUTS / "hpo_baselines.csv")
b = b[b.scope == "all"]
assert (b.n_scorable == m.sum()).all() and (abs(b.model_auroc_same_pairs - auc_model) < 5e-5).all()
cov = m.mean()
assert abs(100 * cov - b.coverage_pct.iloc[0]) < 0.05
sig = ((b.d_ci_lo > 0) | (b.d_ci_hi < 0)).any()


def rng(lo, hi, nd=2):
    return f"{lo:.{nd}f}-{hi:.{nd}f}"


hpo_rng = rng(b.auroc_scorable.min(), b.auroc_scorable.max())
all_vals = b.auroc_all_unscorable_last.round(2).unique()
assert len(all_vals) == 1
rows = [
    [f"HPO-only classifier, pairs with phenotype annotations on both diseases ({100 * cov:.0f}% of pairs)",
     f"{f3(auc_model)} ({'significantly' if sig else 'not significantly'} higher; "
     f"+{rng(b.delta_model_minus.min(), b.delta_model_minus.max())})",
     f"{hpo_rng} (HPO)"],
    ["HPO-only classifier, all pairs", "--",
     f"{all_vals[0]:.2f} (HPO, undefined for {100 * (1 - cov):.0f}% of pairs)"],
]
df = pd.DataFrame(rows, columns=["Check", "DisAtlas AUROC", "Ontology-only baseline"])

EXPECTED = [
    ["HPO-only classifier, pairs with phenotype annotations on both diseases (35% of pairs)",
     "0.69-0.72 (not significantly higher)", "0.69-0.72 (HPO)"],
    ["HPO-only classifier, all pairs", "--", "0.54 (HPO, undefined for 65% of pairs)"],
]
# Manuscript table repeats the HPO baseline range in the DisAtlas column; the Appendix prose
# ("DisAtlas numerically higher by 0.03-0.06 but not significantly so") matches the rebuilt cell.
KNOWN = [(0, "DisAtlas AUROC")]

if __name__ == "__main__":
    write_and_check(df, __file__, "tbl_SI_ontology_baseline", EXPECTED,
                    "Ontology-only (HPO) baselines on IxIDN-Matched. The ORDO-dependent rows of the "
                    "manuscript table are not rebuilt in this repository.", known=KNOWN)
