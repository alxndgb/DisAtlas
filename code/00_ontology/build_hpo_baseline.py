"""HPO-only classifier on IxIDN-Matched (4,364 pairs).

Manuscript: appendix Section SI_ixidn_matched_ontology_baseline, "Ontology-only classifier", and the
two HPO rows of Table SI_ontology_baseline (subset with phenotype annotations on both diseases: AUROC
0.69-0.72, 35% of pairs; all pairs with unscorable pairs ranked last: 0.54).

Inputs : data/00_raw/associations/disease_symptom.parquet, data/00_raw/ontologies/ancestors/hpo_*,
         data/08S_eval/ixidn4364_scores.tsv (pair list, label y, pool, frozen DisAtlas score p)
Outputs: data/00_rebuilt_ontology/ixidn_matched_hpo_scores.csv   per-pair HPO similarities
         data/00_rebuilt_ontology/ixidn_matched_hpo_baseline.csv  AUROCs + 95% CIs
CIs: 1000 two-way disease-clustered Poisson bootstrap replicates, seed 7.
No learned representation and no external service is used.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ontology as O  # noqa: E402
import semsim as SS  # noqa: E402

ROOT = O.REPO_ROOT
OUT = ROOT / "data" / "00_rebuilt_ontology"
PAIRS = ROOT / "data" / "08S_eval" / "ixidn4364_scores.tsv"
SYMPTOMS = ROOT / "data" / "00_raw" / "associations" / "disease_symptom.parquet"
HPO_D_MAX = 25.0
METHODS = {"hpo_resnik_bma": "HPO Resnik-BMA", "hpo_lin_bma": "HPO Lin-BMA",
           "hpo_simgic": "HPO SimGIC", "mechdist_hpo": "MechDist HPO distance, direct"}
N_BOOT, SEED = 1000, 7


def pair_scores(pairs):
    nodes, anc, remap = O.load_ancestors("hpo")
    M = O.distance_matrix(nodes, anc)
    ids = sorted(set(pairs.a) | set(pairs.b))
    S = SS.hpo_scores(ids, pd.read_parquet(SYMPTOMS), anc, remap, M, {t: i for i, t in enumerate(nodes)}, HPO_D_MAX)
    pos = {d: i for i, d in enumerate(ids)}
    ia, ib = pairs.a.map(pos).to_numpy(), pairs.b.map(pos).to_numpy()
    out = pairs[["blind_id"]].copy()
    for c in METHODS:
        out[c] = S[c][ia, ib].astype(np.float64)
    return out


def auc(y, s):
    if y.min() == y.max():
        return np.nan
    r = rankdata(s)
    n1 = y.sum()
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * (len(y) - n1))


def baseline_table(d):
    dis = pd.Index(sorted(set(d.a) | set(d.b)))
    A, B = dis.get_indexer(d.a), dis.get_indexer(d.b)
    rng = np.random.default_rng(SEED)
    W = [rng.poisson(1, len(dis)) for _ in range(N_BOOT)]
    y = d.y.to_numpy()

    def boot(mask, cols, diff):
        idx = np.where(mask)[0]
        res = {c: auc(y[idx], d[c].to_numpy()[idx]) for c in cols}
        bs, bd = {c: [] for c in cols}, []
        for w in W:
            sel = np.repeat(idx, (w[A] * w[B])[idx])
            if len(sel) < 20:
                continue
            v = {c: auc(y[sel], d[c].to_numpy()[sel]) for c in cols}
            for c in cols:
                bs[c].append(v[c])
            bd.append(v[diff[0]] - v[diff[1]])
        out = {c: (res[c], *np.nanpercentile(bs[c], [2.5, 97.5])) for c in cols}
        out["delta"] = (res[diff[0]] - res[diff[1]], *np.nanpercentile(bd, [2.5, 97.5]))
        return out

    rows = []
    for scope, sm in [("all", np.ones(len(d), bool)), ("test-test", (d.pool == "test").to_numpy())]:
        for c, lab in METHODS.items():
            m = d[c].notna().to_numpy() & sm
            r = boot(m, [c, "p"], ("p", c))
            d["_f"] = d[c].fillna(d[c].min() - 1)  # unscorable pairs ranked last
            rf = boot(sm, ["_f", "p"], ("p", "_f"))
            rows.append(dict(
                scope=scope, baseline=lab, column=c, n_scorable=int(m.sum()), n_pos=int(y[m].sum()),
                coverage_pct=round(100 * m.sum() / sm.sum(), 1),
                auroc_scorable=r[c][0], ci_lo=r[c][1], ci_hi=r[c][2],
                model_auroc_same_pairs=r["p"][0], m_ci_lo=r["p"][1], m_ci_hi=r["p"][2],
                delta_model_minus=r["delta"][0], d_ci_lo=r["delta"][1], d_ci_hi=r["delta"][2],
                auroc_all_unscorable_last=rf["_f"][0], all_ci_lo=rf["_f"][1], all_ci_hi=rf["_f"][2],
                model_auroc_all=rf["p"][0], delta_all=rf["delta"][0],
                dall_ci_lo=rf["delta"][1], dall_ci_hi=rf["delta"][2]))
    d.drop(columns="_f", inplace=True)
    return pd.DataFrame(rows).round(4)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(PAIRS, sep="\t")
    assert len(pairs) == 4364 and pairs.y.sum() == 2182
    sc = pair_scores(pairs)
    sc.to_csv(OUT / "ixidn_matched_hpo_scores.csv", index=False)
    tab = baseline_table(pairs[["blind_id", "a", "b", "y", "p", "pool"]].merge(sc, on="blind_id"))
    tab.to_csv(OUT / "ixidn_matched_hpo_baseline.csv", index=False)
    a = tab[tab.scope == "all"]
    print(f"coverage {a.coverage_pct.iloc[0]}% | subset AUROC {a.auroc_scorable.min():.2f}-{a.auroc_scorable.max():.2f}"
          f" | all-pairs AUROC {a.auroc_all_unscorable_last.min():.2f}-{a.auroc_all_unscorable_last.max():.2f}")


if __name__ == "__main__":
    main()
