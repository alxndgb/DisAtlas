"""Ontology-only classifier and ORDO/HPO concordance checks on IxIDN-Matched.

Appendix sections: "Dependence between ORDO and HPO" and "Ontology-only classifier and concordance
checks" (Table SI_ontology_baseline).

Inputs (repository-relative):
  data/08S_eval/ixidn4364_scores.tsv                      DisAtlas score p, label y, pool (stored)
  data/00_rebuilt_ontology/ixidn_matched_hpo_scores.csv   HPO baselines (build_hpo_baseline.py)
  data/09S_ontology_baseline/ixidn4364_pairs_scored.csv   MONDO SimGIC label-only control (stored)
  data/00_raw/ontologies/ordo/, mondo.obo, crosswalk      ORDO / MONDO max-LCA distances (max_lca.py)
Outputs (data/00_rebuilt_ontology/):
  ixidn_matched_ordo_mondo_distance.tsv   blind_id, ordo_dist, mondo_dist (-1 = unmapped)
  ontology_dependency_baselines.csv       ontology-only classifiers (scorable pairs; all pairs, unscorable last)
  ontology_dependency_correlations.csv    Spearman correlation of ORDO / MONDO / HPO similarity
  ontology_dependency_discordant.csv      least-concordant ORDO-vs-HPO subsets, no-annotation subsets
  ontology_dependency_orphanet.csv        pairs inside / outside the Orphanet classification
  ontology_dependency_annotation_source.csv  share of MechDist HPO annotations on Orphanet-keyed diseases
                                          (input: data/00_raw/ontologies/ordo/mechdist_hpo_annotation_counts.csv)
CIs: 1000 two-way disease-clustered Poisson bootstrap replicates, seed 7.

    python code/00_ontology/build_ontology_dependency.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import max_lca  # noqa: E402
from ontology import REPO_ROOT  # noqa: E402

OUT = REPO_ROOT / "data" / "00_rebuilt_ontology"
SCORES = REPO_ROOT / "data" / "08S_eval" / "ixidn4364_scores.tsv"
HPO = OUT / "ixidn_matched_hpo_scores.csv"
MONDO_SIMGIC = REPO_ROOT / "data" / "09S_ontology_baseline" / "ixidn4364_pairs_scored.csv"
N = 4364


def ordo_mondo_table(ms):
    od = max_lca.pair_distances(ms.a, ms.b)
    od.insert(0, "blind_id", ms.blind_id.values)
    od.to_csv(OUT / "ixidn_matched_ordo_mondo_distance.tsv", sep="\t", index=False)
    return od


def annotation_source():
    """Share of MechDist disease-HPO annotations that sit on Orphanet-keyed (ORPHA:) diseases."""
    c = pd.read_csv(max_lca.ORDO_DIR / "mechdist_hpo_annotation_counts.csv")
    orpha = c.disease_id.str.startswith("ORPHA:")
    out = pd.DataFrame([dict(n_diseases=len(c), n_hpo_annotations=int(c.n_hpo_annotations.sum()),
                             pct_annotations_on_orphanet_diseases=round(100 * c.n_hpo_annotations[orpha].sum()
                                                                        / c.n_hpo_annotations.sum(), 1),
                             pct_diseases_orphanet_keyed=round(100 * orpha.mean(), 1))])
    out.to_csv(OUT / "ontology_dependency_annotation_source.csv", index=False)
    return out


def auc(y, s):
    if y.min() == y.max():
        return np.nan
    r = rankdata(s); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def main():
    annotation_source()
    ms = pd.read_csv(SCORES, sep="\t")
    od = ordo_mondo_table(ms).replace(-1, np.nan)
    hp = pd.read_csv(HPO)
    mg = pd.read_csv(MONDO_SIMGIC)[["blind_id", "mondo_simgic"]]
    d = (ms[["blind_id", "a", "b", "y", "p", "pool"]]
         .merge(hp[["blind_id", "hpo_simgic", "hpo_resnik_bma", "hpo_lin_bma", "mechdist_hpo"]], on="blind_id")
         .merge(mg, on="blind_id").merge(od, on="blind_id"))
    assert len(d) == N and d.y.sum() == N // 2
    d["neg_ordo_dist"] = -d.ordo_dist; d["neg_mondo_dist"] = -d.mondo_dist   # higher = more similar
    TT = (d.pool == "test").values                                              # held-out test-test pairs
    dis = pd.Index(sorted(set(d.a) | set(d.b))); A = dis.get_indexer(d.a); B = dis.get_indexer(d.b)
    rng = np.random.default_rng(7); W = [rng.poisson(1, len(dis)) for _ in range(1000)]

    def boot(mask, cols, diff=None):
        idx = np.where(mask)[0]; y = d.y.values
        res = {c: auc(y[idx], d[c].values[idx]) for c in cols}
        bs = {c: [] for c in cols}; bd = []
        for w in W:
            sel = np.repeat(idx, (w[A] * w[B])[idx])
            if len(sel) < 20:
                continue
            v = {c: auc(y[sel], d[c].values[sel]) for c in cols}
            for c in cols:
                bs[c].append(v[c])
            if diff:
                bd.append(v[diff[0]] - v[diff[1]])
        out = {c: (res[c], *np.nanpercentile(bs[c], [2.5, 97.5])) for c in cols}
        if diff:
            out["delta"] = (res[diff[0]] - res[diff[1]], *np.nanpercentile(bd, [2.5, 97.5]))
        return out

    # 1. ontology-only classifiers: scorable subset, and all pairs with unscorable ranked last
    base = {"hpo_resnik_bma": "HPO Resnik-BMA (Phenomizer)", "hpo_lin_bma": "HPO Lin-BMA", "hpo_simgic": "HPO SimGIC",
            "mechdist_hpo": "MechDist HPO distance, direct", "neg_ordo_dist": "ORDO max-LCA distance",
            "neg_mondo_dist": "MONDO max-LCA distance", "mondo_simgic": "MONDO SimGIC (label-only control)"}
    rows = []
    for scope, sm in [("all", np.ones(N, bool)), ("test-test", TT)]:
        for c, lab in base.items():
            m = d[c].notna().values & sm
            r = boot(m, [c, "p"], diff=("p", c))
            d["_f"] = d[c].fillna(d[c].min() - 1)          # unscorable ranked last
            rf = boot(sm, ["_f", "p"], diff=("p", "_f"))
            rows.append(dict(scope=scope, baseline=lab, column=c, n_scorable=int(m.sum()), n_pos=int(d.y[m].sum()),
                             coverage_pct=round(100 * m.sum() / sm.sum(), 1),
                             auroc_scorable=r[c][0], ci_lo=r[c][1], ci_hi=r[c][2],
                             model_auroc_same_pairs=r["p"][0], m_ci_lo=r["p"][1], m_ci_hi=r["p"][2],
                             delta_model_minus=r["delta"][0], d_ci_lo=r["delta"][1], d_ci_hi=r["delta"][2],
                             auroc_all_unscorable_last=rf["_f"][0], all_ci_lo=rf["_f"][1], all_ci_hi=rf["_f"][2],
                             model_auroc_all=rf["p"][0], delta_all=rf["delta"][0], dall_ci_lo=rf["delta"][1],
                             dall_ci_hi=rf["delta"][2]))
    pd.DataFrame(rows).round(4).to_csv(OUT / "ontology_dependency_baselines.csv", index=False)

    # 2. ORDO vs HPO-derived similarity (Spearman), pooled and within class
    cr = []
    for x, yv in [("neg_ordo_dist", "mechdist_hpo"), ("neg_ordo_dist", "hpo_simgic"), ("neg_mondo_dist", "mechdist_hpo"),
                  ("neg_mondo_dist", "hpo_simgic"), ("neg_ordo_dist", "neg_mondo_dist"), ("neg_ordo_dist", "p"),
                  ("mechdist_hpo", "p"), ("hpo_simgic", "p")]:
        for grp, g in [("pooled", d), ("IxIDN+", d[d.y == 1]), ("IxIDN-", d[d.y == 0])]:
            g = g[[x, yv]].dropna()
            rho, pv = spearmanr(g[x], g[yv]) if len(g) > 10 else (np.nan, np.nan)
            cr.append(dict(x=x, y=yv, subset=grp, n=len(g), spearman_rho=rho, p_value=pv))
    pd.DataFrame(cr).round(4).to_csv(OUT / "ontology_dependency_correlations.csv", index=False)

    # 3. least-concordant ORDO-vs-HPO subsets, and pairs no ontology classifier can score
    dc = []
    for hcol in ["mechdist_hpo", "hpo_simgic"]:
        for ocol in ["neg_ordo_dist", "neg_mondo_dist"]:
            m = d[[hcol, ocol]].notna().all(1).values; sub = d[m]
            disc = np.abs(rankdata(sub[hcol]) / len(sub) - rankdata(sub[ocol]) / len(sub))
            for q, name in [(0.5, "least concordant half"), (2 / 3, "least concordant third")]:
                mm = np.zeros(N, bool); mm[np.where(m)[0][disc >= np.quantile(disc, q)]] = True
                r = boot(mm, ["p", hcol, ocol], diff=("p", hcol))
                dc.append(dict(hpo=hcol, ontology=ocol, subset=name, n=int(mm.sum()), n_pos=int(d.y[mm].sum()),
                               model=r["p"][0], m_lo=r["p"][1], m_hi=r["p"][2], hpo_auroc=r[hcol][0],
                               ontology_auroc=r[ocol][0], delta_model_minus_hpo=r["delta"][0], d_lo=r["delta"][1],
                               d_hi=r["delta"][2]))
    for name, mm in [("no HPO annotation (HPO classifier undefined)", d.hpo_simgic.isna().values),
                     ("no HPO and no ORDO mapping", (d.hpo_simgic.isna() & d.ordo_dist.isna()).values),
                     ("no HPO, held-out test-test", d.hpo_simgic.isna().values & TT)]:
        r = boot(mm, ["p"])
        dc.append(dict(hpo="-", ontology="-", subset=name, n=int(mm.sum()), n_pos=int(d.y[mm].sum()),
                       model=r["p"][0], m_lo=r["p"][1], m_hi=r["p"][2]))
    pd.DataFrame(dc).round(4).to_csv(OUT / "ontology_dependency_discordant.csv", index=False)

    # 4. pairs inside vs outside the Orphanet classification (stratified-quintile check of the
    #    staging script is not reported in the manuscript and is omitted; it consumes no bootstrap draws)
    orph = d.ordo_dist.notna().values
    oc = []
    for name, mm in [("all pairs", np.ones(N, bool)), ("both diseases in ORDO", orph),
                     ("at least one disease outside ORDO", ~orph), ("outside ORDO, held-out test-test", ~orph & TT)]:
        r = boot(mm, ["p"])
        oc.append(dict(subset=name, n=int(mm.sum()), n_pos=int(d.y[mm].sum()), n_neg=int((1 - d.y[mm]).sum()),
                       model=r["p"][0], m_lo=r["p"][1], m_hi=r["p"][2]))
    pd.DataFrame(oc).round(4).to_csv(OUT / "ontology_dependency_orphanet.csv", index=False)
    print("diseases", len(dis), "test-test pairs", int(TT.sum()), "ORDO-mapped pairs", int(orph.sum()))


if __name__ == "__main__":
    main()
