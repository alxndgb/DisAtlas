"""Ontology-only baselines and concordance checks on IxIDN-Matched (Supplementary Table SI_ontology_baseline).

DisAtlas is scored by the embedding cosine (column ``c``: mean of the five per-seed cosines, seeds 4707-4711) from the
gene-anchored score table; no decoder is involved. Every other column is an ontology-only similarity or distance.

Inputs
  data/zenodo/model_outputs/scores/ixidn_matched_cosine_gene_anchored.parquet   blind_id, a, b, y, pool, c
  data/derived/ontology_baseline/ixidn4364_hpo_similarity.tsv                    HPO Resnik/Lin/SimGIC, MechDist HPO, MONDO SimGIC
  data/derived/ixidn_robustness/ixidn4364_ontology_distance.tsv                  ORDO / MONDO max-LCA distance (-1 = undefined)
Outputs (data/derived/ontology_baseline/)
  ontology_baselines.csv  ontology_correlations.csv  ontology_discordant.csv  ontology_orphanet.csv

"All pairs" baselines rank unscorable pairs last. 95% CIs: disease-cluster bootstrap with the same draws as the shared
subgroup bootstrap (subgroups_shared_bootstrap_B2000.tsv): B = 2,000, numpy default_rng(4707), each draw resamples
the diseases with replacement and weights a pair by the product of its two diseases' counts; percentile intervals.
Values are written at 6 decimals so that the table script rounds once.

    python -m disatlas.ontology_baseline [--out DIR]
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import roc_auc_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
N = 4364
B_DRAWS, SEED = 2000, 4707
SCORE = "c"
BASE = {"hpo_resnik_bma": "HPO Resnik-BMA (Phenomizer)", "hpo_lin_bma": "HPO Lin-BMA", "hpo_simgic": "HPO SimGIC",
        "mechdist_hpo": "MechDist HPO distance, direct", "neg_ordo_dist": "ORDO max-LCA distance (ORDON definition)",
        "neg_mondo_dist": "MONDO max-LCA distance (MONDON definition)", "mondo_simgic": "MONDO SimGIC (label-only control)"}


def load(root=ROOT):
    c = pd.read_parquet(os.path.join(root, "data/zenodo/model_outputs/scores/ixidn_matched_cosine_gene_anchored.parquet"),
                        columns=["blind_id", "a", "b", "y", "pool", SCORE])
    h = pd.read_csv(os.path.join(root, "data/derived/ontology_baseline/ixidn4364_hpo_similarity.tsv"), sep="\t")
    o = pd.read_csv(os.path.join(root, "data/derived/ixidn_robustness/ixidn4364_ontology_distance.tsv"), sep="\t").replace(-1, np.nan)
    d = c.merge(h, on="blind_id").merge(o, on="blind_id")
    assert len(d) == N and d.y.sum() == N // 2 and d[SCORE].notna().all()
    d["neg_ordo_dist"] = -d.ordo_dist; d["neg_mondo_dist"] = -d.mondo_dist   # higher = more similar
    return d


class Boot:
    """Disease-cluster bootstrap; one fixed set of draws shared by every subset and score."""

    def __init__(self, d, n_draws=B_DRAWS, seed=SEED):
        dis = pd.Index(sorted(set(d.a) | set(d.b))); self.A = dis.get_indexer(d.a); self.B = dis.get_indexer(d.b)
        rg = np.random.default_rng(seed); n = len(dis)
        self.W = [np.bincount(rg.integers(0, n, n), minlength=n) for _ in range(n_draws)]
        self.d = d

    def __call__(self, mask, cols, diff=None):
        idx = np.where(mask)[0]; y = self.d.y.values[idx]; X = {c: self.d[c].values[idx] for c in cols}
        res = {c: roc_auc_score(y, X[c]) for c in cols}
        bs = {c: [] for c in cols}; bd = []
        for w in self.W:
            ww = (w[self.A[idx]] * w[self.B[idx]]).astype(float); q = ww > 0
            if q.sum() < 20 or y[q].min() == y[q].max(): continue
            v = {c: roc_auc_score(y[q], X[c][q], sample_weight=ww[q]) for c in cols}
            for c in cols: bs[c].append(v[c])
            if diff: bd.append(v[diff[0]] - v[diff[1]])
        out = {c: (res[c], *np.percentile(bs[c], [2.5, 97.5])) for c in cols}
        if diff: out["delta"] = (res[diff[0]] - res[diff[1]], *np.percentile(bd, [2.5, 97.5]))
        return out


def build(d, boot, s=SCORE):
    TT = (d.pool == "test").values                                              # held-out test-test pairs
    rows = []
    for scope, sm in [("all", np.ones(N, bool)), ("test-test", TT)]:
        for c, lab in BASE.items():
            m = d[c].notna().values & sm
            r = boot(m, [c, s], diff=(s, c))
            d["_f"] = d[c].fillna(d[c].min() - 1)                               # unscorable ranked last
            rf = boot(sm, ["_f", s], diff=(s, "_f"))
            rows.append(dict(scope=scope, baseline=lab, column=c, n_scorable=int(m.sum()), n_pos=int(d.y[m].sum()),
                             coverage_pct=round(100 * m.sum() / sm.sum(), 1),
                             auroc_scorable=r[c][0], ci_lo=r[c][1], ci_hi=r[c][2],
                             model_auroc_same_pairs=r[s][0], m_ci_lo=r[s][1], m_ci_hi=r[s][2],
                             delta_model_minus=r["delta"][0], d_ci_lo=r["delta"][1], d_ci_hi=r["delta"][2],
                             auroc_all_unscorable_last=rf["_f"][0], all_ci_lo=rf["_f"][1], all_ci_hi=rf["_f"][2],
                             model_auroc_all=rf[s][0], delta_all=rf["delta"][0], dall_ci_lo=rf["delta"][1], dall_ci_hi=rf["delta"][2]))
    d.drop(columns="_f", inplace=True)
    cr = []
    for x, yv in [("neg_ordo_dist", "mechdist_hpo"), ("neg_ordo_dist", "hpo_simgic"), ("neg_mondo_dist", "mechdist_hpo"),
                  ("neg_mondo_dist", "hpo_simgic"), ("neg_ordo_dist", "neg_mondo_dist"), ("neg_ordo_dist", s), ("mechdist_hpo", s), ("hpo_simgic", s)]:
        for grp, g in [("pooled", d), ("IxIDN+", d[d.y == 1]), ("IxIDN-", d[d.y == 0])]:
            g = g[[x, yv]].dropna(); rho, pv = spearmanr(g[x], g[yv])
            cr.append(dict(x=x, y=yv, subset=grp, n=len(g), spearman_rho=rho, p_value=pv))
    dc = []
    for hcol in ["mechdist_hpo", "hpo_simgic"]:
        for ocol in ["neg_ordo_dist", "neg_mondo_dist"]:
            m = d[[hcol, ocol]].notna().all(1).values; sub = d[m]
            disc = np.abs(rankdata(sub[hcol]) / len(sub) - rankdata(sub[ocol]) / len(sub))
            for qq, name in [(0.5, "least concordant half"), (2 / 3, "least concordant third")]:
                mm = np.zeros(N, bool); mm[np.where(m)[0][disc >= np.quantile(disc, qq)]] = True
                r = boot(mm, [s, hcol, ocol], diff=(s, hcol))
                dc.append(dict(hpo=hcol, ontology=ocol, subset=name, n=int(mm.sum()), n_pos=int(d.y[mm].sum()),
                               model=r[s][0], m_lo=r[s][1], m_hi=r[s][2], hpo_auroc=r[hcol][0], ontology_auroc=r[ocol][0],
                               delta_model_minus_hpo=r["delta"][0], d_lo=r["delta"][1], d_hi=r["delta"][2]))
    for name, mm in [("no HPO annotation (HPO classifier undefined)", d.hpo_simgic.isna().values),
                     ("no HPO and no ORDO mapping", (d.hpo_simgic.isna() & d.ordo_dist.isna()).values),
                     ("no HPO, held-out test-test", d.hpo_simgic.isna().values & TT)]:
        r = boot(mm, [s]); dc.append(dict(hpo="-", ontology="-", subset=name, n=int(mm.sum()), n_pos=int(d.y[mm].sum()),
                                          model=r[s][0], m_lo=r[s][1], m_hi=r[s][2]))
    orph = d.ordo_dist.notna().values; oc = []
    for name, mm in [("all pairs", np.ones(N, bool)), ("both diseases in ORDO", orph), ("at least one disease outside ORDO", ~orph),
                     ("outside ORDO, held-out test-test", ~orph & TT)]:
        r = boot(mm, [s]); oc.append(dict(subset=name, n=int(mm.sum()), n_pos=int(d.y[mm].sum()), n_neg=int((1 - d.y[mm]).sum()),
                                          model=r[s][0], m_lo=r[s][1], m_hi=r[s][2]))
    return {"ontology_baselines.csv": pd.DataFrame(rows), "ontology_correlations.csv": pd.DataFrame(cr),
            "ontology_discordant.csv": pd.DataFrame(dc), "ontology_orphanet.csv": pd.DataFrame(oc)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=os.path.join(ROOT, "data/derived/ontology_baseline"))
    a = ap.parse_args(); d = load(); os.makedirs(a.out, exist_ok=True)
    for f, t in build(d, Boot(d)).items():
        t.round(6).to_csv(os.path.join(a.out, f), index=False); print(f, len(t))


if __name__ == "__main__":
    main()
