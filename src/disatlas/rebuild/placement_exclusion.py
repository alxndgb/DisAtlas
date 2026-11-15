"""Rebuild data/derived/placement_confidence/exclusion_sensitivity.tsv (Supplementary Figure SI_placement_confidence;
SI placement-confidence paragraph: 105 critical-tier diseases, AUROC 0.881 on their 1,140 pairs, kept AUROC 0.817,
difference -0.017 (95% CI -0.034 to -0.0003)).

    python -m disatlas.rebuild.placement_exclusion OUT.tsv

Scores: five-seed mean cosine; IxIDN-Matched from data/derived/fig6_comparison/per_run_scores.tsv (column DisAtlas),
the two other pair sets from data/zenodo/model_outputs/cosine_matrix. Pair sets: all 4,364 IxIDN-Matched pairs, the
test split (pool == "test"), the MechDist link-prediction test split and the 906-pair tuning split (staged by
scripts/maintainer/stage_placement_confidence.py). Flags: data/derived/placement_confidence/disease_placement_flags.tsv.gz.
Rules: C5 unstable, C6 ill-defined, C5 or C6, Critical tier (>= 3 issue classes, scored diseases), Critical tier
counting class 5 only when a flag other than UNSTABLE carries it.
Per benchmark (authors' s_exclusion_sensitivity.run): diseases sorted; default_rng(0); 1,000 multinomial
disease-cluster draws (pair weight w_a * w_b) for the delta interval, then per rule 500 class-matched random removals.
Test: tests/test_rebuild_placement_exclusion.py.
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
PC = ROOT / "data/derived/placement_confidence"
CM = ROOT / "data/zenodo/model_outputs/cosine_matrix"


def wauc(y, p, w=None):
    y = np.asarray(y); p = np.asarray(p); w = np.ones(len(y)) if w is None else np.asarray(w, float)
    u, inv = np.unique(p, return_inverse=True)
    wp = np.bincount(inv, w * (y == 1), len(u)); wn = np.bincount(inv, w * (y == 0), len(u)); cn = np.cumsum(wn) - wn
    return float(((cn + 0.5 * wn) * wp).sum() / (wp.sum() * wn.sum()))


def run(d, sets, dclass, nboot=1000, nnull=500, seed=0):
    rng = np.random.default_rng(seed)
    dis = np.array(sorted(set(d.a) | set(d.b))); ix = {x: i for i, x in enumerate(dis)}
    ia = d.a.map(ix).values; ib = d.b.map(ix).values; y = d.y.values; p = d.p.values
    full = wauc(y, p)
    W = rng.multinomial(len(dis), np.ones(len(dis)) / len(dis), size=nboot).astype(float)
    PW = W[:, ia] * W[:, ib]
    bfull = np.array([wauc(y, p, w) if w.sum() > 0 else np.nan for w in PW])
    cl = pd.Series(dclass).reindex(dis).fillna("NA").values
    out = []
    for name, S in sets.items():
        exc = np.isin(dis, list(S)); keep = ~(exc[ia] | exc[ib])
        sub = wauc(y[keep], p[keep]); tch = wauc(y[~keep], p[~keep]) if (~keep).sum() > 10 and len(set(y[~keep])) > 1 else np.nan
        bs = np.array([wauc(y[keep], p[keep], w[keep]) if (w[keep] * (y[keep] == 1)).sum() > 0 and (w[keep] * (y[keep] == 0)).sum() > 0 else np.nan for w in PW])
        dl = bs - bfull; nd = []
        for _ in range(nnull):
            m = np.zeros(len(dis), bool)
            for c in np.unique(cl):
                pool = np.where(cl == c)[0]; k = int((exc & (cl == c)).sum())
                if k: m[rng.choice(pool, k, replace=False)] = True
            kk = ~(m[ia] | m[ib]); nd.append(wauc(y[kk], p[kk]) - full if kk.sum() > 0 and len(set(y[kk])) > 1 else np.nan)
        nd = np.array(nd)
        out.append(dict(exclusion=name, n_dis_excluded=int(exc.sum()), n_dis_total=len(dis), n_pairs_kept=int(keep.sum()), n_pairs_total=len(d),
                        pos_kept=int(y[keep].sum()), neg_kept=int((1 - y[keep]).sum()), auroc_full=full, auroc_kept=sub,
                        delta=sub - full, delta_lo=np.nanpercentile(dl, 2.5), delta_hi=np.nanpercentile(dl, 97.5),
                        auroc_touching=tch, null_delta_median=np.nanmedian(nd), null_delta_lo=np.nanpercentile(nd, 2.5),
                        null_delta_hi=np.nanpercentile(nd, 97.5), pct_null_ge_obs=float(np.nanmean(nd >= sub - full))))
    return pd.DataFrame(out)


def rules():
    Q = pd.read_csv(PC / "disease_placement_flags.tsv.gz", sep="\t").set_index("disease_id")
    sc = Q.in_disatlas_scored.astype(bool); c5f = ["flag_UNSTABLE", "flag_REPAIR_DEPENDENT", "flag_LOW_DEGREE", "flag_ANTIHUB"]
    assert (Q[c5f].astype(bool).any(axis=1) == Q.C5_fragile.astype(bool)).all()
    only = Q.C5_fragile.astype(bool) & Q.flag_UNSTABLE.astype(bool) & ~Q[c5f[1:]].astype(bool).any(axis=1)
    crit = set(Q.index[(Q.n_issue_classes >= 3) & sc]); crit_ev = set(Q.index[((Q.n_issue_classes - only.astype(int)) >= 3) & sc])
    assert len(crit) == 892, len(crit)
    sets = {"C5 unstable": set(Q.index[Q.C5_fragile == 1]), "C6 ill-defined": set(Q.index[Q.C6_entity == 1]),
            "C5 or C6": set(Q.index[(Q.C5_fragile == 1) | (Q.C6_entity == 1)]), "Critical tier": crit, "Critical (evidence-only)": crit_ev}
    return sets, Q.disease_class.to_dict()


def cosine(d):
    C = np.load(CM / "cosine_ens5_13021.npy", mmap_mode="r"); ix = {x: i for i, x in enumerate(json.load(open(CM / "ids_13021.json")))}
    d = d.copy(); d["p"] = np.asarray(C[d.a.map(ix).values, d.b.map(ix).values], float); return d


def benchmarks():
    S = pd.read_csv(ROOT / "data/derived/fig6_comparison/per_run_scores.tsv", sep="\t"); assert len(S) == 4364
    S = S.rename(columns={"DisAtlas": "p"})[["a", "b", "y", "pool", "p"]]
    return {"IxIDN-Matched 4,364 (all pools)": S, "IxIDN-Matched test-test (held out)": S[S.pool == "test"],
            "test_dm (held-out MechDist test split)": cosine(pd.read_csv(PC / "test_dm_pairs.tsv.gz", sep="\t")),
            "IxIDN-910 (tune split)": cosine(pd.read_csv(PC / "ixidn910_pairs.tsv", sep="\t"))}


def compute():
    sets, dclass = rules(); out = []
    for nm, d in benchmarks().items():
        T = run(d[["a", "b", "y", "p"]].reset_index(drop=True), sets, dclass); T.insert(0, "benchmark", nm)
        T["pct_pairs_kept"] = 100 * T.n_pairs_kept / T.n_pairs_total; out.append(T)
    return pd.concat(out, ignore_index=True)


if __name__ == "__main__":
    compute().to_csv(sys.argv[1], sep="\t", index=False); print("wrote", sys.argv[1])
