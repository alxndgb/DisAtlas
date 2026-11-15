"""Rebuild the pre-registered event-time summary tables (Figure 7, Supplementary Figures S19-S22, Tables S13-S16) from the staged
per-event percentile table:
    data/derived/si_event_time/{deltas, point_estimates_all, point_estimates_test, merge_and_matching_census}.tsv and
    data/derived/fig7_prospective/{decision_event_time.json, merge_and_matching_census.tsv}.

    python -m disatlas.rebuild.event_time boot  POOL SHARD B BOOTDIR     # one bootstrap shard (pool all|test, SHARD 0-9, B = 200)
    python -m disatlas.rebuild.event_time summary BOOTDIR OUTDIR          # point estimates, deltas with 95% CIs, census, decision

Input: data/zenodo/event_time/event_percentiles.parquet (one row per pool x decision year x drug x disease x method: score, the six
percentile analyses P1 and S1-S5, coverage flag, merge / matching counters; cleared first-trial triplets as (ChEMBL id, MONDO id,
first-trial year); no trial identifiers, no dates) and event_top10_precision.parquet (secondary endpoint). The code is the authors'
pre-registered et03_boot.py / et04_summary.py with the file reading replaced by the staged table; statistics are unchanged:
  * E1 = median over drugs of the drug's median event percentile, E2 = share of drugs whose best event percentile is <= 10;
  * drug-clustered paired bootstrap, cluster = drug over all decision years, RNG numpy default_rng(4707 + 1009 * shard + (pool == 'test')),
    10 shards x 200 = 2,000 replicates per pool (the seed is per shard, so the shard layout must not change);
  * 95% interval = 2.5 / 97.5 percentiles of the paired differences DisAtlas (disatlas_ens5) minus comparator;
  * decision: pooled E1 against weighted Adamic-Adar, checks {P1, S3} x {test, all}, SUPPORTED 4/4, MIXED 1-3/4, NOT_SUPPORTED 0/4.
The 24 yearly slice models x 5 seeds behind the percentiles are not retrained here (percentiles are inputs).
Also rebuilt: abs_percentiles (bootstrap rows of fig7_prospective/abs_percentiles.tsv), coverage.tsv (coverage()) and table_pooled_top10_precision.tsv (top10_precision_table()).
Test: tests/test_rebuild_event_time.py.
"""
import sys, os, glob, json
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
ZEN = ROOT / "data/zenodo/event_time"
AN = ["pct_p1", "pct_s1", "pct_s2", "pct_s3", "pct_s4", "pct_s5"]
REF, DECISIVE = "disatlas_ens5", "adamic_adar_w"
EK = ["y", "compound_id", "event_v3"]
BOOT_USE = ["y", "compound_id", "event_v3", "method", "score", "covered", "n_events_drug_year"] + AN


def load_pool(pool, zen=ZEN):
    X = pd.read_parquet(zen / "event_percentiles.parquet"); X = X[X.pool == pool].drop(columns="pool").reset_index(drop=True)
    PR = pd.read_parquet(zen / "event_top10_precision.parquet"); PR = PR[PR.pool == pool].drop(columns="pool").reset_index(drop=True)
    return X, PR


def subset_flags(X):
    F = X[X.method == X.method.iloc[0]][EK].copy()
    F["excl_2020_2021"] = ~F.y.isin([2020, 2021])
    F = F.merge(X[X.method == X.method.iloc[0]][EK + ["n_events_drug_year"]], on=EK)
    F["single_event"] = F.pop("n_events_drug_year") == 1
    for sm in ("gobp_simgic", "gene_jaccard", "hpo_simgic"):
        if sm in set(X.method):
            F = F.merge(X[X.method == sm][EK + ["covered"]].rename(columns={"covered": f"covered_{sm}"}), on=EK)
    return F


def subsets(F): return [c for c in F.columns if c not in EK]


def wmed(v, w):        # median of v with integer multiplicities w (bootstrap counts), nan rows dropped
    ok = ~np.isnan(v) & (w > 0); v, w = v[ok], w[ok]
    if not len(v): return np.nan
    o = np.argsort(v, kind="stable"); v, w = v[o], w[o]; c = np.cumsum(w); tot = c[-1]
    lo = v[np.searchsorted(c, (tot + 1) // 2 if tot % 2 else tot // 2)]
    hi = v[np.searchsorted(c, tot // 2 + 1)] if tot % 2 == 0 else lo
    return 0.5 * (lo + hi)


def wmean(v, w):
    ok = ~np.isnan(v) & (w > 0); return float((v[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else np.nan


def boot_shard(pool, sh, B=200, zen=ZEN):
    X, PR = load_pool(pool, zen); X = X[BOOT_USE]
    AAF = X[X.method == "adamic_adar"].set_index(EK).score.gt(0).rename("aa_pos")
    X = X.join(AAF, on=EK); assert X.aa_pos.notna().all()
    F = subset_flags(X); X = X.merge(F, on=EK, how="left"); assert X[subsets(F)].notna().all().all()
    meths = sorted(X.method.unique())
    keys = {m: set(map(tuple, X[X.method == m][EK].values)) for m in meths}
    assert all(keys[m] == keys[meths[0]] for m in meths), "event sets differ across methods"
    drugs = np.array(sorted(X.compound_id.unique())); nd = len(drugs); years = sorted(X.y.unique())
    def summ(G, by): return G.groupby(by)[AN].median(), G.groupby(by)[AN].min()
    cube = {}
    for m in meths:
        G = X[X.method == m]
        med, mn = summ(G, "compound_id")
        cube[("pooled", m)] = (med.reindex(drugs).values, (mn.reindex(drugs).values <= 10).astype(float))
        for lab, sel in (("pooled_aa_pos", G.aa_pos.astype(bool)), ("pooled_aa_zero", ~G.aa_pos.astype(bool))):
            med, mn = summ(G[sel], "compound_id"); a = med.reindex(drugs).values; e = (mn.reindex(drugs).values <= 10).astype(float); e[np.isnan(a)] = np.nan
            cube[(lab, m)] = (a, e)
        for lab in subsets(F):
            sel = G[lab].astype(bool); med, mn = summ(G[sel], "compound_id"); a = med.reindex(drugs).values
            e = (mn.reindex(drugs).values <= 10).astype(float); e[np.isnan(a)] = np.nan; cube[("pooled_" + lab, m)] = (a, e)
        for yy in years:
            Gy = G[G.y == yy]; med, mn = summ(Gy, "compound_id")
            a = med.reindex(drugs).values; b = (mn.reindex(drugs).values <= 10).astype(float); b[np.isnan(a)] = np.nan
            cube[(int(yy), m)] = (a, b)
    pc = {}
    for m in meths:
        for sc, lab in (("all", "prec_s4"), (None, "prec_s3")):
            g = PR[(PR.method == m) & ((PR.pop_stratum == "all") if sc else (PR.pop_stratum != "all"))]
            pc[(lab, m)] = g.groupby("compound_id").prec10.mean().reindex(drugs).values
    rng = np.random.default_rng(4707 + 1009 * sh + (pool == "test")); rows = []
    for b in range(B):
        w = np.bincount(rng.integers(0, nd, nd), minlength=nd).astype(float)
        for (scope, m), (A, E) in cube.items():
            for j, an in enumerate(AN):
                rows.append((sh * B + b, scope, an, m, wmed(A[:, j], w), wmean(E[:, j], w)))
        for (lab, m), v in pc.items():
            rows.append((sh * B + b, "pooled", lab, m, wmed(v, w), wmean(v, w)))
    return pd.DataFrame(rows, columns=["rep", "scope", "analysis", "method", "E1", "E2"])


def summarize(boot_dir, zen=ZEN, pools=("test", "all")):
    res, cens, PT = [], [], {}
    for POOL in pools:
        X, PR = load_pool(POOL, zen)
        BT = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{boot_dir}/etboot_{POOL}_sh*.tsv.gz"))]) if boot_dir else None
        if BT is not None: assert BT.rep.nunique() == 2000, (POOL, BT.rep.nunique())
        E = X[X.method == REF]
        for yy, g in [("pooled", E)] + [(int(k), v) for k, v in E.groupby("y")]:
            cens.append(dict(pool=POOL, scope=yy, n_events=len(g), n_drugs=g.compound_id.nunique(),
                             p1_merge0=int((g.p1_merge_level == 0).sum()), p1_merge1=int((g.p1_merge_level == 1).sum()),
                             p1_merge2=int((g.p1_merge_level == 2).sum()), p1_small=int(g.p1_small.sum()),
                             s2_widen=int(g.s2_widen.sum()), s2_short=int(g.s2_short.sum()), n_cand_median=float(g.n_cand.median())))
        pt = {}
        AAF = X[X.method == "adamic_adar"].set_index(EK).score.gt(0).rename("aa_pos")
        X = X.join(AAF, on=EK)
        F = subset_flags(X); X = X.merge(F, on=EK, how="left")
        for m, G in X.groupby("method"):
            for yy, g in ([("pooled", G), ("pooled_aa_pos", G[G.aa_pos]), ("pooled_aa_zero", G[~G.aa_pos])]
                          + [("pooled_" + lab, G[G[lab].astype(bool)]) for lab in subsets(F)]
                          + [(str(int(k)), v) for k, v in G.groupby("y")]):
                med = g.groupby("compound_id")[AN].median(); mn = g.groupby("compound_id")[AN].min()
                for an in AN: pt[(yy, an, m)] = (float(med[an].median()), float((mn[an] <= 10).mean()))
                if yy == "pooled":
                    pt[("pooled", "pct_s2_drawsd_median", m)] = (float(g.sd_s2.median()), np.nan)
        for m in X.method.unique():
            for lab, sel in (("prec_s4", PR.pop_stratum == "all"), ("prec_s3", PR.pop_stratum != "all")):
                v = PR[sel & (PR.method == m)].groupby("compound_id").prec10.mean(); pt[("pooled", lab, m)] = (float(v.median()), float(v.mean()))
        if BT is not None:
            BT["scope"] = BT.scope.astype(str)
            piv = BT.set_index(["rep", "scope", "analysis", "method"])[["E1", "E2"]]; ref = piv.xs(REF, level="method")
            for m in sorted(X.method.unique()):
                if m == REF or m.startswith("disatlas_s"): continue
                dlt = (ref - piv.xs(m, level="method")).groupby(level=["scope", "analysis"])
                q = dlt.quantile([.025, .975]).unstack()
                for (sc, an), r in q.iterrows():
                    a, b = pt.get((sc, an, REF), (np.nan, np.nan)), pt.get((sc, an, m), (np.nan, np.nan))
                    res.append(dict(pool=POOL, scope=sc, analysis=an, comparator=m, disatlas_E1=a[0], comp_E1=b[0], dE1=a[0] - b[0],
                                    dE1_lo=r[("E1", .025)], dE1_hi=r[("E1", .975)], disatlas_E2=a[1], comp_E2=b[1], dE2=a[1] - b[1],
                                    dE2_lo=r[("E2", .025)], dE2_hi=r[("E2", .975)]))
        PT[POOL] = pd.DataFrame([dict(pool=POOL, scope=k[0], analysis=k[1], method=k[2], E1=v[0], E2=v[1]) for k, v in pt.items()])
    return dict(deltas=pd.DataFrame(res), point=PT, census=pd.DataFrame(cens))


def decision(R):
    CHECKS = [(an, p) for an in ("pct_p1", "pct_s3") for p in ("test", "all")]
    d = R[(R.scope == "pooled") & (R.comparator == DECISIVE)].set_index(["analysis", "pool"])
    assert all(k in d.index for k in CHECKS) and not d.index.duplicated().any(), "decision rows missing or duplicated"
    ci = {f"{an}|{p}": [float(d.loc[(an, p), "dE1_lo"]), float(d.loc[(an, p), "dE1_hi"])] for an, p in CHECKS}
    passes = {k: bool(v[1] < 0) for k, v in ci.items()}; n_pass = sum(passes.values())
    verdict = "SUPPORTED" if n_pass == 4 else ("MIXED" if n_pass >= 1 else "NOT_SUPPORTED")
    return dict(rule="pooled E1, DisAtlas_ens5 - weighted AA, 95% CI upper bound < 0, for each of {P1, S3} x {test, all}; "
                     "SUPPORTED = 4/4, MIXED = 1-3/4 (name the passing checks, show all 4 CIs), NOT_SUPPORTED = 0/4",
                ci=ci, passes=passes, n_pass=n_pass, verdict=verdict)


def coverage(zen=ZEN):
    """Coverage of the similarity comparators (share of events whose drug and disease the comparator can score)."""
    X = pd.read_parquet(zen / "event_percentiles.parquet", columns=["pool", "method", "compound_id", "covered"])
    X = X[X.method.isin(["gene_jaccard", "gobp_simgic", "hpo_simgic"])]
    g = X.groupby(["pool", "method"]).agg(n_events=("covered", "size"), n_covered=("covered", "sum"), n_drugs=("compound_id", "nunique")).reset_index()
    g["n_covered"] = g.n_covered.astype(int); g["frac_covered"] = (g.n_covered / g.n_events).round(4)
    return g


COMP_ORDER = ["adamic_adar_w", "adamic_adar", "strength", "degree", "popularity", "gobp_simgic", "gene_jaccard", "hpo_simgic"]


def top10_precision_table(S):
    """Pooled top-10% precision of DisAtlas against each comparator: medians / means over drugs with drug-clustered 95% intervals of the
    paired differences (4 d.p.)."""
    D, f, rows = S["deltas"], (lambda v: f"{v:+.4f}"), []
    D = D[D.scope.astype(str) == "pooled"]
    for pool in ("test", "all"):
        for an in ("prec_s4", "prec_s3"):
            for cp in COMP_ORDER:
                r = D[(D.pool == pool) & (D.analysis == an) & (D.comparator == cp)].iloc[0]
                rows.append(dict(pool=pool, analysis=an, comparator=cp, DisAtlas_median=round(r.disatlas_E1, 4), comp_median=round(r.comp_E1, 4),
                                 dmedian_CI=f"{f(r.dE1)} [{f(r.dE1_lo)}, {f(r.dE1_hi)}]", DisAtlas_mean=round(r.disatlas_E2, 4), comp_mean=round(r.comp_E2, 4),
                                 dmean_CI=f"{f(r.dE2)} [{f(r.dE2_lo)}, {f(r.dE2_hi)}]"))
    return pd.DataFrame(rows)


def abs_percentiles(boot_dir, zen=ZEN, pools=("test", "all")):
    """Absolute median event percentile (E1) with its bootstrap 95% interval for DisAtlas, weighted Adamic-Adar and disease popularity,
    by decision year and pooled (rows with method != chance_expected of fig7_prospective/abs_percentiles.tsv; the chance_expected rows
    need the per-event comparison-set sizes and are not rebuilt)."""
    S = summarize(boot_dir, zen, pools); out = []
    for pool in pools:
        BT = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{boot_dir}/etboot_{pool}_sh*.tsv.gz"))]); BT["scope"] = BT.scope.astype(str)
        q = BT.groupby(["scope", "analysis", "method"]).E1.quantile([.025, .975]).unstack().reset_index().rename(columns={.025: "lo", .975: "hi"})
        P = S["point"][pool].copy(); P["scope"] = P.scope.astype(str)
        d = P.merge(q, on=["scope", "analysis", "method"]); d = d[d.method.isin(["adamic_adar_w", REF, "popularity"]) & d.analysis.isin(["pct_p1", "pct_s3"])]
        d = d[d.scope.isin([str(y) for y in range(1996, 2022)] + ["pooled"])]
        d["n_rep"] = BT.rep.nunique(); out.append(d[["pool", "scope", "analysis", "method", "E1", "lo", "hi", "n_rep"]])
    return pd.concat(out, ignore_index=True)


if __name__ == "__main__":
    if sys.argv[1] == "boot":
        pool, sh, B, out = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]; os.makedirs(out, exist_ok=True)
        T = boot_shard(pool, sh, B); T.to_csv(f"{out}/etboot_{pool}_sh{sh:02d}.tsv.gz", sep="\t", index=False); print(pool, sh, B, len(T))
    else:
        bd, out = sys.argv[2], sys.argv[3]; os.makedirs(out, exist_ok=True); S = summarize(bd)
        for p, t in S["point"].items(): t.to_csv(f"{out}/point_estimates_{p}.tsv", sep="\t", index=False)
        S["deltas"].to_csv(f"{out}/deltas.tsv", sep="\t", index=False); S["census"].to_csv(f"{out}/merge_and_matching_census.tsv", sep="\t", index=False)
        json.dump(decision(S["deltas"]), open(f"{out}/decision_event_time.json", "w"), indent=1)
        coverage().to_csv(f"{out}/coverage.tsv", sep="\t", index=False); top10_precision_table(S).to_csv(f"{out}/table_pooled_top10_precision.tsv", sep="\t", index=False); abs_percentiles(bd).to_csv(f"{out}/abs_percentiles_bootstrap_rows.tsv", sep="\t", index=False); print("wrote", out)
