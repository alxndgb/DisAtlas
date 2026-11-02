"""Supplementary Table GO-release: DisAtlas advantage over weighted Adamic-Adar with current and slice-dated GO annotations.

Post hoc robustness check on the 2008 and 2016 slices (decision years 2009 and 2017), common-node event set.
Inputs: data/zenodo/si_bio/go_release_event_percentiles.tsv.gz (event-level matched percentiles of the original models,
weighted Adamic-Adar, and the current-release and dated-release GO models) and data/derived/si_go_release/ (annotation
counts, common-node census, eculizumab case summary). Drug-cluster bootstrap, B = 2000, numpy default_rng(4707), the same
draws for every arm so differences are paired. Positive = DisAtlas better (E1 = weighted AA - DisAtlas median percentile;
E2 = (DisAtlas - weighted AA) share of drugs with the best new trial in the top 10%, x 100).
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, write_numbers  # noqa: E402
import numpy as np, pandas as pd

B, SEED, YEARS = 2000, 4707, (2009, 2017)
AN = {"pct_p1": "Popularity and weighted-degree groups", "pct_s3": "Popularity groups", "pct_s4": "No grouping"}
POOL = {"test": "Test-pool diseases", "all": "All diseases"}
EV = pd.read_csv(ZENODO_SIBIO / "go_release_event_percentiles.tsv.gz", sep="\t")
SRC = DERIVED / "si_go_release"
A = pd.read_csv(SRC / "grm_annotation_counts.tsv", sep="\t").set_index(["arm", "slice"])
C = pd.read_csv(SRC / "grm_common_node_census.tsv", sep="\t").set_index("y")
NMO = pd.read_csv(SRC / "grm_nmo_S2008_summary.tsv", sep="\t").set_index("arm")

def wmed(v, w):
    ok = ~np.isnan(v) & (w > 0); v, w = v[ok], w[ok]
    if not len(v): return np.nan
    o = np.argsort(v, kind="stable"); v, w = v[o], w[o]; c = np.cumsum(w); tot = c[-1]
    lo = v[np.searchsorted(c, (tot + 1) // 2 if tot % 2 else tot // 2)]
    hi = v[np.searchsorted(c, tot // 2 + 1)] if tot % 2 == 0 else lo
    return 0.5 * (lo + hi)
def wmean(v, w):
    ok = ~np.isnan(v) & (w > 0); return float((v[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else np.nan

res = []
for pool in ("all", "test"):
    X = {a: g.set_index(["y", "compound_id", "event_v3"]).sort_index() for a, g in EV[EV.pool == pool].groupby("arm")}
    k0 = X["orig"].index
    assert all(X[a].index.equals(k0) for a in X), "event sets differ"
    drugs = np.array(sorted(k0.get_level_values("compound_id").unique())); nd = len(drugs)
    rng = np.random.default_rng(SEED); W = np.ones((B + 1, nd)); W[1:] = [np.bincount(rng.integers(0, nd, nd), minlength=nd) for _ in range(B)]
    for scope in ["pooled"] + list(YEARS):
        est = {}
        for a in X:
            g = X[a].reset_index(); g = g if scope == "pooled" else g[g.y == scope]
            med = g.groupby("compound_id")[list(AN)].median().reindex(drugs).values
            e2 = (g.groupby("compound_id")[list(AN)].min().reindex(drugs).values <= 10).astype(float); e2[np.isnan(med)] = np.nan
            nev = int(len(g)); n_dr = int((~np.isnan(med[:, 0])).sum()) if a == "orig" else None
            if a == "orig": N_DR, N_EV = n_dr, nev
            for j, an in enumerate(AN): est[(a, an)] = np.array([[wmed(med[:, j], w), wmean(e2[:, j], w)] for w in W])
        for an in AN:
            adv = {}
            for a in ("orig", "gafnow", "gafrel"):
                d = est[("aaw", an)] - est[(a, an)]; d[:, 1] *= -100.0; adv[a] = d
            adv["change"] = adv["gafrel"] - adv["gafnow"]
            for q, d in adv.items():
                for k, ep in enumerate(("E1", "E2")):
                    res.append(dict(pool=pool, scope=str(scope), analysis=an, endpoint=ep, quantity=q, estimate=d[0, k],
                                    ci_lo=np.nanpercentile(d[1:, k], 2.5), ci_hi=np.nanpercentile(d[1:, k], 97.5), n_drugs=N_DR, n_events=N_EV))
R = pd.DataFrame(res)
R.to_csv(HERE / "SI_go_release_advantage_full.tsv", sep="\t", index=False, float_format="%.6g")

def fmt(x, d):
    s = f"{x:+.{d}f}"; s = "+" + s[1:] if float(s) == 0 else s
    return s.replace("-", "\u2212")
def cell(r, d): return f"{fmt(r.estimate, d)} [{fmt(r.ci_lo, d)}, {fmt(r.ci_hi, d)}]"
P = R[R.scope == "pooled"].set_index(["analysis", "pool", "endpoint", "quantity"])
rows, nums = [], {}
for an in ("pct_p1", "pct_s3"):
    for pool in ("test", "all"):
        for ep, lab, d in (("E1", "E1 (percentile points)", 2), ("E2", "E2 (percentage points)", 1)):
            g = {q: P.loc[(an, pool, ep, q)] for q in ("gafnow", "gafrel", "change")}
            rows.append([AN[an], POOL[pool], lab, f"{int(g['gafnow'].n_drugs):,}", f"{int(g['gafnow'].n_events):,}",
                         cell(g["gafnow"], d), cell(g["gafrel"], d), cell(g["change"], d)])
            nums.setdefault(an, {}).setdefault(pool, {})[ep] = {q: [round(float(g[q][k]), 6) for k in ("estimate", "ci_lo", "ci_hi")] for q in g} | \
                {"n_drugs": int(g["gafnow"].n_drugs), "n_events": int(g["gafnow"].n_events)}
cur, dat = NMO.loc["gafnow_acc"], NMO.loc["gafrel_acc"]
assert int(dat.n_scored) == 0
rows.append(["Eculizumab case: PNH \u2192 NMO, 2008 slice", "All diseases", "Percentile of NMO among matched candidate diseases, mean of five seeds",
             "\u2014", "\u2014", f"{cur.mean_pct_cm:.2f} (range {cur.min_pct_cm:.2f}\u2013{cur.max_pct_cm:.2f})", "not assessed (NMO absent from the graph)", "not assessed"])
COLS = ["Comparison set", "Diseases", "Endpoint", "Drugs", "New trials", "Current GO release", "GO release dated to the slice", "Change (dated \u2212 current)"]
T = pd.DataFrame(rows, columns=COLS); T.to_csv(HERE / "SI_go_release_table.tsv", sep="\t", index=False)
with open(HERE / "SI_go_release_table.md", "w", encoding="utf-8") as fh:
    fh.write("| " + " | ".join(COLS) + " |\n|:--|:--|:--|--:|--:|--:|--:|--:|\n")
    for r in rows: fh.write("| " + " | ".join(r) + " |\n")
# change cells whose interval excludes zero, any pool / control / endpoint, pooled years
ch = R[(R.quantity == "change") & (R.scope == "pooled") & (R.analysis != "pct_s4")]
e1_orig = R[(R.quantity.isin(["gafnow", "gafrel"])) & (R.endpoint == "E1") & (R.scope == "pooled") & (R.analysis != "pct_s4")]
N = dict(table=nums, B=B,
         n_change_ci_excl0=int(((ch.ci_lo > 0) | (ch.ci_hi < 0)).sum()),
         e1_advantage_ci_excl0_all=bool((e1_orig.ci_lo > 0).all()),
         common_nodes={str(y): int(C.loc[y, "common_nodes"]) for y in YEARS},
         slice_disease_genes_2008=int(A.loc[("gafrel_acc", "S2008"), "slice_disease_genes"]),
         annotated_2008={"dated": int(A.loc[("gafrel_acc", "S2008"), "slice_disease_genes_annotated"]),
                         "current": int(A.loc[("gafnow_acc", "S2008"), "slice_disease_genes_annotated"])},
         nmo_current={"mean": float(cur.mean_pct_cm), "min": float(cur.min_pct_cm), "max": float(cur.max_pct_cm), "n_seeds": int(cur.n_seeds)},
         nmo_dated_scored=int(dat.n_scored), n_seeds=int(cur.n_seeds))
write_numbers(__file__, N, stem="SI_go_release_table")
print(T.to_string())
