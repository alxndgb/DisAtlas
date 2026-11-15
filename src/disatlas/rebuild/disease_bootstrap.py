"""Rebuild data/derived/ixidn_robustness/subgroups_shared_bootstrap_B2000.tsv (Supplementary Figure SI_llm_subgroups /
Table S9-S10 area: disease-cluster bootstrap AUROC of the node-disjoint and in-graph DisAtlas cosine readouts and three LLM
systems on the 4,357 pairs the node-disjoint model can score; scopes all / test / non_oncology / oncology) from released scores.

    python -m disatlas.rebuild.disease_bootstrap OUTDIR

Inputs: data/derived/fig6_comparison/per_run_scores.tsv (DisAtlas = five-seed mean cosine; LLM runs), node_disjoint_per_seed.tsv
(five-seed node-disjoint cosine; the 7 pairs it cannot score are dropped, not ranked last) and
data/derived/ixidn_robustness/ixidn4364_per_pair_scores.tsv.gz (either_is_cancer; LLM run means).
Bootstrap (authors' llm_subgroups_boot.py boot_auc): diseases indexed in order of appearance in the row-wise (a, b) sequence of
each scope, default_rng(4707), per draw a disease-cluster resample (pair weight w_a * w_b), draws without both classes skipped;
the same draws are shared by all systems (paired). dis_minus_sys = node-disjoint-referenced difference with the sign flipped
(system minus node-disjoint is negated), 95% intervals are 2.5/97.5 percentiles, p_two = 2 * min(P(d <= 0), P(d >= 0)).
Systems: nd = node-disjoint cosine (reference), p = in-graph DisAtlas cosine.
Test: tests/test_rebuild_disease_bootstrap.py.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[3]
SYS = ["nd", "p", "Claude Sci", "Claude Chat", "Gemini Pro 3.1"]


def boot_auc(D, cols, ref, B=2000, seed=4707):
    y = D.y.values; nodes = pd.unique(D[["a", "b"]].values.ravel()); ix = {n: i for i, n in enumerate(nodes)}
    ai, bi = D.a.map(ix).values, D.b.map(ix).values; rg = np.random.default_rng(seed); bs = {c: [] for c in cols}
    for _ in range(B):
        w = np.bincount(rg.integers(0, len(nodes), len(nodes)), minlength=len(nodes)); ww = (w[ai] * w[bi]).astype(float); q = ww > 0
        if len(set(y[q])) < 2: continue
        for c in cols: bs[c].append(roc_auc_score(y[q], D[c].values[q], sample_weight=ww[q]))
    out = {}
    for c in cols:
        b = np.array(bs[c]); r = dict(nb=len(b), auc=roc_auc_score(y, D[c].values), lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5))
        if c != ref:
            d = b - np.array(bs[ref]); r.update(d=r["auc"] - roc_auc_score(y, D[ref].values), dlo=np.percentile(d, 2.5), dhi=np.percentile(d, 97.5),
                                                p_two=float(min(1, 2 * min((d <= 0).mean(), (d >= 0).mean()))))
        out[c] = r
    return out


def load():
    DV = ROOT / "data/derived"
    S = pd.read_csv(DV / "fig6_comparison/per_run_scores.tsv", sep="\t"); assert len(S) == 4364
    N = pd.read_csv(DV / "fig6_comparison/node_disjoint_per_seed.tsv", sep="\t")
    P = pd.read_csv(DV / "ixidn_robustness/ixidn4364_per_pair_scores.tsv.gz", sep="\t")
    assert (P.blind_id.values == S.blind_id.values).all()
    S["nd"] = S.blind_id.map(dict(zip(N.blind_id, N.iloc[:, 2:].mean(axis=1))))
    S = S.rename(columns={"DisAtlas": "p"})
    for k in ("Claude Sci", "Claude Chat", "Gemini Pro 3.1"): S[k] = P[k].values
    S["either_is_cancer"] = P.either_is_cancer.values
    return S


def compute(B=2000, scopes=None):
    S = load(); SN = S[S.nd.notna()].reset_index(drop=True); assert len(SN) == 4357
    onc = SN.either_is_cancer.values == 1
    SC = {"all": np.ones(len(SN), bool), "test": (SN.pool == "test").values, "non_oncology": ~onc, "oncology": onc}
    rows = []
    for sc, m in SC.items():
        if scopes and sc not in scopes: continue
        r = boot_auc(SN[m].reset_index(drop=True), SYS, "nd", B=B)
        for c in SYS:
            rows.append(dict(scope=sc, system=c, n_pairs=int(m.sum()), n_pos=int(SN.y[m].sum()), B=B, n_valid=r[c]["nb"], auroc=r[c]["auc"], lo=r[c]["lo"], hi=r[c]["hi"],
                             dis_minus_sys=-r[c].get("d", np.nan), dms_lo=-r[c].get("dhi", np.nan), dms_hi=-r[c].get("dlo", np.nan), p_two=r[c].get("p_two", np.nan)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    o = Path(sys.argv[1]); o.mkdir(parents=True, exist_ok=True)
    sc = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    T = compute(scopes=sc); T.to_csv(o / ("subgroups_shared_bootstrap_B2000" + ("_" + "_".join(sc) if sc else "") + ".tsv"), sep="\t", index=False); print("wrote", o)
