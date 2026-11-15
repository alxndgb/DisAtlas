"""Rebuild data/derived/si_heur_vs_ai/heuristic_minus_ai_B2000.tsv (Supplementary Table S7, SI_heuristic_vs_ai):
training-free graph scores minus each AI system, paired AUROC differences on the 4,357 node-disjoint IxIDN-Matched
pairs and the 1,093-pair test split. Port of the authors' revision/response/_ai_comparison/heur_vs_ai.py.

Scores. Heuristics (held-out diseases placed from 55 feature neighbours, or keeping their own edges): staged
data/derived/si_baselines/baseline_pair_scores.tsv.gz columns ho_tfidf_svd_place55, ho_aa_place55, ho_cn_place55,
ho_tfidf_graph_svd_cosine. AI systems: mean over the available runs of the released comparator scores (Claude Chat, Claude Sci, Gemini: the stored
per-pair means data/derived/ixidn_robustness/ixidn4364_per_pair_scores.tsv.gz; Phylo, TxGNN: recomputed from
data/derived/fig6_comparison/per_run_scores.tsv, one column per system|run); a pair a system did not score (all runs
missing) is ranked last (minimum - 1). Pairs: the 4,357 pairs with a node-disjoint DisAtlas score (node_disjoint_per_seed.tsv),
in the row order of per_run_scores.tsv (the order fixes which Poisson draw belongs to which disease); test = pool 'test' (1,093).
Bootstrap: per subset one rng default_rng(0).poisson(1.0, (B, n_diseases)) with diseases in order of first appearance
over (a, b); pair weight = w_a * w_b; weighted roc_auc_score on pairs with weight > 0; the same B = 2,000 draws for every
score and system; 95% percentile intervals; differences paired on the draws.

    python -m disatlas.rebuild.si_heur_vs_ai OUT.tsv      (several minutes on 1 CPU)
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[3]
D = ROOT / "data/derived"
B = 2000
H = {"TF-IDF graph SVD, placed": "ho_tfidf_svd_place55", "Adamic-Adar, placed": "ho_aa_place55",
     "Common neighbours, placed": "ho_cn_place55", "TF-IDF graph SVD, own edges": "ho_tfidf_graph_svd_cosine"}
A = {"Claude Chat": "Claude-Chat Opus 5.5", "Claude Sci": "Claude-Sci Opus 5.5", "Gemini Pro 3.1": "Gemini Pro 3.1",
     "Phylo Biomni": "Phylo", "TxGNN": "TxGNN published ckpt"}


def pair_table():
    P = pd.read_csv(D / "fig6_comparison/per_run_scores.tsv", sep="\t")
    N = pd.read_csv(D / "fig6_comparison/node_disjoint_per_seed.tsv", sep="\t")
    seeds = [c for c in N.columns if "|s" in c]
    P = P.merge(pd.DataFrame({"blind_id": N.blind_id, "nd": N[seeds].mean(axis=1).where(N[seeds].notna().all(axis=1))}),
                on="blind_id", how="left", validate="1:1")
    S = pd.read_csv(D / "si_baselines/baseline_pair_scores.tsv.gz", sep="\t").set_index("blind_id")
    df = P[P.nd.notna()].reset_index(drop=True); assert len(df) == 4357
    for k, c in H.items(): df[k] = S.loc[df.blind_id, c].to_numpy(float)
    L = pd.read_csv(D / "ixidn_robustness/ixidn4364_per_pair_scores.tsv.gz", sep="\t").set_index("blind_id")
    for k, g in A.items():
        if k in L.columns:      # LLM run means as stored by the authors (the 1-ulp summation differences of a recomputed mean move ties)
            df[k] = L.loc[df.blind_id, k].to_numpy(float)
        else:
            df[k] = df[[c for c in P.columns if c.split("|")[0] == g]].mean(axis=1, skipna=True)
    for k in list(H) + list(A):
        s = df[k].astype(float); df[k] = s.fillna(s.min() - 1)
    return df


def build(nboot=B):
    df = pair_table(); rows = []
    for sub, m in (("all_4357", np.ones(len(df), bool)), ("test_1093", (df.pool == "test").to_numpy())):
        d = df[m].reset_index(drop=True); y = d.y.to_numpy()
        assert sub != "test_1093" or len(d) == 1093, len(d)
        dis = pd.Index(pd.unique(pd.concat([d.a, d.b]))); ia, ib = dis.get_indexer(d.a), dis.get_indexer(d.b)
        W = np.random.default_rng(0).poisson(1.0, size=(nboot, len(dis))).astype(float); PW = W[:, ia] * W[:, ib]
        cols = list(H) + list(A)
        au = {k: roc_auc_score(y, d[k]) for k in cols}
        bs = {k: np.array([roc_auc_score(y[w > 0], d[k].to_numpy()[w > 0], sample_weight=w[w > 0]) for w in PW]) for k in cols}
        for h in H:
            for a in A:
                g = bs[h] - bs[a]
                rows.append(dict(subset=sub, heuristic=h, system=a, heuristic_auroc=au[h], system_auroc=au[a], diff=au[h] - au[a],
                                 ci_lo=np.percentile(g, 2.5), ci_hi=np.percentile(g, 97.5)))
        for k in cols:
            rows.append(dict(subset=sub, heuristic=k, system="(self)", heuristic_auroc=au[k], ci_lo=np.percentile(bs[k], 2.5), ci_hi=np.percentile(bs[k], 97.5)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    build().to_csv(sys.argv[1], sep="\t", index=False, float_format="%.6f")
