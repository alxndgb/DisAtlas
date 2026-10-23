"""Matched shortcut baselines and the rarity-stratum contrast DisAtlas minus matched shortcuts.

Manuscript: Results (matched shortcut baselines stay close to chance; rarity stratification, Figure 5c),
appendix section "Shortcut audit" and the rare/common summary table ("vs. matched shortcuts" column).

Matched shortcut baselines are the four non-degree pair features, each used directly as a predictor:
  label_mondo_simgic   MONDO SimGIC between the two labels; the 7 pairs without a MONDO term are dropped
                       from its single-score AUROC (4,357 pairs)
  label_name_tfidf     TF-IDF cosine of the two disease names
  gene_count_logsum    log1p(n_genes_a) + log1p(n_genes_b), tier 1-2 causal genes (column lgp)
  gene_coverage_min    the gene-coverage matching key (AUROC 0.500 by construction)
Combined matched shortcuts: within each stratum, out-of-fold probabilities of a logistic model on the four
features (missing values filled with the stratum median, z-scored; 5-fold stratified CV, shuffled,
random_state 0). DisAtlas minus combined: point difference of AUROCs; 95% CI from a paired disease-cluster
bootstrap (1,000 draws, seed 4707; each draw resamples diseases and weights a pair by the product of its two
diseases' counts).

Inputs (repo-relative):
  data/01_rebuilt/ixidn_matched.tsv          IxIDN-Matched pairs with matching features
  data/08S_eval/ixidn4364_scores.tsv         DisAtlas (full) ensemble probability p and rarity stratum
Outputs (next to this script):
  matched_shortcut_baselines.tsv             single-score AUROC per baseline and stratum
  rarity_vs_matched_shortcuts.tsv            DisAtlas, combined matched shortcuts and the paired contrast
    python figures/main/rarity_vs_matched_shortcuts.py
"""
import os
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

FEATURES = {"label_mondo_simgic": "mondo_simgic", "label_name_tfidf": "name_tfidf_cos",
            "gene_count_logsum": "lgp", "gene_coverage_min": "genecov"}
STRATA = ["rare-rare", "rare-common", "common-common", "all"]

B = pd.read_csv(os.path.join(REPO_ROOT, "data/01_rebuilt/ixidn_matched.tsv"), sep="\t",
                usecols=["a", "b", "label"] + list(FEATURES.values()))
S = pd.read_csv(os.path.join(REPO_ROOT, "data/08S_eval/ixidn4364_scores.tsv"), sep="\t",
                usecols=["a", "b", "y", "p", "rare_common_pair"])
key = lambda d: [frozenset(x) for x in zip(d.a, d.b)]
B["k"] = key(B); S["k"] = key(S)
X = B.merge(S.drop(columns=["a", "b"]), on="k", how="left", validate="one_to_one")
assert len(X) == 4364 and X.p.notna().all() and (X.label == X.y).all()


def stack_oof(d):
    Z = d[list(FEATURES.values())].astype(float); Z = Z.fillna(Z.median()); Z = (Z - Z.mean()) / (Z.std() + 1e-9)
    return cross_val_predict(LogisticRegression(max_iter=1000), Z.values, d.y.values,
                             cv=StratifiedKFold(5, shuffle=True, random_state=0), method="predict_proba")[:, 1]


def boot_delta(d, s1, s2, B_=1000, seed=4707):
    nodes = pd.unique(d[["a", "b"]].values.ravel()); idx = {n: i for i, n in enumerate(nodes)}
    ai = d.a.map(idx).values; bi = d.b.map(idx).values; y = d.y.values; r = np.random.default_rng(seed); out = []
    for _ in range(B_):
        w = np.bincount(r.integers(0, len(nodes), len(nodes)), minlength=len(nodes)); ww = w[ai] * w[bi]; m = ww > 0
        if len(np.unique(y[m])) < 2: continue
        out.append(roc_auc_score(y[m], s1[m], sample_weight=ww[m]) - roc_auc_score(y[m], s2[m], sample_weight=ww[m]))
    return np.percentile(out, [2.5, 97.5])


single, rows = [], []
for st in STRATA:
    d = (X if st == "all" else X[X.rare_common_pair == st]).reset_index(drop=True)
    for nm_, c_ in FEATURES.items():
        ok = d[c_].notna()
        single.append(dict(score=nm_, stratum=st, n=int(ok.sum()), n_pos=int(d.y[ok].sum()), n_neg=int((1 - d.y[ok]).sum()),
                           auroc=roc_auc_score(d.y[ok], d[c_][ok])))
    stk = stack_oof(d); a_ = roc_auc_score(d.y, d.p); s_ = roc_auc_score(d.y, stk); lo, hi = boot_delta(d, d.p.values, stk)
    rows.append(dict(stratum=st, n_pos=int(d.y.sum()), n_neg=int((1 - d.y).sum()), disatlas_auroc=a_,
                     combined_matched_shortcuts_auroc=s_, delta=a_ - s_, delta_ci_lo=lo, delta_ci_hi=hi,
                     significant=bool(lo > 0 or hi < 0)))
SB, RT = pd.DataFrame(single), pd.DataFrame(rows)
SB.to_csv(os.path.join(HERE, "matched_shortcut_baselines.tsv"), sep="\t", index=False, float_format="%.6f")
RT.to_csv(os.path.join(HERE, "rarity_vs_matched_shortcuts.tsv"), sep="\t", index=False, float_format="%.6f")
al = SB[SB.stratum == "all"].set_index("score").auroc
assert round(al["label_mondo_simgic"], 4) == 0.5418 and al.min() == 0.5
print("matched shortcut AUROC range (all pairs): %.3f-%.3f" % (al.min(), al.max()))
print(RT.round(3).to_string(index=False))
