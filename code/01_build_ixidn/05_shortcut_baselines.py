"""Step 5 — matched shortcut baselines on IxIDN-Matched.

Manuscript: Results (matched shortcut baselines built from ontology or annotation similarity stay
close to chance) and the appendix section on IxIDN-Matched shortcut baselines; rarity-stratum
"vs. matched shortcuts" contrast in the appendix rare/common table.

Four non-degree shortcut scores, each a single pair feature used directly as a predictor:
  label_mondo_simgic   MONDO SimGIC between the two labels. Seven pairs (all involving one disease
                       without a SimGIC value) have no score and are dropped for this AUROC (n = 4,357).
  label_name_tfidf     TF-IDF cosine of the two disease names
  gene_count_logsum    log1p(n_genes_a) + log1p(n_genes_b), tier 1-2 causal genes
  gene_coverage_min    genecov, the gene-coverage matching key (0.500 by construction)

Shortcut stack: 5-fold out-of-fold logistic regression on the four scores (missing values filled
with the stratum median, then standardised; StratifiedKFold(5, shuffle, random_state=0)), fitted
separately within each stratum. The DisAtlas-minus-stack AUROC delta gets a 95% disease-cluster
bootstrap interval (1,000 resamples, seed 4707; a pair is weighted by the product of its two
diseases' resample counts).

The supplementary common-common benchmark (IxIDN-Matched plus 309 matched negatives) is evaluated
too when CCPLUS_ROWS (pair_uid -> row of the frozen matching features) is present; its DisAtlas
scores are read from tables/inputs/ixidn_matched_ccplus_scores.tsv.

Inputs: data/01_rebuilt/ixidn_matched.tsv; data/08S_eval/ixidn4364_scores.tsv;
data/00_raw/llm_annotations/disease_labels.parquet; frozen matching features and causal genes.
Outputs (data/01_rebuilt/):
  ixidn_matched_shortcut_baselines.tsv   score, stratum, n, n_pos, n_neg, auroc
  ixidn_matched_shortcut_stack.tsv       benchmark, scope, stratum, n, auroc, stack_auroc,
                                         delta, delta_lo, delta_hi, p_le0
"""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ixidn_lib as L  # noqa: E402

LABELS = os.path.join(L.RAW, "llm_annotations", "disease_labels.parquet")
SCORES = os.path.join(L.REPO_ROOT, "data", "08S_eval", "ixidn4364_scores.tsv")
CCPLUS_SCORES = os.path.join(L.REPO_ROOT, "tables", "inputs", "ixidn_matched_ccplus_scores.tsv")
CCPLUS_ROWS = os.path.join(L.FROZEN, "ccplus_pair_rows.tsv.gz")

SCORE_COLS = {"label_mondo_simgic": "mondo_simgic", "label_name_tfidf": "name_tfidf_cos",
              "gene_count_logsum": "lgp", "gene_coverage_min": "genecov"}
STACK = list(SCORE_COLS.values())
STRATA = ["all", "rare-rare", "rare-common", "common-common"]
MIN_CLASS = 20


def n_genes():
    cg = pd.read_parquet(L.CAUSAL_GENES)
    tier = pd.to_numeric(cg.tier.astype(str).str.extract(r"(\d)")[0], errors="coerce")
    return cg[tier <= 2].groupby("disease_id").gene_symbol.nunique()


def rarity(df):
    rare = pd.read_parquet(LABELS).set_index("disease_id").is_rare.astype(float)
    ra, rb = df.a.map(rare).values, df.b.map(rare).values
    assert not (np.isnan(ra).any() or np.isnan(rb).any()), "disease without a rare/common label"
    return np.where((ra == 1) & (rb == 1), "rare-rare", np.where((ra == 0) & (rb == 0), "common-common", "rare-common"))


def stack_oof(d):
    X = d[STACK].astype(float)
    X = X.fillna(X.median())
    X = (X - X.mean()) / (X.std() + 1e-9)
    cv = StratifiedKFold(5, shuffle=True, random_state=0)
    return cross_val_predict(LogisticRegression(max_iter=1000), X.values, d.y.values, cv=cv,
                             method="predict_proba")[:, 1]


def cluster_boot_delta(d, s1, s2, n_boot=1000, seed=4707):
    nodes = pd.unique(d[["a", "b"]].values.ravel())
    idx = {n: i for i, n in enumerate(nodes)}
    ai, bi, y = d.a.map(idx).values, d.b.map(idx).values, d.y.values
    x1, x2 = np.asarray(s1, float), np.asarray(s2, float)
    rng, out = np.random.default_rng(seed), []
    for _ in range(n_boot):
        w = np.bincount(rng.integers(0, len(nodes), len(nodes)), minlength=len(nodes))
        ww = w[ai] * w[bi]
        m = ww > 0
        if len(np.unique(y[m])) < 2:
            continue
        out.append(roc_auc_score(y[m], x1[m], sample_weight=ww[m]) - roc_auc_score(y[m], x2[m], sample_weight=ww[m]))
    out = np.asarray(out)
    lo, hi = np.percentile(out, [2.5, 97.5])
    return lo, hi, float((out <= 0).mean())


def single_scores(B):
    strata = [("all", np.ones(len(B), bool))]
    strata += [(f"pool={p}", (B.pool == p).values) for p in ("tune", "cross", "test")]
    strata += [(f"rarity={r}", (B.stratum == r).values) for r in STRATA[1:]]
    rows = []
    for name, col in SCORE_COLS.items():
        for sn, m in strata:
            d = B[m & B[col].notna().values]
            rows.append((name, sn, len(d), int(d.y.sum()), int((d.y == 0).sum()), roc_auc_score(d.y, d[col])))
    return pd.DataFrame(rows, columns=["score", "stratum", "n", "n_pos", "n_neg", "auroc"])


def stack_contrast(X, bench):
    rows = []
    for scope, D1 in [("All pairs", X), ("Test split pairs", X[X.pool == "test"])]:
        for st in STRATA:
            d = (D1 if st == "all" else D1[D1.stratum == st]).copy()
            if min(d.y.sum(), (d.y == 0).sum()) < MIN_CLASS:
                continue
            d["stk"] = stack_oof(d)
            a, sa = roc_auc_score(d.y, d.p), roc_auc_score(d.y, d.stk)
            lo, hi, ple0 = cluster_boot_delta(d, d.p, d.stk)
            rows.append((bench, scope, st, len(d), a, sa, a - sa, lo, hi, ple0))
    return rows


def ccplus_frame(ng):
    rows = pd.read_csv(CCPLUS_ROWS, sep="\t")
    S = pd.read_csv(CCPLUS_SCORES, sep="\t")
    assert (rows.pair_uid.values == S.pair_uid.values).all()
    F = pd.read_parquet(L.MATCH_FEATURES).set_index("row").loc[rows.row]
    X = S[["pair_uid", "label", "stratum", "source", "pool"]].copy()
    X["p"] = S.disatlas_ens5.values
    for c in ["a", "b", "mondo_simgic", "name_tfidf_cos", "genecov"]:
        X[c] = F[c].values
    X["lgp"] = np.log1p(X.a.map(ng).fillna(0)) + np.log1p(X.b.map(ng).fillna(0))
    X["y"] = X.label
    return X


def main():
    B = pd.read_csv(L.out_path("ixidn_matched.tsv"), sep="\t")
    s = pd.read_csv(SCORES, sep="\t")
    p = dict(zip(map(frozenset, zip(s.a, s.b)), s.p))
    B["p"] = [p[frozenset(x)] for x in zip(B.a, B.b)]
    B["y"] = B.label
    B["stratum"] = rarity(B)

    single = single_scores(B)
    single.to_csv(L.out_path("ixidn_matched_shortcut_baselines.tsv"), sep="\t", index=False)
    print(single[single.stratum == "all"].to_string(index=False))

    rows = stack_contrast(B, "IxIDN-Matched")
    if os.path.exists(CCPLUS_ROWS):
        X = ccplus_frame(n_genes())
        assert (X.source == "IxIDN-Matched").sum() == len(B)
        rows += stack_contrast(X, "IxIDN-Matched-CC+")
    else:
        print(f"{os.path.relpath(CCPLUS_ROWS, L.REPO_ROOT)} not found: supplementary benchmark skipped")
    out = pd.DataFrame(rows, columns=["benchmark", "scope", "stratum", "n", "auroc", "stack_auroc",
                                      "delta", "delta_lo", "delta_hi", "p_le0"])
    out.to_csv(L.out_path("ixidn_matched_shortcut_stack.tsv"), sep="\t", index=False)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
