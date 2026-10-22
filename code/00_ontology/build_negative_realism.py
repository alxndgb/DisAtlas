"""Negative-set realism and hard-negative stress test (appendix: Negative-set realism, Figure SI_negrealism).

A. ORDO max-LCA distance for random ORDO disease pairs and for IxIDN-Matched positives / negatives.
B. DisAtlas AUROC of all positives against only the negatives that are as close / similar as a typical
   positive (organ system, ORDO clinical area, MONDO / ORDO distance, HPO phenotype, GO function).
C. Within-stratum AUROC where both classes share an organ system / ORDO clinical area, with HPO Resnik-BMA,
   GO Resnik gene-BMA and MONDO SimGIC on the same pairs (unscorable ranked last).
Each AUROC on all 4,364 pairs and on test-split pairs; 95% CIs from 1000 disease-cluster bootstrap
replicates. One generator (seed 4707) is shared by the random-pair draw and the bootstraps, in that order.

Inputs (repository-relative):
  data/08S_eval/ixidn4364_scores.tsv                       DisAtlas score p, label y, pool (stored)
  data/01_rebuilt/ixidn_matched.tsv                        shared_mondo_organsys, shared_ordo_area, mondo_simgic
  data/00_rebuilt_ontology/ixidn_matched_hpo_scores.csv    HPO Resnik-BMA (build_hpo_baseline.py)
  figures/main/data/ixidn_matched_go_resnik.tsv.gz         GO Resnik gene-BMA (stored)
  data/00_rebuilt_ontology/ixidn_matched_ordo_mondo_distance.tsv  (build_ontology_dependency.py)
  data/00_raw/ontologies/ordo/                             ORDO ancestors, random-pair disease list
Outputs (data/00_rebuilt_ontology/):
  negative_realism_ordo_distance.tsv.gz      set, ordo_dist          (figure panel a)
  negative_realism_hard_negative_auroc.tsv   figure panel b rows
  negative_realism_summary.tsv               medians and shared-class shares quoted in the appendix

    python code/00_ontology/build_negative_realism.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import max_lca  # noqa: E402
from ontology import REPO_ROOT  # noqa: E402

OUT = REPO_ROOT / "data" / "00_rebuilt_ontology"
SCORES = REPO_ROOT / "data" / "08S_eval" / "ixidn4364_scores.tsv"
BENCH = REPO_ROOT / "data" / "01_rebuilt" / "ixidn_matched.tsv"
HPO = OUT / "ixidn_matched_hpo_scores.csv"
GO = REPO_ROOT / "figures" / "main" / "data" / "ixidn_matched_go_resnik.tsv.gz"
ONTODIST = OUT / "ixidn_matched_ordo_mondo_distance.tsv"
SEED, B, N_RANDOM = 4707, 1000, 50000
ALL, TEST = "All 4,364 pairs", "Test-split pairs"


def pair_table():
    S = pd.read_csv(SCORES, sep="\t")[["blind_id", "a", "b", "y", "pool", "p"]]
    K = pd.read_csv(BENCH, sep="\t")[["a", "b", "shared_mondo_organsys", "shared_ordo_area", "mondo_simgic"]]
    key = lambda d: [tuple(sorted(x)) for x in zip(d.a, d.b)]   # benchmark pairs are unordered
    S["k"], K["k"] = key(S), key(K)
    Y = (S.merge(K.drop(columns=["a", "b"]), on="k").drop(columns="k")
         .merge(pd.read_csv(HPO)[["blind_id", "hpo_resnik_bma"]], on="blind_id")
         .merge(pd.read_csv(GO, sep="\t")[["blind_id", "go_resnik_genebma"]], on="blind_id")
         .merge(pd.read_csv(ONTODIST, sep="\t"), on="blind_id"))
    assert len(Y) == 4364 and Y.y.sum() == 2182
    for c in ("ordo_dist", "mondo_dist"):
        Y.loc[Y[c] < 0, c] = np.nan
    return Y


def neg_rule(T):
    P = T[T.y == 1]
    med = lambda c: P[c].median()
    close = {"mondo": T.mondo_dist <= med("mondo_dist"), "ordo": T.ordo_dist <= med("ordo_dist"),
             "hpo": T.hpo_resnik_bma >= med("hpo_resnik_bma"), "go": T.go_resnik_genebma >= med("go_resnik_genebma")}
    return {
        "All negatives (reference)": pd.Series(True, index=T.index),
        "Share an organ system with each other": T.shared_mondo_organsys.eq(1),
        "Share an ORDO clinical area": T.shared_ordo_area.eq(1),
        "MONDO distance <= positive median (ontology-close)": close["mondo"],
        "ORDO distance <= positive median (ontology-close)": close["ordo"],
        "HPO Resnik-BMA >= positive median (phenotype-similar)": close["hpo"],
        "GO Resnik gene-BMA >= positive median (function-similar)": close["go"],
        "Close or similar on >= 2 of the 4 axes above": sum(v.astype(int) for v in close.values()) >= 2,
    }


def main():
    rng = np.random.default_rng(SEED)
    Y = pair_table()

    # A. distributions (the random-pair draw comes first on the shared generator)
    rd = max_lca.random_ordo_pairs(rng, N_RANDOM)
    L = pd.concat([pd.DataFrame({"set": "random_disease_pairs", "ordo_dist": rd}),
                   pd.DataFrame({"set": "ixidn_matched_positives", "ordo_dist": Y.ordo_dist[Y.y == 1]}),
                   pd.DataFrame({"set": "ixidn_matched_negatives", "ordo_dist": Y.ordo_dist[Y.y == 0]})]).dropna()
    L.to_csv(OUT / "negative_realism_ordo_distance.tsv.gz", sep="\t", index=False,
             compression={"method": "gzip", "mtime": 0})          # byte-stable output
    summ = []
    for name, v, cls in [("Random ORDO disease pairs", rd, None), ("IxIDN-Matched positives", Y.ordo_dist[Y.y == 1], 1),
                         ("IxIDN-Matched negatives", Y.ordo_dist[Y.y == 0], 0)]:
        v = np.asarray(v, float); v = v[~np.isnan(v)]
        r = {"set": name, "n_scorable": len(v), "median": np.median(v), "p10": np.percentile(v, 10),
             "p90": np.percentile(v, 90)}
        if cls is not None:
            r["pct_share_ordo_area"] = 100 * Y.shared_ordo_area[Y.y == cls].eq(1).mean()
            r["pct_share_organ_system"] = 100 * Y.shared_mondo_organsys[Y.y == cls].eq(1).mean()
        summ.append(r)
    pd.DataFrame(summ).to_csv(OUT / "negative_realism_summary.tsv", sep="\t", index=False, float_format="%.4g")

    # B / C. AUROC with disease-cluster bootstrap
    def boot(df, col="p"):
        nodes = pd.unique(df[["a", "b"]].values.ravel()); idx = {n: i for i, n in enumerate(nodes)}
        ai = df.a.map(idx).values; bi = df.b.map(idx).values; y = df.y.values; p = df[col].values
        out = []
        for _ in range(B):
            w = np.bincount(rng.integers(0, len(nodes), len(nodes)), minlength=len(nodes)); ww = w[ai] * w[bi]
            k = ww > 0
            if len(np.unique(y[k])) < 2:
                continue
            out.append(roc_auc_score(y[k], p[k], sample_weight=ww[k]))
        return np.percentile(out, [2.5, 97.5])

    rows = []
    for scope, T in ((ALL, Y), (TEST, Y[Y.pool == "test"])):
        for name, keep in neg_rule(T).items():
            d = T[(T.y == 1) | ((T.y == 0) & keep.fillna(False))]
            n0 = int((d.y == 0).sum())
            if n0 < 20:
                continue
            lo, hi = boot(d)
            rows.append({"scope": scope, "subset": name, "n_pos": int(d.y.sum()), "n_neg": n0,
                         "pct_neg_kept": 100 * n0 / (T.y == 0).sum(), "auroc_disatlas": roc_auc_score(d.y, d.p),
                         "ci_lo": lo, "ci_hi": hi})
        for name, col in (("Both classes share an organ system", "shared_mondo_organsys"),
                          ("Both classes share an ORDO clinical area", "shared_ordo_area")):
            d = T[T[col].eq(1)]
            lo, hi = boot(d)
            r = {"scope": scope, "subset": name, "n_pos": int(d.y.sum()), "n_neg": int((d.y == 0).sum()),
                 "pct_neg_kept": 100 * (d.y == 0).sum() / (T.y == 0).sum(), "auroc_disatlas": roc_auc_score(d.y, d.p),
                 "ci_lo": lo, "ci_hi": hi}
            for b, lab in (("hpo_resnik_bma", "hpo_resnik_bma"), ("go_resnik_genebma", "go_resnik_gene_bma"),
                           ("mondo_simgic", "mondo_simgic")):   # unscorable ranked last
                r[f"auroc_{lab}"] = roc_auc_score(d.y, d[b].fillna(-1e9))
            rows.append(r)
    cols = ["scope", "subset", "n_pos", "n_neg", "pct_neg_kept", "auroc_disatlas", "ci_lo", "ci_hi",
            "auroc_hpo_resnik_bma", "auroc_go_resnik_gene_bma", "auroc_mondo_simgic"]
    R = pd.DataFrame(rows).reindex(columns=cols).round(4)
    R.to_csv(OUT / "negative_realism_hard_negative_auroc.tsv", sep="\t", index=False)
    print(R[["scope", "subset", "n_neg", "auroc_disatlas", "ci_lo", "ci_hi"]].to_string())


if __name__ == "__main__":
    main()
