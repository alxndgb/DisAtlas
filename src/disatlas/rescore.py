"""Score the benchmarks with an ensemble by cosine: IxIDN-Matched (4,364 pairs), IxIDN+/IxIDN- and the FDA
repurposing dev controls. The pair score is the cosine of the two disease embeddings, per seed (c_s<seed>) and the
5-seed mean (c). Writes per-pair scores and benchmark_auroc.json to --out.

    python -m disatlas.rescore --out build/scores                 # model from disatlas.model.final_model_dir()
    python -m disatlas.rescore --model DIR --out DIR [--all-pairs]
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

from disatlas import paths

BM = paths.RAW / "benchmarks"


def auroc(y, s):
    from scipy.stats import rankdata
    y, s = np.asarray(y, bool), np.asarray(s, float)
    r = rankdata(s)
    n1, n0 = y.sum(), (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def benchmarks():
    m = pd.read_csv(BM / "ixidn_matched" / "ixidn_matched_pairs.tsv", sep="\t")
    pos = pd.read_parquet(BM / "ixidn" / "ixidn_plus_pairs.parquet").assign(y=1)
    neg = pd.read_parquet(BM / "ixidn" / "ixidn_minus_pairs.parquet").assign(y=0)
    f = pd.read_csv(BM / "fda_repurposing" / "fda_repurposing_controls.tsv", sep="\t")
    f = f[f.split.eq("dev") & f.in_model].assign(y=lambda d: d.fda_verified_pair.astype(int))
    return {"ixidn_matched": m, "ixidn_pm": pd.concat([pos, neg], ignore_index=True), "fda_dev": f}


def aurocs(name, d, p):
    ok = ~np.isnan(p)
    out = {name: auroc(d.y[ok], p[ok])}
    if name == "ixidn_matched":
        for col, val in (("split", "test"), ("pair_class", "rare-rare"), ("pair_class", "rare-common")):
            k = (d[col] == val).values & ok
            out[f"{name}|{val}"] = auroc(d.y[k], p[k])
    return out


def score(model, out):
    from disatlas.model import load_ensemble
    m = load_ensemble(str(model))
    os.makedirs(out, exist_ok=True)
    res = {}
    for name, d in benchmarks().items():
        P = m.score_pairs(d.disease_a, d.disease_b, per_seed=True)
        t = d[[c for c in ("pair_id", "disease_a", "disease_b", "y") if c in d]].copy()
        for s, row in zip(m.seeds, P):
            t[f"c_s{s}"] = row
            res.update({f"{k}|s{s}": v for k, v in aurocs(name, d, row).items()})
        t["c"] = P.mean(0)
        res.update(aurocs(name, d, t.c.values))
        if name == "ixidn_matched":
            t.drop(columns=["disease_a", "disease_b", "y"]).to_csv(
                os.path.join(out, "ixidn_matched_disatlas_cosine.tsv"), sep="\t", index=False)
        else:
            t.to_parquet(os.path.join(out, f"{name}_scores.parquet"), index=False)
    json.dump(res, open(os.path.join(out, "benchmark_auroc.json"), "w"), indent=1, sort_keys=True)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--all-pairs", action="store_true",
                    help="also write the all-pairs cosine matrix to <out>/cosine_matrix (disatlas.pair_matrix)")
    a = ap.parse_args()
    from disatlas.model import final_model_dir
    model = a.model or final_model_dir()
    res = score(model, a.out)
    print(json.dumps({k: round(v, 4) for k, v in res.items() if "|s" not in k}, indent=1))
    if a.all_pairs:
        from disatlas import pair_matrix
        out = os.path.join(a.out, "cosine_matrix")
        pair_matrix.build(out, str(model))
        if (pair_matrix.SHIPPED / "cosine_ens5_13021.npy").exists():
            print(json.dumps(pair_matrix.compare(out)))


if __name__ == "__main__":
    main()
