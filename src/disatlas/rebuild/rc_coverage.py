"""Rebuild the rare/common neighbour-coverage tables of Supplementary Figure SI_rare_common
(data/derived/si_rare_common/: per_disease_first_other_class.tsv.gz, coverage_by_k_summary.tsv,
enrichment_other_class_frac10.tsv, rc_coverage_meta.json). Port of the authors'
revision/v32_rewrite/headline_swap/code/rc_enrichment_B.py (RELEASE=v2.1.0).

Neighbours: the 50 nearest diseases by cosine on the five-seed embedding (seeds 4707-4711, each L2-normalised,
concatenated / sqrt(5); equals the mean of the five per-seed cosines). Labels: is_rare. Nulls (1,000 label permutations
each, one rng stream default_rng(4707) consumed in the order uniform, degree, annotation): uniform; within training-
degree deciles; within annotation bins (5 quantile bins of genes x 5 of HPO terms among annotated diseases, 0 = none).
Inputs: data/zenodo/model/embedding_s47xx.npz, per_disease_first_other_class.tsv.gz (disease order, is_rare,
train_degree of the authors' tune-split training graph), rc_annotation_counts.tsv.gz (staged gene / HPO counts).

    python -m disatlas.rebuild.rc_coverage OUTDIR        (about 10 min on 2 CPUs)
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
D = ROOT / "data/derived/si_rare_common"
SEEDS = range(4707, 4712); KS = list(range(1, 51)); KMAX = 50; NPERM = 1000; SEED = 4707




def neighbours(ids):
    """50 nearest neighbours (descending five-seed mean cosine, float64) of every disease in `ids`."""
    from disatlas.model import load_ensemble
    M = load_ensemble()
    X = np.hstack([M.emb[s][[M.idx[s][i] for i in ids]].astype(np.float64) for s in M.seeds]) / np.sqrt(len(M.seeds))
    n = len(ids); NB = np.zeros((n, KMAX), np.int32)
    for s in range(0, n, 2000):
        Sm = X[s:s + 2000] @ X.T
        Sm[np.arange(Sm.shape[0]), np.arange(s, s + Sm.shape[0])] = -np.inf
        top = np.argpartition(-Sm, KMAX, axis=1)[:, :KMAX]
        o = np.argsort(-np.take_along_axis(Sm, top, 1), axis=1); NB[s:s + 2000] = np.take_along_axis(top, o, 1)
    return NB, X


def qbin(v, q=5):
    b = np.zeros(len(v), int); nz = v > 0
    if nz.any(): b[nz] = 1 + pd.qcut(v[nz], q, labels=False, duplicates="drop")
    return b


def build(nperm=NPERM):
    P = pd.read_csv(D / "per_disease_first_other_class.tsv.gz", sep="\t")
    A = pd.read_csv(D / "rc_annotation_counts.tsv.gz", sep="\t"); assert (A.disease_id.values == P.disease_id.values).all()
    ids = P.disease_id.astype(str).tolist(); rare = P.is_rare.to_numpy(int); deg = P.train_degree.to_numpy(); n = len(ids)
    NB, _ = neighbours(ids); rng = np.random.default_rng(SEED)
    deg_bin = pd.qcut(deg, 10, labels=False, duplicates="drop"); ann_bin = qbin(A.n_genes.to_numpy()) * 10 + qbin(A.n_hpo.to_numpy())

    def first_other(lab):
        oth = lab[NB] != lab[:, None]; return np.where(oth.any(1), oth.argmax(1) + 1, KMAX + 1)

    def perm_within(groups):
        lab = rare.copy()
        for g in np.unique(groups):
            w = np.where(groups == g)[0]; lab[w] = rare[rng.permutation(w)]
        return lab

    fo = first_other(rare); oth10 = (rare[NB[:, :10]] != rare[:, None]).sum(1)
    PD = pd.DataFrame(dict(disease_id=ids, is_rare=rare, first_other_class_rank=fo, n_other_class_in_10=oth10, train_degree=deg))
    rows = []
    for nm, grp in (("observed", None), ("null_uniform", np.zeros(n, int)), ("null_degree", deg_bin), ("null_annotation", ann_bin)):
        for r in range(1 if grp is None else nperm):
            lab = rare if grp is None else perm_within(grp); f = first_other(lab)
            for a, v in (("rare", 1), ("common", 0)):
                fa = f[lab == v]; o10 = (lab[NB[:, :10]] != lab[:, None]).sum(1)[lab == v] / 10
                rows.append(dict(model=nm, rep=r, anchor_class=a, frac_other10=float(o10.mean()), **{f"k{k}": float((fa <= k).mean()) for k in KS}))
    Dr = pd.DataFrame(rows); S = []
    for (m, a), g in Dr.groupby(["model", "anchor_class"]):
        for k in KS:
            v = g[f"k{k}"].values; S.append(dict(model=m, anchor_class=a, k=k, mean=v.mean(), lo=np.percentile(v, 2.5), hi=np.percentile(v, 97.5)))
    S = pd.DataFrame(S); obs = S[S.model == "observed"].set_index(["anchor_class", "k"])["mean"]
    for m in ("null_uniform", "null_degree", "null_annotation"):
        for a in ("rare", "common"):
            for k in KS:
                v = Dr[(Dr.model == m) & (Dr.anchor_class == a)][f"k{k}"].values
                S.loc[(S.model == m) & (S.anchor_class == a) & (S.k == k), "p_lower"] = (1 + (v <= obs[(a, k)]).sum()) / (nperm + 1)
    E = []
    for (m, a), g in Dr.groupby(["model", "anchor_class"]):
        v = g.frac_other10.values; o = Dr[(Dr.model == "observed") & (Dr.anchor_class == a)].frac_other10.iloc[0]
        E.append(dict(model=m, anchor_class=a, mean_frac_other_in_10=v.mean(), lo=np.percentile(v, 2.5), hi=np.percentile(v, 97.5), observed=o, p_lower=(1 + (v <= o).sum()) / (len(v) + 1)))
    meta = dict(n=n, n_rare=int(rare.sum()), n_common=int(n - rare.sum()), nperm=nperm,
                median_first_common_rank_for_rare=float(np.median(fo[rare == 1])), median_first_rare_rank_for_common=float(np.median(fo[rare == 0])),
                n_rare_no_common_in_50=int((fo[rare == 1] > 50).sum()))
    return PD, S, pd.DataFrame(E), meta


if __name__ == "__main__":
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True); PD, S, E, meta = build()
    PD.to_csv(out / "per_disease_first_other_class.tsv.gz", sep="\t", index=False); S.to_csv(out / "coverage_by_k_summary.tsv", sep="\t", index=False)
    E.to_csv(out / "enrichment_other_class_frac10.tsv", sep="\t", index=False); json.dump(meta, open(out / "rc_coverage_meta.json", "w"), indent=1)
