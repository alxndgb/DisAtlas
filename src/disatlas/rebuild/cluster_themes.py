"""Rebuild the per-cluster annotation tables of the SI cluster-themes table (Supplementary Table S10, K = 55) from the released
cluster labels and public ontologies.

    python -m disatlas.rebuild.cluster_themes OUTDIR

Inputs (all released or public):
  data/zenodo/model_outputs/clusters/disease_clusters_k55.parquet   K = 55 labels of the 13,021 diseases
  data/raw/mondo/{mondo_body_system_13021.tsv, body_system_classes.tsv, mondo_terms_13021.tsv.gz, mondo_terms.tsv}   MONDO 2026-02-06
  data/raw/go/go_bp_terms.tsv.gz, data/raw/specific_channels/gene2go.parquet                                        GO 2026-07-26
  data/raw/hpo/disease_phenotype.parquet and the released disease-gene edges (disatlas.graph.raw_build.gene_edges)
Outputs (names as in data/derived/si_cluster_themes): bodysystem_by_cluster_long, mondo_terms_top5_by_cluster,
mechanism_terms_by_cluster, go_bp_top5_by_cluster, cluster_annotation_k55 (parquet).

The enrichment code is the authors' (anatomy_audit.py): hypergeometric upper tail over the annotated diseases, BH over
clusters x terms, background 10..2000 diseases (MONDO) or 10..500 genes (GO BP), pair statistics from 300,000 random pairs
(default_rng(20261102)) and within-cluster pairs (all pairs, or 400,000 sampled without replacement from the same stream).
Not rebuilt: pairs_module_jaccard_fold and graph_within_frac_specific (module support and the specific-channel graph rest on
licence-restricted evidence; the shipped values are kept in cluster_annotation_k55.parquet).
"""
import json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[3]
MECH_ROOTS = {"MONDO:7770011", "MONDO:0019052", "MONDO:0003778"}
ROOTS = {"MONDO:0000001", "MONDO:0700096", "MONDO:7770006"}
ORGAN_COLS = ["cluster", "category", "n_annot", "k", "share", "fold", "fdr"]


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = p[o] * n / (np.arange(n) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q); out[o] = np.minimum(q, 1); return out


def enrich(Xb, avail, lab, K, min_k=1):
    Xb = Xb.tocsr().astype(np.int32)
    Ntot = int(avail.sum()); bg = np.asarray(Xb[avail].sum(0)).ravel(); res = {}
    for c in range(K):
        m = (lab == c) & avail; n = int(m.sum())
        k = np.asarray(Xb[m].sum(0)).ravel() if n else np.zeros(Xb.shape[1])
        p = hypergeom.sf(k - 1, Ntot, bg, n) if n else np.ones(Xb.shape[1])
        p[k < min_k] = 1.0; res[c] = (n, k, p)
    return res, bg, Ntot


def load_inputs():
    MO = ROOT / "data/zenodo/model_outputs"; RAW = ROOT / "data/raw"
    ids = json.load(open(MO / "cosine_matrix/ids_13021.json"))
    cl = pd.read_parquet(MO / "clusters/disease_clusters_k55.parquet").set_index("disease_id").loc[ids]
    lab = cl.cluster.to_numpy(); K = int(lab.max()) + 1; assert K == 55 and len(ids) == 13021
    M = pd.read_csv(RAW / "mondo/mondo_body_system_13021.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert list(M.disease_id) == ids
    top = pd.read_csv(RAW / "mondo/body_system_classes.tsv", sep="\t")
    mt = pd.read_csv(RAW / "mondo/mondo_terms.tsv", sep="\t", keep_default_na=False)
    mdt = pd.read_csv(RAW / "mondo/mondo_terms_13021.tsv.gz", sep="\t")
    return ids, lab, K, M, top, mt, mdt


def mondo_tables(ids, lab, K, M, top, mt, mdt):
    N = len(ids); idx = {d: i for i, d in enumerate(ids)}
    name = dict(zip(mt.term_id, mt.name)); topids = list(top.mondo_id); topname = dict(zip(top.mondo_id, top.name))
    # body-system classes
    Y = np.zeros((N, len(topids)), bool)
    for i, cat in enumerate(M.categories):
        for t in cat.split("|"):
            if t: Y[i, topids.index(t)] = True
    has = Y.any(1)
    res, bg, Ntot = enrich(sp.csr_matrix(Y), has, lab, K)
    allq = bh(np.concatenate([res[c][2] for c in range(K)])).reshape(K, len(topids)); rows = []
    for c in range(K):
        n, k, p = res[c]; share = k / max(n, 1); fold = share / (bg / Ntot)
        for j, t in enumerate(topids):
            rows.append(dict(cluster=f"c{c}", category=topname[t], n_annot=n, k=int(k[j]), share=share[j], fold=fold[j], fdr=allq[c, j]))
    body = pd.DataFrame(rows)
    # per-disease term sets (mapped diseases only)
    mapped = np.array([bool(m) for m in M.mondo_id])
    dset = defaultdict(set)
    for d, t in zip(mdt.disease_id, mdt.term_id): dset[d].add(t)
    # grouping terms below the 19 classes
    drop = set(topids) | ROOTS
    T = sorted(mt.term_id[~mt.term_id.isin(drop)] if False else [t for t in set(mdt.term_id) if t not in drop])
    cnt = mdt[~mdt.term_id.isin(drop)].term_id.value_counts()
    T = sorted(t for t in cnt.index if 10 <= cnt[t] <= 2000); ti = {t: i for i, t in enumerate(T)}
    sub = mdt[mdt.term_id.isin(ti)]
    XM = sp.csr_matrix((np.ones(len(sub), np.int32), (sub.disease_id.map(idx).to_numpy(), sub.term_id.map(ti).to_numpy())), shape=(N, len(T)))
    resM, bgM, NM = enrich(XM, mapped, lab, K, min_k=3)
    allq = bh(np.concatenate([resM[c][2] for c in range(K)])).reshape(K, len(T)); rows = []
    for c in range(K):
        n, k, p = resM[c]; share = k / max(n, 1); fold = share / (bgM / NM)
        o = np.lexsort((-share, p))[:5]
        for j in o: rows.append(dict(cluster=f"c{c}", term=T[j], name=name[T[j]], k=int(k[j]), n=n, share=share[j], fold=fold[j], fdr=allq[c, j]))
    mondo5 = pd.DataFrame(rows)
    # mechanism terms
    MT = sorted(mt.term_id[mt.under_mech_root.astype(str) == "True"]); mti = {t: i for i, t in enumerate(MT)}
    sub = mdt[mdt.term_id.isin(mti)]
    XMe = sp.csr_matrix((np.ones(len(sub), np.int32), (sub.disease_id.map(idx).to_numpy(), sub.term_id.map(mti).to_numpy())), shape=(N, len(MT)))
    resMe, bgMe, NMe = enrich(XMe, mapped, lab, K, min_k=3)
    allq = bh(np.concatenate([resMe[c][2] for c in range(K)])).reshape(K, len(MT)); rows = []
    for c in range(K):
        n, k, p = resMe[c]; share = k / max(n, 1); fold = share / (bgMe / NMe)
        for j, t in enumerate(MT):
            rows.append(dict(cluster=f"c{c}", term=t, name=name[t], k=int(k[j]), n=n, bg=int(bgMe[j]), N_annot=int(NMe), share=share[j], fold=fold[j], p=p[j], fdr=allq[c, j]))
    mech = pd.DataFrame(rows)
    neo = [t for t, n in name.items() if n == "cancer or benign tumor"]; assert len(neo) == 1
    neo_d = np.array([neo[0] in dset.get(d, ()) for d in ids])
    neo_share = {c: float(neo_d[(lab == c) & mapped].mean()) if ((lab == c) & mapped).any() else 0.0 for c in range(K)}
    return body, mondo5, mech, neo_share


def gene_edges_public():
    sys.path.insert(0, str(ROOT / "src"))
    from disatlas.graph import raw_build
    return raw_build.gene_edges()


def go_tables(ids, lab, K, ge, g2, gt):
    N = len(ids); idx = {d: i for i, d in enumerate(ids)}
    gt = gt.set_index("term_id"); gt_anc = {t: (set(a.split("|")) if isinstance(a, str) and a else set()) for t, a in gt.ancestors.items()}
    g2 = g2[(g2.namespace == "biological_process") & g2.go_id.isin(gt.index)]
    gene_terms = defaultdict(set)
    for g, t in zip(g2.gene, g2.go_id): gene_terms[g] |= gt_anc[t] | {t}
    gsize = defaultdict(int)
    for g, v in gene_terms.items():
        for t in v: gsize[t] += 1
    GT = sorted(t for t, n in gsize.items() if 10 <= n <= 500 and gt.namespace[t] == "biological_process"); gti = {t: i for i, t in enumerate(GT)}
    ag = ge[ge.disease_id.isin(idx)]
    dgo = defaultdict(set)
    for d, g in zip(ag.disease_id, ag.gene_symbol):
        for t in gene_terms.get(g, ()):
            if t in gti: dgo[d].add(gti[t])
    r_, c_ = [], []
    for d, v in dgo.items():
        for j in v: r_.append(idx[d]); c_.append(j)
    XG = sp.csr_matrix((np.ones(len(r_), np.int32), (r_, c_)), shape=(N, len(GT)))
    availG = np.zeros(N, bool); availG[[idx[d] for d in dgo]] = True
    gene_dis = np.zeros(N, bool); gene_dis[[idx[d] for d in ag.disease_id.unique()]] = True
    resG, bgG, NG = enrich(XG, availG, lab, K, min_k=3)
    allq = bh(np.concatenate([resG[c][2] for c in range(K)])).reshape(K, len(GT)); rows = []; cov = {}
    for c in range(K):
        n, k, p = resG[c]; share = k / max(n, 1); fold = share / np.maximum(bgG / NG, 1e-12)
        o = np.lexsort((-share, p))[:5]; cov[c] = round(float(gene_dis[lab == c].mean()), 3)
        for rank, j in enumerate(o):
            rows.append(dict(cluster=f"c{c}", rank=rank + 1, term=GT[j], name=gt.name[GT[j]], k=int(k[j]), n=n, share=share[j], fold=fold[j], fdr=allq[c, j],
                             organ_level=bool(gt.organ_level[GT[j]])))
    return pd.DataFrame(rows), XG, ag, cov


def pair_tables(ids, lab, K, XG, ag, hp):
    N = len(ids); idx = {d: i for i, d in enumerate(ids)}
    def bin_mat(df, cd, cx):
        df = df[df[cd].isin(idx)]; xs = {x: i for i, x in enumerate(df[cx].unique())}
        M = sp.csr_matrix((np.ones(len(df), np.float32), (df[cd].map(idx).to_numpy(), df[cx].map(xs).to_numpy())), shape=(N, len(xs)))
        M.data[:] = 1; return M.tocsr()
    CH = dict(hpo=bin_mat(hp, "disease_id", "hpo_id"), gene=bin_mat(ag, "disease_id", "gene_symbol"), gobp=XG.astype(np.float32))
    rng = np.random.default_rng(20261102)
    def pair_stats(M, ii, jj):
        a = M[ii]; b = M[jj]; inter = np.asarray(a.multiply(b).sum(1)).ravel()
        sa = np.asarray(a.sum(1)).ravel(); sb = np.asarray(b.sum(1)).ravel(); ok = (sa > 0) & (sb > 0)
        jac = inter[ok] / (sa[ok] + sb[ok] - inter[ok]); return ok.sum(), (inter[ok] > 0).mean(), jac.mean()
    ri = rng.integers(0, N, 300000); rj = rng.integers(0, N, 300000); kk = ri != rj; ri, rj = ri[kk], rj[kk]
    base = {ch: dict(zip(("n", "sh", "jac"), pair_stats(M, ri, rj))) for ch, M in CH.items()}
    rows = []
    for c in range(K):
        mem = np.where(lab == c)[0]; iu, ju = np.triu_indices(len(mem), 1)
        if len(iu) > 400000:
            s = rng.choice(len(iu), 400000, replace=False); iu, ju = iu[s], ju[s]
        ii, jj = mem[iu], mem[ju]; r = dict(cluster=f"c{c}")
        for ch, M in CH.items():
            n_ok, sh, jac = pair_stats(M, ii, jj)
            r[f"pairs_{ch}_share_fold"] = round(float(sh) / base[ch]["sh"], 2); r[f"pairs_{ch}_jaccard_fold"] = round(float(jac) / base[ch]["jac"], 2)
        rows.append(r)
    return pd.DataFrame(rows)


def compute(gene_edges=None):
    ids, lab, K, M, top, mt, mdt = load_inputs()
    body, mondo5, mech, neo_share = mondo_tables(ids, lab, K, M, top, mt, mdt)
    ge = gene_edges_public() if gene_edges is None else gene_edges
    g2 = pd.read_parquet(ROOT / "data/raw/specific_channels/gene2go.parquet")
    gt = pd.read_csv(ROOT / "data/raw/go/go_bp_terms.tsv.gz", sep="\t", keep_default_na=False)
    gt["organ_level"] = gt.organ_level.astype(str) == "True"
    go5, XG, ag, cov = go_tables(ids, lab, K, ge, g2, gt)
    hp = pd.read_parquet(ROOT / "data/raw/hpo/disease_phenotype.parquet")
    pr = pair_tables(ids, lab, K, XG, ag, hp)
    ann = pd.DataFrame(dict(cluster=[f"c{c}" for c in range(K)], size=[int((lab == c).sum()) for c in range(K)],
                            neoplasm_share_direct=[neo_share[c] for c in range(K)]))
    ann = ann.merge(pr[["cluster"] + STORED_PAIR_COLS], on="cluster"); ann["gene_coverage"] = [cov[c] for c in range(K)]
    return dict(bodysystem_by_cluster_long=body, mondo_terms_top5_by_cluster=mondo5, mechanism_terms_by_cluster=mech,
                go_bp_top5_by_cluster=go5, cluster_annotation_rebuilt=ann)



STORED_PAIR_COLS = ["pairs_hpo_jaccard_fold", "pairs_hpo_share_fold", "pairs_gene_share_fold", "pairs_gobp_jaccard_fold"]


def classify(BL, ME, G5, M5, A, so=0.60, sm=0.15):
    """Theme of each cluster (rule of tables/SI/table_cluster_themes.py, verbatim logic); returns a Series indexed by cluster."""
    FOLD_ORGAN, FOLD_MECH, FDR, NEO_MIN = 2.0, 10.0, 0.05, 0.40
    order = A.cluster.tolist(); sg = sm + 0.10
    BLo = BL[BL.category != "syndromic disease"]
    organ = BLo.loc[BLo.groupby("cluster", sort=False).share.idxmax()].set_index("cluster")
    G5 = G5.sort_values(["cluster", "rank"])
    gocell = G5[~G5.organ_level.astype(bool)].groupby("cluster", sort=False).head(1).set_index("cluster")
    neo5 = M5[M5.name == "cancer or benign tumor"].set_index("cluster").share
    out = {}
    for c in order:
        o = organ.loc[c]; anat = bool((o.share >= so) and (o.fold >= FOLD_ORGAN) and (o.fdr < FDR))
        g = ME[ME.cluster == c]; ok = g[(g.share >= sm) & (g.fold >= FOLD_MECH) & (g.fdr < FDR)]; m1 = len(ok) > 0
        gg = gocell.loc[c] if c in gocell.index else None
        mg = bool(gg is not None and (gg.share >= sg) and (gg.fold >= FOLD_MECH) and (gg.fdr < FDR))
        isneo = float(neo5.get(c, 0.0)) >= NEO_MIN
        out[c] = ("neoplasm" if isneo and not m1 and not anat else "mixed" if anat and (m1 or mg) else
                  "anatomy" if anat else "process" if (m1 or mg) else "unresolved")
    return pd.Series(out)


if __name__ == "__main__":
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    for k, v in compute().items(): v.to_parquet(out / f"{k}.parquet", index=False)

