"""Maintainer only: stage the frozen model outputs of the manuscript into data/zenodo/model_outputs.

Not part of the referee pipeline (reads the authors' internal run directories). Every output is either a
byte copy of the file the manuscript was computed from or a column/row subset of it (float values are
copied, never recomputed). Writes model_outputs/SHA256SUMS and MANIFEST.tsv.

    DISATLAS_STAGE_SOURCES=sources.json python scripts/maintainer/stage_model_outputs.py --out data/zenodo/model_outputs
"""
import argparse, hashlib, json, os, shutil
import numpy as np
import pandas as pd

SEEDS = (4707, 4708, 4709, 4710, 4711)
# Source run directories are site-specific and not part of the release: a JSON file named by
# $DISATLAS_STAGE_SOURCES maps GA_MIRROR, V32, GA, CR, CG, K5, HY to the authors' run outputs.
_SRC = json.load(open(os.environ["DISATLAS_STAGE_SOURCES"])) if os.environ.get("DISATLAS_STAGE_SOURCES") else {}
GA_MIRROR, V32, GA, CR, CG, K5, HY = (_SRC.get(k, "") for k in ("GA_MIRROR", "V32", "GA", "CR", "CG", "K5", "HY"))
DIS = ("MONDO:", "ORPHANET:", "EFO:", "MESH:", "OMIM:")
PQ = dict(index=False, compression="zstd")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); a = ap.parse_args()
    O = a.out; man = []

    def put(rel, src=None, df=None, note=""):
        p = os.path.join(O, rel); os.makedirs(os.path.dirname(p), exist_ok=True)
        if df is not None:
            df.to_parquet(p, **PQ)
        else:
            shutil.copyfile(src, p)
        man.append(dict(file=rel, bytes=os.path.getsize(p), source=(src or "") if df is None else "subset", note=note))

    # 1. embeddings: gene-anchored primary model (5 seeds); the v3.2 embeddings ship once, under zenodo/model
    ids13 = json.load(open(K5)); assert len(ids13) == 13021
    for s in SEEDS:
        src = os.path.realpath(f"{GA_MIRROR}/node2vec_s{s}.npz"); z = np.load(src)
        assert z["X"].shape == (13021, 32) and z["X"].dtype == np.float32 and set(z["ids"].astype(str)) == set(ids13)
        put(f"embeddings/gene_anchored/embedding_s{s}.npz", src, note="primary model (gene-anchored), ids + X float32")

    # 2. full 13,021 x 13,021 five-seed cosine matrix of the primary model (row order = ids_13021.json)
    put("cosine_matrix/cosine_ens5_13021.npy", f"{GA}/out/matrix/C_cos_ens5_13021.npy", note="float32, mean of 5 per-seed cosines")
    put("cosine_matrix/ids_13021.json", K5, note="row/column order of the matrix")
    put("cosine_matrix/matrix_meta.json", f"{GA}/out/matrix/matrix_meta.json", note="pair count and cosine quantiles")

    # 3. IxIDN-Matched pair scores (4,364 pairs): score columns only, no names, rationale or drug evidence
    keep = ["blind_id", "a", "b", "y", "pool", "split_a", "split_b"]
    sc = [f"c_s{s}" for s in SEEDS] + ["c"] + [f"ndc_s{s}" for s in SEEDS] + ["nd_c"]
    for arm, src in (("gene_anchored", f"{GA}/out/core/ixidn4364_cosine_scores.tsv.gz"),):
        T = pd.read_csv(src, sep="\t"); assert len(T) == 4364 and T.y.sum() == 2182
        put(f"scores/ixidn_matched_cosine_{arm}.parquet", df=T[keep + sc],
            note=f"{arm}: per-seed and 5-seed cosine (c), node-disjoint cosine (nd_c); from {src}")

    # 4. node-disjoint placement cosines per fold and seed (primary model)
    rows = []
    for f in range(5):
        for s in SEEDS:
            z = np.load(f"{GA}/out/nd_place/nd_cos_fold{f}_s{s}.npz", allow_pickle=True)
            rows.append(pd.DataFrame(dict(blind_id=z["blind_id"].astype(str), fold=f, seed=s, held_a=z["held_a"], held_b=z["held_b"],
                                          cos_rmref=z["cos_rmref"].astype(np.float64))))
    put("scores/node_disjoint_placement_cosine_gene_anchored.parquet", df=pd.concat(rows, ignore_index=True),
        note="cosine of the re-placed held-out endpoint per fold/seed (cos_rmref); nd_c = mean over held folds of the seed mean")

    # 5. full IxIDN benchmark pair scores (all pairs, incl. unmatched), primary model
    src = f"{GA}/benchmark/out/pair_scores_full_geneanchored.tsv.gz"; T = pd.read_csv(src, sep="\t")
    put("scores/ixidn_full_pair_scores_gene_anchored.parquet", df=T,
        note="all IxIDN pairs: node-disjoint cosine (ND_xF, ND_rel) and graph baselines (svd, cn, aa, degree)")

    # 6. clustering (Data availability: per-disease cluster table and cluster summary)
    M = pd.read_csv(f"{CG}/themes/members_k55_13021.tsv.gz", sep="\t").rename(columns={"mondo_id": "disease_id"})
    L = np.load(f"{GA}/fig5_redraw/k55ga/labels_k55_13021.npz", allow_pickle=True)
    xy = np.load(f"{GA}/fig5_redraw/k55ga/isomap_cos_15_xy_13021.npy")
    assert len(M) == 13021 and M.cluster.nunique() == 55 and xy.shape == (13021, 2)
    lk = [k for k in L.files]; man.append(dict(file="(check)", bytes=0, source="labels npz keys", note=",".join(lk)))
    pos = {d: i for i, d in enumerate(ids13)}; M["isomap_x"] = xy[M.disease_id.map(pos).values, 0]; M["isomap_y"] = xy[M.disease_id.map(pos).values, 1]
    put("clusters/disease_clusters_k55.parquet", df=M[["disease_id", "label", "cluster", "isomap_x", "isomap_y"]],
        note="K = 55 k-medoids cluster of each of the 13,021 diseases and its Figure 5 map coordinate")
    put("clusters/cluster_summary_k55.parquet", df=pd.read_csv(f"{CG}/themes/themes_k55.csv"),
        note="one row per cluster: size, theme, organ-system / mechanism / GO enrichment")

    # 7. P_full retrain reference: pfull_reference/pfull_ixidn_cosine_5seed.parquet holds the IxIDN-Matched pair
    #    columns (blind_id, a, b, label, pool; same 4,364 rows as hybrid_pairs_ens5.tsv) and the per-seed cosines
    #    Pfull_cos_s<seed> of the authors' P_full retrain. It was staged by hand from that run; no decoder column is shipped.

    pd.DataFrame(man).to_csv(os.path.join(O, "MANIFEST.tsv"), sep="\t", index=False)
    with open(os.path.join(O, "SHA256SUMS"), "w") as f:
        for r in sorted(x["file"] for x in man if x["file"] != "(check)"):
            f.write(f"{sha(os.path.join(O, r))}  {r}\n")
    print("STAGED", len(man))


if __name__ == "__main__":
    main()
