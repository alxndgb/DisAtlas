"""Rebuild data/derived/si_explainability/case_scores_cosine.tsv (Supplementary Figure S26 panel a, DisAtlas column) from the
released cosine matrix.

    python -m disatlas.rebuild.case_scores OUTDIR

Input: data/zenodo/model_outputs/cosine_matrix/cosine_ens5_13021.npy (13,021 x 13,021 mean cosine of the five-seed gene-anchored
ensemble) with ids_13021.json. For each of the three rare-common case pairs: cosine (4 d.p.), percentile of the cosine among the
84,766,710 upper-triangle disease pairs (share of pairs with cosine <= the case cosine, 4 d.p.), and the rank of the common disease in
the rare disease's row (1 = highest, the rare disease itself excluded). Disease labels are the MONDO labels of the released K = 55 cluster table.
Test: tests/test_rebuild_case_scores.py.
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
CASES = [("Case1", "Marfan-TAA", "MONDO:0007947", "MONDO:0005396"),
         ("Case2", "MODY1-T2D", "MONDO:0007452", "MONDO:0005148"),
         ("Case3", "GaucherSapC-PD24", "MONDO:0012517", "MONDO:0859183")]


def compute():
    SC = ROOT / "data/zenodo/model_outputs/cosine_matrix"
    ids = json.load(open(SC / "ids_13021.json")); ix = {d: i for i, d in enumerate(ids)}; n = len(ids)
    C = np.load(SC / "cosine_ens5_13021.npy", mmap_mode="r")
    allc = np.sort(np.asarray(C, np.float64)[np.triu_indices(n, 1)])
    T = pd.read_parquet(ROOT / "data/zenodo/model_outputs/clusters/disease_clusters_k55.parquet", columns=["disease_id", "label"])
    lab = dict(zip(T.disease_id, T.label)); rows = []
    for case, tag, a, b in CASES:
        i, j = ix[a], ix[b]; c = float(C[i, j]); rc = np.asarray(C[i], np.float64).copy(); rc[i] = -np.inf
        rows.append(dict(case=case, tag=tag, id_a=a, label_a=lab[a], id_b=b, label_b=lab[b], cosine=round(c, 4),
                         cosine_pct_all_pairs=round(100 * np.searchsorted(allc, c, "right") / len(allc), 4),
                         rank_b_in_row_a=int((rc > rc[j]).sum() + 1), n_diseases=n, n_pairs=len(allc)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    o = Path(sys.argv[1]); o.mkdir(parents=True, exist_ok=True); compute().to_csv(o / "case_scores_cosine.tsv", sep="\t", index=False); print("wrote", o)
