"""Supplementary Table SI_cluster_controls: positive and negative control disease pairs at K = 55.

Pair names, roles and MONDO ids: data/derived/si_cluster/SI_controls_table.tsv. Everything else is computed from the
released model outputs: same cluster from the K = 55 cosine k-medoids labels
(data/zenodo/model_outputs/clusters/disease_clusters_k55.parquet); pair cosine score from the five-seed ensemble cosine
matrix (data/zenodo/model_outputs/cosine_matrix/cosine_ens5_13021.npy); percentile of that score among all 84,766,710
disease pairs (i < j); best partner rank = the better of the two ranks each disease gives the other among its 13,020
other diseases ordered by cosine (1 = highest).
Outputs next to this script: table_cluster_controls.{md,tsv}, table_cluster_controls_numbers.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MO = ROOT / "data/zenodo/model_outputs"
ids = json.load(open(MO / "cosine_matrix/ids_13021.json")); n = len(ids); assert n == 13021
ix = {d: i for i, d in enumerate(ids)}
C = np.load(MO / "cosine_matrix/cosine_ens5_13021.npy", mmap_mode="r")
K = pd.read_parquet(MO / "clusters/disease_clusters_k55.parquet").set_index("disease_id")
lab = K.loc[ids, "cluster"].to_numpy(); assert len(np.unique(lab)) == 55
T = pd.read_csv(ROOT / "data/derived/si_cluster/SI_controls_table.tsv", sep="\t")
assert len(T) == 17
iu = np.triu_indices(n, 1)
allc = np.sort(np.asarray(C, dtype=np.float32)[iu]); NP = len(allc); assert NP == n * (n - 1) // 2


def rank_of(i, j):
    r = np.asarray(C[i], dtype=np.float32).copy(); r[i] = -np.inf
    return int((r > r[j]).sum()) + 1


rows = []
for t in T.itertuples():
    i, j = ix[t.id_a], ix[t.id_b]
    s = float(C[i, j])
    rows.append(dict(pair=t.pair, role=t.role, tag=t.tag, same="yes" if lab[i] == lab[j] else "no", score=s,
                     pct=100.0 * np.searchsorted(allc, np.float32(s), side="right") / NP,
                     best_rank=min(rank_of(i, j), rank_of(j, i))))
R = pd.DataFrame(rows)
R = R.assign(o=(R.role != "positive").astype(int)).sort_values(["o", "score"], ascending=[True, False]).drop(columns="o")
cols = ["Control pair", "Role", "Same cluster at K = 55", "DisAtlas pair cosine score", "Percentile among all pairs",
        "Best partner rank"]
out = [[r.pair, r.role, r.same, f"{r.score:.3f}", f"{r.pct:.1f}", f"{r.best_rank:,}"] for r in R.itertuples()]
pd.DataFrame(out, columns=cols).to_csv(HERE / "table_cluster_controls.tsv", sep="\t", index=False)
md = ["| " + " | ".join(cols) + " |", "|:--|:--|:--|--:|--:|--:|"] + ["| " + " | ".join(r) + " |" for r in out]
(HERE / "table_cluster_controls.md").write_text("\n".join(md) + "\n", encoding="utf-8")
pos, neg = R[R.role == "positive"], R[R.role == "negative"]
N = dict(n_pairs_all=int(NP), median_cosine_all_pairs=round(float(np.median(allc)), 3),
         n_pos=len(pos), n_neg=len(neg), pos_same_cluster=int((pos.same == "yes").sum()),
         neg_separated=int((neg.same == "no").sum()), pos_top100=int((pos.best_rank <= 100).sum()),
         pos_best_rank_range=[int(pos.best_rank.min()), int(pos.best_rank.max())],
         neg_best_rank_range=[int(neg.best_rank.min()), int(neg.best_rank.max())],
         controls=R[["tag", "role", "same", "score", "pct", "best_rank"]].round(4).to_dict("records"))
(HERE / "table_cluster_controls_numbers.json").write_text(json.dumps(N, indent=1) + "\n", encoding="utf-8")
print(json.dumps({k: v for k, v in N.items() if k != "controls"}))
