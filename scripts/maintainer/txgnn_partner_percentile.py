"""Maintainer-only: TxGNN column of the rare-common documented-partner recovery (SI, rare-common section).

Not needed to reproduce any result: its output ships as data/derived/si_rare_common/partner_percentile_per_pair_txgnn.tsv
and is read by figures/SI/si_cluster_rarecommon_numbers.py. Rerunning it needs the published TxGNN Explorer checkpoint,
the PrimeKG full_graph split and the TxGNN benchmark scoring helpers (score_pairs.py: idmap, mondo_lookup, map_pairs,
distmult, CKPT, SPLIT, ET), none of which are redistributed here; point DISATLAS_TXGNN_SCORER at the directory holding
score_pairs.py and run it in a TxGNN environment (torch).

Scoring is the IxIDN-Matched TxGNN comparator unchanged: DistMult disease_disease decoder of the published checkpoint,
MONDO ids mapped to PrimeKG disease nodes. A pair is scored only if both diseases map to nodes with non-zero embeddings
and do not collapse onto one node; every other pair is ranked last (score below every mapped pair; ties share the
mid-rank). Validation: the IxIDN-Matched benchmark pairs are re-scored and must equal the shipped TxGNN comparator
column (same unmapped pattern, max |diff| < 1e-6).

Candidates and percentile formula are those of partner_percentile_per_pair.tsv: the anchor's score against every
disease of the other class in data/derived/si_rare_common/partner_universe.tsv (disease_id in score-matrix order,
disease_class), percentile = (#below + 0.5 * (#ties - 1)) / (n_candidates - 1).

usage: DISATLAS_TXGNN_SCORER=<dir> python scripts/maintainer/txgnn_partner_percentile.py OUT_DIR
"""
import os, sys, pickle, json
import numpy as np, pandas as pd, torch

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
L = os.path.join(REPO, "data", "derived", "si_rare_common")
BEN = os.path.join(REPO, "data", "derived", "fig6_comparison", "per_run_scores.tsv")
sys.path.insert(0, os.environ["DISATLAS_TXGNN_SCORER"])
import score_pairs as SPm  # noqa: E402  (main() is guarded; only helpers are used)
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)

df = pd.concat([pd.read_csv(os.path.join(SPm.SPLIT, f"{s}.csv")) for s in ["train", "valid", "test"]])
dis = SPm.idmap(df, "disease"); lookup = SPm.mondo_lookup(dis)
sd = torch.load(os.path.join(SPm.CKPT, "model.pt"), map_location="cpu")
etypes = sorted({tuple(r) for r in df[["x_type", "relation", "y_type"]].drop_duplicates().values})
assert len(etypes) == sd["w_rels"].shape[0]
w = sd["w_rels"].numpy()[etypes.index(SPm.ET)]
emb = pickle.load(open(os.path.join(SPm.CKPT, "node_emb.pkl"), "rb"))
E = emb["disease"]; E = E.numpy() if hasattr(E, "numpy") else np.asarray(E)
assert E.shape[0] == len(set(dis.values()))
del df


def score(a_ids, b_ids):
    """the benchmark scoring path for a list of MONDO pairs; NaN where not mapped"""
    pr = pd.DataFrame({"a": a_ids, "b": b_ids})
    ia, ib, coll = SPm.map_pairs(pr, lookup)
    ok = (ia >= 0) & (ib >= 0)
    ea, eb = E[np.where(ok, ia, 0)], E[np.where(ok, ib, 0)]
    ok &= (np.abs(ea).sum(1) > 0) & (np.abs(eb).sum(1) > 0)
    ok &= ~coll
    return np.where(ok, SPm.distmult(ea, w, eb), np.nan), ok


# validation against the shipped benchmark TxGNN column
B = pd.read_csv(BEN, sep="\t", usecols=["a", "b", "TxGNN published ckpt|DistMult"])
s_b, ok_b = score(B.a.values, B.b.values); ref = B["TxGNN published ckpt|DistMult"].values
val = dict(n_bench=len(B), n_mapped_rescored=int(ok_b.sum()), n_mapped_benchmark=int(np.isfinite(ref).sum()),
           same_nan_pattern=bool((np.isfinite(ref) == ok_b).all()), max_abs_diff=float(np.nanmax(np.abs(s_b - ref))))
assert val["same_nan_pattern"] and val["max_abs_diff"] < 1e-6, val

U = pd.read_csv(os.path.join(L, "partner_universe.tsv"), sep="\t")
uid = U.disease_id.tolist(); ucls = U.disease_class.values; ui = {d: i for i, d in enumerate(uid)}
PR = pd.read_csv(os.path.join(L, "partner_percentile_per_pair.tsv"), sep="\t")
FLOOR = -1.0  # DistMult score is a sigmoid in (0, 1): every unmapped pair sits below every mapped pair


def pct(row, cand, partner):
    s = row[cand]; v = row[partner]; return ((s < v).sum() + 0.5 * ((s == v).sum() - 1)) / (len(cand) - 1)


recs = []
for _, r in PR.iterrows():
    cl_ = "common" if r.direction.startswith("rare anchor") else "rare"
    a, p_ = ui[r.anchor], ui[r.partner]; cand = np.where((ucls == cl_) & (np.arange(len(uid)) != a))[0]
    s, ok = score([r.anchor] * len(uid), uid)
    recs.append(dict(direction=r.direction, anchor=r.anchor, partner=r.partner,
                     TxGNN=pct(np.where(ok, s, FLOOR), cand, p_), TxGNN_pair_mapped=bool(ok[p_])))
T = pd.DataFrame(recs)
T.to_csv(os.path.join(OUT, "partner_percentile_per_pair_txgnn.tsv"), sep="\t", index=False)
S = pd.read_csv(os.path.join(L, "partner_percentile_per_pair_txgnn.tsv"), sep="\t")
val.update(n_rows=len(T), n_pair_unmapped=int((~T.TxGNN_pair_mapped).sum()),
           equals_shipped=bool(T[["direction", "anchor", "partner", "TxGNN_pair_mapped"]].equals(S[["direction", "anchor", "partner", "TxGNN_pair_mapped"]])
                               and np.allclose(T.TxGNN, S.TxGNN, rtol=0, atol=1e-12)))
print(json.dumps(val, indent=1))
