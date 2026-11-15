"""Remaining SI-biology numbers: single-component compositionality ranks,
single-gene limitation census, mechanism-class positive controls, approved-repurposing benchmark composition.
Writes SI_bio_gap_numbers.json next to this script. (Supplementary Figure S26 case scores: si_explainability.py.)"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; _sys.path.insert(0, str(HERE.parents[1] / "src"))
import numpy as np, pandas as pd
from disatlas.sibio_common import DERIVED, model_dir, write_numbers
ZEN = HERE.parents[1] / "data" / "zenodo"
derived = lambda t: DERIVED / t
r3 = lambda x: round(float(x), 3)
N = {}
# ---- compositionality: rank of the target when ranked from one component alone (same metric as the barycentre rank)
REL = model_dir(); nodes = open(REL / "nodes.txt").read().split(); idn = {d: i for i, d in enumerate(nodes)}
XS = []
for s in [4707, 4708, 4709, 4710, 4711]:
    X = np.load(REL / "embeddings" / f"node2vec_s{s}.npz", allow_pickle=True)["X"].astype(np.float64)
    XS.append(X / np.linalg.norm(X, axis=1, keepdims=True))
def bary_rank(target, comps, exclude=None):
    """Rank of target by cosine to the (per-seed) mean of comps; every disease in `exclude` (default comps) is removed."""
    S_ = [idn[c] for c in comps]; E_ = [idn[c] for c in (exclude or comps)]
    sim = np.mean([X @ (X[S_].mean(0) / np.linalg.norm(X[S_].mean(0))) for X in XS], 0)
    keep = np.ones(len(nodes), bool); keep[E_] = False
    return int((sim[keep] > sim[idn[target]]).sum()) + 1
EX = {"MCTD": ("MONDO:0005854", ["MONDO:0007915", "MONDO:0005100", "MONDO:0010030", "MONDO:0019434", "MONDO:0019127"]),
      "thymoma": ("MONDO:0006456", ["MONDO:0009688", "MONDO:0016367"])}
N["compositionality_single"] = {k: {cc: bary_rank(t, [cc], exclude=comps) for cc in comps} | {"barycentre": bary_rank(t, comps)} for k, (t, comps) in EX.items()}
# ---- limitation: single-gene pairs census. Best partner rank = the better of the two ranks each disease gives the other
# among all other diseases ordered by 5-seed mean cosine (1 = highest); beyond 100 for both diseases = not recovered.
# The census file supplies the pairs and the UMBRELLA / HUB / layer flags; its best_rank column is not used.
C = pd.read_csv(derived("si_limitations") / "single_gene_pairs_census.tsv.gz", sep="\t")
def cosine_best_rank(Xs, idx, A, B):
    a = np.array([idx[x] for x in A]); b = np.array([idx[x] for x in B])
    M = np.mean([X @ X.T for X in Xs], 0); np.fill_diagonal(M, -np.inf)
    s = M[a, b]
    return np.minimum((M[a] > s[:, None]).sum(1), (M[b] > s[:, None]).sum(1)) + 1
best = cosine_best_rank(XS, idn, C.a, C.b)
fail = best > 100; broad_only = C.in_broad & ~C.in_spec
N["single_gene"] = dict(n_pairs=len(C), fail_pct=round(100 * fail.mean()),
                        umbrella_fail_pct=[round(100 * fail[C.UMBRELLA].mean()), round(100 * fail[~C.UMBRELLA].mean())],
                        hub_fail_pct=[round(100 * fail[C.HUB].mean()), round(100 * fail[~C.HUB].mean())],
                        broad_only_share_of_failures_pct=round(100 * broad_only[fail].mean()))
# gene-anchored arm (data/zenodo/model_outputs/embeddings/gene_anchored): same census, same rank rule
GA = ZEN / "model_outputs" / "embeddings" / "gene_anchored"
ga_nodes = [str(x) for x in np.load(GA / "embedding_s4707.npz", allow_pickle=True)["ids"]]
ga_idx = {d: i for i, d in enumerate(ga_nodes)}
GXS = []
for s_ in [4707, 4708, 4709, 4710, 4711]:
    X = np.load(GA / f"embedding_s{s_}.npz", allow_pickle=True)["X"].astype(np.float64)
    GXS.append(X / np.linalg.norm(X, axis=1, keepdims=True))
ga_best = cosine_best_rank(GXS, ga_idx, C.a, C.b)
TOP100 = [round(100 * float((best <= 100).mean()), 1), round(100 * float((ga_best <= 100).mean()), 1)]
# ---- mechanism-class positive controls (frozen degree strata)
M = pd.read_csv(derived("si_mechanism_timelines") / "mechanism_class_controls.tsv", sep="\t")
N["mechanism_classes"] = dict(n_classes=len(M), n_cases=int(M.n_cases.sum()), n_recovered=int(M.k_recovered_majority.sum()),
                              median_pct_range=[float(M.median_case_pct_cm.min()), float(M.median_case_pct_cm.max())])
# ---- approved-repurposing benchmark composition + channel-comparison universe
A = pd.read_csv(derived("si_fda") / "benchmark_pairs_split.tsv", sep="\t")
N["fda_benchmark"] = dict(n_pairs=len(A), n_diseases=len(set(A.mondo_a) | set(A.mondo_b)),
                          **{f"n_{k}": int(v) for k, v in A.stratum.value_counts().items()},
                          **{f"n_split_{k}": int(v) for k, v in A.split.value_counts().items()},
                          n_hard_complement=int(((A.stratum == "hard") & (A.split == "complement_six")).sum()))
U = pd.read_csv(derived("si_fda") / "channel_arm_universe.tsv", sep="\t").set_index("arm").n_nodes
N["fda_channel_universe"] = int(U["full@broad_off_nodes"]); assert U["full@broad_off_nodes"] == U["spec_off@broad_off_nodes"]
# rare-shared-gene link ablation (limitations): five-seed retrain summary, data/derived/si_limitations/README.md
GAS = pd.read_csv(DERIVED / "si_limitations" / "gene_anchor_ablation_summary.tsv", sep="\t").set_index("quantity").value
# fda_lost and ixidn_delta are copied from the precomputed summary (decoder-era readout; the final SI value, 3 of 159,
# is not reproducible from release files: see the fda_lost record in tests/parity/sibio_limitations.json)
N["gene_anchor_ablation"] = dict(top100_pct=TOP100,
    targets=[int(GAS.target_pairs_score_up), int(GAS.target_pairs_n)], fda_lost=[int(GAS.fda_dev_pairs_lost), int(GAS.fda_dev_pairs_n)],
    ixidn_delta=[round(GAS.ixidn_auroc_delta, 3), round(GAS.ixidn_auroc_delta_lo, 3), round(GAS.ixidn_auroc_delta_hi, 3)])
write_numbers(__file__, N, "SI_bio_gap")
print(N)
