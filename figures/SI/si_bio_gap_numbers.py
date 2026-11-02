"""Remaining SI-biology numbers: explainability case scores, single-component compositionality ranks,
single-gene limitation census, mechanism-class positive controls, approved-repurposing benchmark composition.
Writes SI_bio_gap_numbers.json and SI_explainability_case_scores.tsv next to this script."""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; _sys.path.insert(0, str(HERE.parents[1] / "src"))
import numpy as np, pandas as pd
from disatlas.sibio_common import DERIVED, model_dir, write_numbers
ZEN = HERE.parents[1] / "data" / "zenodo"
derived = lambda t: DERIVED / t
r3 = lambda x: round(float(x), 3)
N = {}
# ---- explainability: three literature pairs on the released ensemble (all-pairs score file)
z = np.load(ZEN / "scores" / "disatlas_all_pairs_ensemble.npz"); ids = list(z["ids"]); p = z["p_upper"]; n = len(ids)
ix = {d: i for i, d in enumerate(ids)}; ps = np.sort(p)
def pair(a, b):
    i, j = sorted((ix[a], ix[b])); return float(p[i * n - i * (i + 1) // 2 + (j - i - 1)])
CASES = [("Case1", "MONDO:0007947", "Marfan syndrome", "MONDO:0005396", "thoracic aortic aneurysm"),
         ("Case2", "MONDO:0007452", "maturity-onset diabetes of the young type 1", "MONDO:0005148", "type 2 diabetes mellitus"),
         ("Case3", "MONDO:0012517", "Gaucher disease due to saposin C deficiency", "MONDO:0859183", "Parkinson disease 24, autosomal dominant")]
rows = []
for c, a, la, b, lb in CASES:
    v = pair(a, b); rows.append(dict(case=c, id_a=a, label_a=la, id_b=b, label_b=lb, score=round(v, 4),
                                     percentile_all_pairs=round(100 * np.searchsorted(ps, v) / len(ps), 2), n_diseases=n))
T = pd.DataFrame(rows); T.to_csv(HERE / "SI_explainability_case_scores.tsv", sep="\t", index=False)
N["explainability"] = {r.case: dict(score=r3(r.score), pct=round(r.percentile_all_pairs, 1)) for r in T.itertuples()}
N["explainability"]["n_diseases"] = n
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
# ---- limitation: single-gene pairs census (best rank beyond 100 for both diseases = not recovered)
C = pd.read_csv(derived("si_limitations") / "single_gene_pairs_census.tsv.gz", sep="\t")
fail = C.best_rank > 100; broad_only = C.in_broad & ~C.in_spec
N["single_gene"] = dict(n_pairs=len(C), fail_pct=round(100 * fail.mean()),
                        umbrella_fail_pct=[round(100 * fail[C.UMBRELLA].mean()), round(100 * fail[~C.UMBRELLA].mean())],
                        hub_fail_pct=[round(100 * fail[C.HUB].mean()), round(100 * fail[~C.HUB].mean())],
                        broad_only_share_of_failures_pct=round(100 * broad_only[fail].mean()))
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
N["gene_anchor_ablation"] = dict(top100_pct=[round(100 * GAS.classic_top100_released, 1), round(100 * GAS.classic_top100_arm, 1)],
    targets=[int(GAS.target_pairs_score_up), int(GAS.target_pairs_n)], fda_lost=[int(GAS.fda_dev_pairs_lost), int(GAS.fda_dev_pairs_n)],
    ixidn_delta=[round(GAS.ixidn_auroc_delta, 3), round(GAS.ixidn_auroc_delta_lo, 3), round(GAS.ixidn_auroc_delta_hi, 3)])
write_numbers(__file__, N, "SI_bio_gap")
print(N)
