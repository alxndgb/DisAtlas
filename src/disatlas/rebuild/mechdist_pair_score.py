"""Rebuild data/derived/si_mechdist_pair_score/ (Supplementary Figure S6, pair score and the aHUS-PNH example).

    python -m disatlas.rebuild.mechdist_pair_score OUTDIR

* ixidn_matched_cosine_all.tsv  : y, pool, c = 5-seed mean embedding cosine of every IxIDN-Matched pair
                                  (data/derived/fig6_comparison/per_run_scores.tsv, columns DisAtlas|s*).
* schematic_cd_control_cosines.tsv : cosine of the six schematic pairs (released 13,021 x 13,021 ensemble cosine matrix,
                                  3 decimals).
* s6_inputs.json : Youden threshold on the tune split and sensitivity / specificity on the test split (from the cosines);
                  shared-neighbour count, direct-edge flag and hub membership of the aHUS-PNH pair in the training graph
                  (data/raw/training_graph, benchmark and FDA pairs stripped), and the 12 drawn neighbours (all in the shared set).
  Not reproduced: training_graph.hub_deg (shipped 3,036; the staged graph gives 3,006 for the same hub). The 30 missing
  edges of the 'neoplasm' hub are not recoverable from the released inputs (the figure's graph was stripped slightly
  differently from the staged baseline graph). Table label (c).
Test: tests/test_rebuild_mechdist_pair_score.py
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_curve

ROOT = Path(__file__).resolve().parents[3]
REF = ROOT / "data/derived/si_mechdist_pair_score"
MO = ROOT / "data/zenodo/model_outputs/cosine_matrix"


def cosines():
    R = pd.read_csv(ROOT / "data/derived/fig6_comparison/per_run_scores.tsv", sep="\t")
    c = R[[x for x in R.columns if x.startswith("DisAtlas|s")]].mean(axis=1)
    return pd.DataFrame({"y": R.y, "pool": R.pool, "c": c})


def threshold_sens_spec(C):
    t = C[C.pool == "tune"]; fpr, tpr, thr = roc_curve(t.y, t.c); thr = float(thr[(tpr - fpr).argmax()])
    te = C[C.pool == "test"]; pred = te.c >= thr
    return thr, float((pred & (te.y == 1)).sum() / (te.y == 1).sum()), float(((~pred) & (te.y == 0)).sum() / (te.y == 0).sum())


def graph_counts(pair):
    G = pd.read_parquet(ROOT / "data/raw/training_graph/training_graph_edges.parquet")
    deg = pd.concat([G.u, G.v]).value_counts(); a, b = pair
    nb = lambda x: set(G.v[G.u == x]) | set(G.u[G.v == x])
    A, Bn = nb(a), nb(b); sh = A & Bn; hub = deg[list(sh)].idxmax()
    return dict(n_shared=len(sh), hub_deg=int(deg[hub]), direct_AB=bool(b in A), hub_in_shared=bool(hub in sh)), sh


def schematic(ref):
    ids = json.load(open(MO / "ids_13021.json")); ix = {d: i for i, d in enumerate(ids)}
    M = np.load(MO / "cosine_ens5_13021.npy", mmap_mode="r"); out = ref.copy()
    out["cosine_score"] = [round(float(M[ix[a], ix[b]]), 3) for a, b in zip(ref.id_a, ref.id_b)]
    return out


def build(out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    O = json.load(open(REF / "s6_inputs.json")); C = cosines()
    C.to_csv(out / "ixidn_matched_cosine_all.tsv", sep="\t", index=False, float_format="%.6g")
    pd.read_csv(REF / "schematic_cd_control_cosines.tsv", sep="\t").pipe(schematic).to_csv(out / "schematic_cd_control_cosines.tsv", sep="\t", index=False, float_format="%.3f")
    thr, se, sp = threshold_sens_spec(C); T, sh = graph_counts(O["pair"])
    json.dump(dict(pair=O["pair"], training_graph=T, figure_shared_in_shared_set=bool(set(O["figure_shared"]) <= sh), threshold_youden_tune=round(thr, 6),
                   test_sensitivity=se, test_specificity=sp), open(out / "s6_inputs_rebuilt.json", "w"), indent=1)


if __name__ == "__main__":
    build(sys.argv[1])
