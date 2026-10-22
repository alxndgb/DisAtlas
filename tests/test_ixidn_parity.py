"""Parity of the rebuilt IxIDN benchmark (code/01_build_ixidn/) with the manuscript and the frozen
reference tables. Run the build steps first (code/01_build_ixidn/run_all.sh), then pytest."""
import os

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
OUT = os.path.join(ROOT, "data", "01_rebuilt")
FROZEN = os.path.join(ROOT, "data", "00_raw", "frozen")
STAGED_MATCHED = os.path.join(ROOT, "data", "06S_final_benchmark", "ixidn_matched_4364.tsv")
SCORES = os.path.join(ROOT, "data", "08S_eval", "ixidn4364_scores.tsv")
NET = os.path.join(ROOT, "data", "03_network")

pytestmark = pytest.mark.skipif(not os.path.exists(os.path.join(OUT, "ixidn_matched.tsv")),
                                reason="run code/01_build_ixidn/run_all.sh first")


def pairs(df, a="a", b="b"):
    return set(map(frozenset, zip(df[a], df[b])))


def test_ixidn_pos():
    e = pd.read_parquet(os.path.join(OUT, "ixidn_pos_edges.parquet"))
    assert len(e) == 43112
    cv = pd.read_parquet(os.path.join(OUT, "disease_coinvestigated_drugs.parquet"))
    assert len(cv) == 285773


def test_ixidn_neg():
    n = pd.read_parquet(os.path.join(OUT, "ixidn_neg_edges.parquet"))
    assert len(n) == 335557
    assert (n.d_sp >= 3).all()
    ev = pd.read_parquet(os.path.join(OUT, "ixidn_eval_pairs.parquet"))
    fz = pd.read_parquet(os.path.join(FROZEN, "ixidn_eval_pairs.parquet"))
    assert len(ev) == len(fz) == 378669 and set(map(tuple, ev.values)) == set(map(tuple, fz.values))


def test_training_strip():
    st = pd.read_parquet(os.path.join(OUT, "ixidn_train_strip_pairs.parquet"))
    fz = pd.read_parquet(os.path.join(FROZEN, "pinned_train_strip_pairs.parquet"))
    core = st[st.source == "ixidn"][["a", "b", "class"]]
    assert set(map(tuple, core.values)) == set(map(tuple, fz.values))
    e = pd.read_parquet(os.path.join(NET, "mechdist_network_edges.parquet.gzip"))
    nodes = pd.read_parquet(os.path.join(NET, "mechdist_network_nodes.parquet.gzip"))
    E = pairs(e, "disease_a", "disease_b")
    S = pairs(st)
    G = nx.Graph()
    G.add_nodes_from(nodes.disease_id)
    G.add_edges_from(tuple(x) for x in E - S)
    assert len(E & S) == 6245
    assert G.number_of_edges() == 325837
    assert sum(1 for v in G if G.degree(v) == 0) == 2


def test_ixidn_matched_equals_staged():
    a = pd.read_csv(os.path.join(OUT, "ixidn_matched.tsv"), sep="\t")
    b = pd.read_csv(STAGED_MATCHED, sep="\t")
    assert len(a) == 4364 and int(a.label.sum()) == 2182
    for c in a.columns:
        x, y = a[c], b[c]
        if x.dtype.kind == "f" or y.dtype.kind == "f":
            assert np.allclose(x.astype(float), y.astype(float), equal_nan=True, rtol=0, atol=1e-12), c
        else:
            assert (x.astype(str) == y.astype(str)).all(), c


def test_headline_auroc():
    a = pd.read_csv(os.path.join(OUT, "ixidn_matched.tsv"), sep="\t")
    s = pd.read_csv(SCORES, sep="\t")
    p = dict(zip(map(frozenset, zip(s.a, s.b)), s.p))
    a["p"] = [p[frozenset(x)] for x in zip(a.a, a.b)]
    assert round(roc_auc_score(a.label, a.p), 3) == 0.765
    t = a[a.pool == "test"]
    assert round(roc_auc_score(t.label, t.p), 3) == 0.782

SHORTCUTS = {  # manuscript appendix, IxIDN-Matched non-degree shortcut baselines (AUROC, 3 d.p.)
    "label_mondo_simgic": 0.542, "label_name_tfidf": 0.506, "gene_count_logsum": 0.504, "gene_coverage_min": 0.500,
}


def test_shortcut_baselines():
    s = pd.read_csv(os.path.join(OUT, "ixidn_matched_shortcut_baselines.tsv"), sep="\t")
    s = s[s.stratum == "all"].set_index("score")
    assert set(s.index) == set(SHORTCUTS)
    for k, v in SHORTCUTS.items():
        assert round(s.auroc[k], 3) == v, (k, s.auroc[k])
    assert s.n["label_mondo_simgic"] == 4357  # 7 pairs without a SimGIC value dropped
    assert round(s.auroc["label_mondo_simgic"], 4) == 0.5418


def test_shortcut_stack_rare_rare():
    # non-degree stack; corrects the manuscript's +0.121 (+0.059 to +0.179), which used a stack with two degree scores
    s = pd.read_csv(os.path.join(OUT, "ixidn_matched_shortcut_stack.tsv"), sep="\t")
    r = s[(s.benchmark == "IxIDN-Matched") & (s.scope == "All pairs") & (s.stratum == "rare-rare")].iloc[0]
    assert (round(r.delta, 3), round(r.delta_lo, 3), round(r.delta_hi, 3)) == (0.180, 0.122, 0.227)
