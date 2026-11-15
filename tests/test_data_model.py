"""Parity tests for the released data, training graph and ensemble (need data/zenodo; skipped otherwise)."""
import os

import numpy as np
import pandas as pd
import pytest

Z = os.path.join("data", "zenodo")
needs_zenodo = pytest.mark.skipif(not os.path.isdir(os.path.join(Z, "graph")), reason="Zenodo bundle not present")


def test_raw_gene_edges_public_identifiers():
    """The raw disease-gene tables carry public identifiers only (no blinded tokens); the only
    knowledge-graph rows are the 8,343 released Habcd edges."""
    import glob
    for f in glob.glob("data/raw/*/disease_gene*.parquet"):
        g = pd.read_parquet(f)
        col = "gene_id" if "gene_id" in g.columns else [c for c in g.columns if "gene" in c][0]
        assert not g[col].astype(str).str.match(r"^(T|M)_[0-9a-f]{12}$").any(), f
        assert not any("bikg" in c.lower() for c in g.columns), f
    assert len(pd.read_parquet("data/raw/bikg_habcd/habcd_bikg_gene_edges.parquet")) == 8343


def test_benchmark_counts():
    m = pd.read_csv("data/raw/benchmarks/ixidn_matched/ixidn_matched_pairs.tsv", sep="\t")
    assert len(m) == 4364 and m.y.sum() == 2182
    f = pd.read_csv("data/raw/benchmarks/fda_repurposing/fda_repurposing_controls.tsv", sep="\t")
    assert len(f) == 1440 and not f.split.eq("holdout").any()
    assert (f.split.eq("dev") & f.fda_verified_pair & f.in_model).sum() == 159


@needs_zenodo
def test_graph_build_stats():
    """The manuscript training graph is not released (its specific channel was built from BIKG inputs); its
    aggregate build statistics ship in graph/graph_build_stats.json. The rebuildable training graph is the P_full
    graph (tests/test_pfull_graph.py), and its retrain is judged by tests/test_pfull_parity.py."""
    import json
    st = json.load(open(os.path.join(Z, "graph", "graph_build_stats.json")))
    u, s = st["union"], st["strip"]
    assert (u["edges"], u["nodes"]) == (498155, 13021) and u["union_edges"] - u["fda_edges_removed"] == u["edges"]
    assert u["broad"] + u["specific"] - u["both"] == u["edges"]
    assert s["removed"] == s["removed_plus"] + s["removed_minus"] and s["train_edges"] == u["edges"] - s["removed"]


@needs_zenodo
def test_ensemble_reproduces_ixidn_matched_scores():
    """The DisAtlas pair score is the mean over seeds 4707-4711 of the embedding cosine (no decoder); it reproduces
    the shipped per-seed and ensemble IxIDN-Matched cosines of the final (gene-anchored) model."""
    from disatlas.model import load_ensemble, SEEDS
    m = load_ensemble()
    s = pd.read_parquet(os.path.join(Z, "model_outputs", "scores", "ixidn_matched_cosine_gene_anchored.parquet"))
    per = m.score_pairs(s.a, s.b, per_seed=True)
    for k, sd in enumerate(SEEDS):
        assert np.nanmax(np.abs(per[k] - s[f"c_s{sd}"].values)) < 1e-5
    assert np.nanmax(np.abs(m.score_pairs(s.a, s.b) - s.c.values)) < 1e-5


@needs_zenodo
def test_ensemble_reproduces_cosine_matrix_sample():
    """Random sample of the shipped 13,021 x 13,021 ensemble cosine matrix (memory-mapped)."""
    import json
    from disatlas.model import load_ensemble
    d = os.path.join(Z, "model_outputs", "cosine_matrix")
    C = np.load(os.path.join(d, "cosine_ens5_13021.npy"), mmap_mode="r")
    ids = json.load(open(os.path.join(d, "ids_13021.json")))
    ids = ids["ids"] if isinstance(ids, dict) else ids
    rng = np.random.default_rng(0); i = rng.integers(0, len(ids), 2000); j = rng.integers(0, len(ids), 2000)
    p = load_ensemble().score_pairs([ids[x] for x in i], [ids[x] for x in j])
    assert np.nanmax(np.abs(p - np.asarray(C[i, j]))) < 1e-5