"""Parity tests for the released data, training graph and ensemble (need data/zenodo; skipped otherwise)."""
import os

import numpy as np
import pandas as pd
import pytest

Z = os.path.join("data", "zenodo")
needs_zenodo = pytest.mark.skipif(not os.path.isdir(os.path.join(Z, "graph")), reason="Zenodo bundle not present")


def test_blinding_tokens_well_formed():
    from disatlas.blinding import is_blinded
    dg = pd.read_parquet("data/raw/channels/broad/disease_gene.parquet")
    b = dg.target.map(is_blinded)
    assert b.any() and dg.target[b].str.fullmatch(r"T_[0-9a-f]{12}").all()
    assert not dg.target[~b].str.startswith(("T_", "M_")).any()
    # restricted-source rows carry only group-level source labels
    assert dg.source_group[b].str.contains("literature|tier_rescue|restricted").all()


def test_benchmark_counts():
    m = pd.read_csv("data/raw/benchmarks/ixidn_matched/ixidn_matched_pairs.tsv", sep="\t")
    assert len(m) == 4364 and m.y.sum() == 2182
    f = pd.read_csv("data/raw/benchmarks/fda_repurposing/fda_repurposing_controls.tsv", sep="\t")
    assert len(f) == 1440 and not f.split.eq("holdout").any()
    assert (f.split.eq("dev") & f.fda_verified_pair & f.in_model).sum() == 159


@needs_zenodo
def test_graph_rebuild_identical():
    from disatlas.graph import build
    arm, tr, _, stats = build.build("data")
    assert len(arm) == 498155 and pd.concat([arm.u, arm.v]).nunique() == 13021
    build.check_against(arm, tr, os.path.join(Z, "graph"))


@needs_zenodo
def test_ensemble_reproduces_ixidn_matched_scores():
    from disatlas.model import load_ensemble
    m = load_ensemble()
    b = pd.read_csv("data/raw/benchmarks/ixidn_matched/ixidn_matched_pairs.tsv", sep="\t")
    s = pd.read_csv(os.path.join(Z, "scores", "ixidn_matched_disatlas_full_scores.tsv"), sep="\t")
    p = m.score_pairs(b.disease_a, b.disease_b)
    assert np.nanmax(np.abs(p - s.p_ensemble.values)) < 1e-5
