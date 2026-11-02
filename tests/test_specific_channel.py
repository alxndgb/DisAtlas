"""Parity of the specific channel rebuilt from its raw evidence (disatlas.channels.specific)."""
import os

import numpy as np
import pandas as pd
import pytest

from disatlas.channels import specific

CH = os.path.join("data", "raw", "channels")
Z = os.path.join("data", "zenodo")
try:
    specific.raw_path("data", "gene_set_counts.parquet")
    HAVE_RAW = True
except FileNotFoundError:
    HAVE_RAW = False
needs_raw = pytest.mark.skipif(not HAVE_RAW, reason="specific raw evidence not present (Zenodo)")


@pytest.fixture(scope="module")
def rebuilt():
    return specific.build("data")


@needs_raw
def test_links_identical_to_shipped(rebuilt):
    ref = pd.read_parquet(os.path.join(CH, "specific", "disease_module_links.parquet"))
    assert len(rebuilt) == len(ref) == 424662
    assert list(rebuilt.columns) == list(ref.columns)
    for c in ("disease_id", "module_id", "module_order", "scope", "n_fam", "families_group"):
        assert (rebuilt[c].values == ref[c].values).all(), c
    assert np.array_equal(rebuilt.docs_max.values, ref.docs_max.values)


@needs_raw
def test_projection_counts(rebuilt):
    from disatlas.graph import build
    bn = set(pd.read_parquet(os.path.join(CH, "broad", "broad_nodes.parquet")).disease_id)
    spec, _, s = build.specific_channel(rebuilt, bn)
    assert (s["links"], s["diseases"], s["modules"]) == (422804, 5198, 6746)
    assert (s["raw_edges"], s["edges"], s["nodes"]) == (191298, 191019, 5194)


@needs_raw
@pytest.mark.skipif(not os.path.isdir(os.path.join(Z, "graph")), reason="Zenodo bundle not present")
def test_training_graph_from_rebuilt_links(rebuilt):
    from disatlas.graph import build
    bm = os.path.join("data", "raw", "benchmarks")
    broad = pd.read_parquet(os.path.join(CH, "broad", "broad_edges.parquet"))
    bn = set(pd.read_parquet(os.path.join(CH, "broad", "broad_nodes.parquet")).disease_id)
    classes = pd.read_csv(os.path.join(CH, "specific", "mechanism_node_classes.tsv"), sep="\t")
    fda = pd.read_csv(os.path.join(bm, "fda_repurposing", "fda_repurposing_pairs.tsv"), sep="\t")
    pos = pd.read_parquet(os.path.join(bm, "ixidn", "ixidn_plus_pairs.parquet"))
    neg = pd.read_parquet(os.path.join(bm, "ixidn", "ixidn_minus_pairs.parquet"))
    spec, kept, _ = build.specific_channel(rebuilt, bn)
    un, _ = build.union_graph(broad, spec, fda)
    arm, _, _ = build.arm_weights(un, kept, classes)
    tr, _ = build.strip_benchmark(arm, pos, neg)
    build.check_against(arm, tr, os.path.join(Z, "graph"))
