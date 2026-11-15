"""P_full training graph: released specific edges + Habcd broad graph -> exact rebuild (counts, edge set, w_arm bits).

The rebuild from the released inputs (licence-restricted disease-gene pairs not redistributed) is checked against
tests/pfull_rebuild_signature.json; the graph of record of the released model (tests/pfull_graph_signature.json)
is checked against the shipped model configuration.

The broad input is data/zenodo/graph/graph_edges.parquet (warm; skipped without the bundle). The same broad edge set
is rebuilt from data/raw by the opt-in test in tests/test_raw_graph.py, which also checks it against the P_full input
signature, so both steps together are a rebuild from released files.
"""
import hashlib
import json

import pandas as pd
import pytest

from disatlas import paths
from disatlas.graph import pfull_build as P

REC = json.load(open(paths.ROOT / "tests" / "pfull_graph_signature.json"))
SIG = json.load(open(paths.ROOT / "tests" / "pfull_rebuild_signature.json"))["warm"]


def test_specific_edges_file():
    assert hashlib.sha256(P.SPEC.read_bytes()).hexdigest() == open(str(P.SPEC) + ".sha256").read().split()[0] == REC["specific_edges_sha256"]
    s = pd.read_parquet(P.SPEC)
    assert list(s.columns) == ["u", "v", "w"] and len(s) == SIG["inputs"]["spec_edges"] == REC["inputs"]["spec_edges"]
    assert (s.u < s.v).all() and not s.duplicated(["u", "v"]).any() and (s.w > 0).all()


def test_tier1_input():
    t = P.tier1_table()
    assert len(t) == SIG["inputs"]["tier1_rows"] and P.set_sha(t) == P.EXPECT["tier1_set"]


@pytest.fixture(scope="module")
def rebuilt():
    if not P.BROAD.is_file():
        pytest.skip("data/zenodo/graph/graph_edges.parquet not present")
    return P.build(check=True)


def test_counts_and_signatures(rebuilt):
    G, st = rebuilt
    for k in ("union_edges", "fda_edges_removed", "fda_spec_edges_removed", "edges", "nodes", "broad", "spec", "both",
              "ga_added", "ga_reweighted"):
        assert st[k] == SIG[k], k
    assert st["inputs"] == SIG["inputs"]
    assert st["edgeset_sha256"] == SIG["edgeset_sha256"] == P.EXPECT["out_edgeset"][st["inputs"]["broad_edgeset_sha256"]]


def test_w_arm_bit_identical(rebuilt):
    G, st = rebuilt
    assert list(G.columns) == P.COLS
    assert st["w_arm_sha256"] == SIG["w_arm_sha256"]
    assert hashlib.sha256(pd.util.hash_pandas_object(G[["u", "v"]], index=False).values.tobytes()).hexdigest() == SIG["row_order_sha256"]


def test_record_signature():
    assert (REC["inputs"]["broad_edgeset_sha256"], REC["inputs"]["tier1_set_sha256"], REC["edgeset_sha256"]) == \
        (P.RECORD["broad_edgeset"], P.RECORD["tier1_set"], P.RECORD["out_edgeset"])
    assert REC["inputs"]["broad_edgeset_sha256"] == SIG["inputs"]["broad_edgeset_sha256"]
    cfg = paths.ZENODO / "model" / "disatlas_v3.2_ensemble" / "model_config.json"
    if not cfg.is_file():
        pytest.skip("data/zenodo/model not present")
    f = json.load(open(cfg))["final_model_training_graph"]
    assert (f["n_edges"], f["n_nodes"], f["edgeset_sha256"], f["w_arm_sha256"]) == \
        (REC["edges"], REC["nodes"], REC["edgeset_sha256"], REC["w_arm_sha256"])
