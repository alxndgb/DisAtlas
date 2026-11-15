"""Cold mode: the training graph rebuilt into build/graph equals the release rebuild of P_full edge for edge
(edge set, row order and walk weights; signature in tests/pfull_rebuild_signature.json: "warm" with the shipped
broad graph, "cold" with DISATLAS_RAWGRAPH=1)."""
import hashlib
import json

import pandas as pd
import pytest

from disatlas import paths
from disatlas.graph import pfull_build

pytestmark = pytest.mark.skipif(paths.MODE != "cold", reason="cold mode only (DISATLAS_DATA_MODE=cold)")
G = paths.BUILD / "graph"
REB = json.load(open(paths.ROOT / "tests" / "pfull_rebuild_signature.json"))


def test_cold_graph_identical():
    if not (G / "pfull_graph_edges.parquet").exists():
        pytest.skip(f"{G} not built (graph stage not run)")
    g = pd.read_parquet(G / "pfull_graph_edges.parquet")
    s = json.load(open(G / "pfull_graph_stats.json"))
    SIG = next(REB[k] for k in ("warm", "cold") if REB[k]["inputs"] == s["inputs"])
    assert (s["edges"], s["nodes"]) == (len(g), pd.concat([g.u, g.v]).nunique()) == (SIG["edges"], SIG["nodes"])
    assert pfull_build.edgeset_sha256(g.u, g.v) == pfull_build.EXPECT["out_edgeset"][s["inputs"]["broad_edgeset_sha256"]] == SIG["edgeset_sha256"]
    assert pfull_build.w_arm_sha256(g) == SIG["w_arm_sha256"]
    assert hashlib.sha256(pd.util.hash_pandas_object(g[["u", "v"]], index=False).values.tobytes()).hexdigest() \
        == SIG["row_order_sha256"]
