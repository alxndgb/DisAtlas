"""Cold mode: channels and training graph rebuilt into build/ equal the shipped ones edge for edge."""
import json

import pandas as pd
import pytest

from disatlas import paths

pytestmark = pytest.mark.skipif(paths.MODE != "cold", reason="cold mode only (DISATLAS_DATA_MODE=cold)")
CH, G = paths.BUILD / "channels", paths.BUILD / "graph"
OUTPUTS = {"broad": ["broad_edges.parquet", "broad_nodes.parquet"], "specific": ["disease_module_links.parquet"]}


def read(p):
    return pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p, sep="\t")


@pytest.mark.parametrize("track", ["broad", "specific"])
def test_channel_identical(track):
    if not (CH / track).is_dir():
        pytest.skip(f"{CH / track} not built (channel stage not run)")
    for name in OUTPUTS[track]:  # the channel outputs; the other files in data/raw/channels are their inputs
        ref = paths.RAW / "channels" / track / name
        r, n = read(ref), read(CH / track / ref.name)
        assert list(r.columns) == list(n.columns), ref.name
        key = list(r.columns)
        r, n = r.sort_values(key).reset_index(drop=True), n.sort_values(key).reset_index(drop=True)
        pd.testing.assert_frame_equal(r, n, check_dtype=False, check_exact=False, rtol=0, atol=1e-12)


def test_graph_counts_and_edges():
    from disatlas.graph.build import check_against
    s = json.load(open(G / "graph_build_stats.json"))
    assert (s["broad"]["nodes"], s["broad"]["edges"]) == (13021, 332082)
    assert (s["specific"]["links"], s["specific"]["edges"]) == (422804, 191019)
    assert (s["union"]["edges"], s["strip"]["train_edges"]) == (498155, 488103)
    ref = json.load(open(paths.ZENODO / "graph" / "graph_build_stats.json"))
    assert s["train_edgeset_sha256"] == ref["train_edgeset_sha256"]
    check_against(pd.read_parquet(G / "union_graph_edges.parquet"),
                  pd.read_parquet(G / "training_graph_edges.parquet"), str(paths.ZENODO / "graph"))
