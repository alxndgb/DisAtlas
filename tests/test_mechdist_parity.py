"""Parity of the MechDist rebuild (code/03_network) with the production network (Methods, mechanism
distance; Results, MechDist network: 13,021 diseases, 332,082 edges).

Run the build first:
    python code/03_network/mechdist_distances.py && python code/03_network/build_mechdist_network.py
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT = REPO_ROOT / "data" / "03_rebuilt"
REF = OUT / "reference"
S5 = OUT / "work" / "dis_dist_all_df_with_crank.parquet"

pytestmark = pytest.mark.skipif(not S5.exists(), reason="run code/03_network/mechdist_distances.py first")


@pytest.fixture(scope="module")
def ref():
    return json.load(open(REF / "mechdist_parity_reference.json"))


def s5_digest(df):
    h = hashlib.sha256()
    h.update("\n".join(df["pair"]).encode())
    for c in ["dist_GO", "dist_HPO", "px_HPA", "px_DOM", "crank"]:
        h.update(np.ascontiguousarray(df[c].to_numpy("float64")).tobytes())
    return h.hexdigest()


def test_s3_thresholds(ref):
    got = json.load(open(OUT / "S3_thresholds.json"))
    for ch in ("go", "hpo"):
        assert got[ch]["n"] == ref["s3"][ch]["n"]
        assert got[ch]["mean"] == pytest.approx(ref["s3"][ch]["mean"], rel=0, abs=1e-15)
        assert got[ch]["frac_le_cut"] == ref["s3"][ch]["frac_le_cut"]


def test_s5_frame_bit_identical(ref):
    df = pd.read_parquet(S5, columns=["pair", "dist_GO", "dist_HPO", "px_HPA", "px_DOM", "crank"])
    assert len(df) == ref["s5_rows"]
    assert s5_digest(df) == ref["s5_sha256"]


def test_s5_budget(ref):
    assert json.load(open(OUT / "S5_px_budget.json"))["quotas"] == ref["quotas"]


def test_network_edges_identical(ref):
    e = pd.read_parquet(OUT / "mechdist_network_edges.parquet")
    r = pd.read_parquet(REF / "mechdist_network_edges_reference.parquet")
    assert len(e) == ref["n_edges"]
    key = ["disease_a", "disease_b"]
    m = e.merge(r, on=key, how="outer", suffixes=("", "_ref"), indicator=True)
    assert (m["_merge"] == "both").all()
    assert (m["channel_class"] == m["channel_class_ref"]).all()
    assert (m["origin"] == m["origin_ref"]).all()


def test_network_nodes_identical(ref):
    n = pd.read_parquet(OUT / "mechdist_network_nodes.parquet")
    r = pd.read_parquet(REF / "mechdist_network_nodes_reference.parquet")
    assert len(n) == ref["n_nodes"]
    pd.testing.assert_frame_equal(n.reset_index(drop=True), r.reset_index(drop=True))
