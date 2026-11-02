"""Broad channel rebuild parity: warm (shipped pair distances) and cold (dist_HPO recomputed from hp.obo).

Warm needs data/zenodo/channels_raw/broad (skipped otherwise); cold also needs DISATLAS_COLD=1 (slow, ~10 min).
"""
import os

import pandas as pd
import pytest

from disatlas.channels import broad

has_pairs = pytest.mark.skipif(not os.path.isdir(broad.ZEN), reason="Zenodo channels_raw/broad not present")
cold_opt_in = pytest.mark.skipif(os.environ.get("DISATLAS_COLD") != "1", reason="set DISATLAS_COLD=1 to run")
QUOTAS = {"hpo_only": 39_726, "hpo+px": 160_274, "all": 200_000, "go+px": 200_000, "go+hpo": 200_000, "go_only": 200_000}


def _check(cold):
    e, q = broad.build(cold)
    ref = pd.read_parquet("data/raw/channels/broad/broad_edges.parquet")
    assert len(e) == 332_082 and len(set(e.u) | set(e.v)) == 13_021
    assert {k: q[k] for k in QUOTAS} == QUOTAS
    pd.testing.assert_frame_equal(e, ref.sort_values(["u", "v"]).reset_index(drop=True))


@has_pairs
def test_broad_warm():
    _check(cold=False)


@has_pairs
@cold_opt_in
def test_broad_cold():
    _check(cold=True)


def test_fda_controls_removed_from_broad():
    from disatlas.graph.build import pair_key
    b = pd.read_parquet("data/raw/channels/broad/broad_edges.parquet")
    f = pd.read_csv("data/raw/benchmarks/fda_repurposing/fda_repurposing_pairs.tsv", sep="\t")
    hit = pair_key(b.u, b.v).isin(set(pair_key(f.disease_a, f.disease_b))).values
    kept = b[~hit]
    assert int(hit.sum()) == 848 and len(kept) == 331_234
    assert len(set(b.u) | set(b.v)) == 13_021 and len(set(kept.u) | set(kept.v)) == 13_020


def test_crank_ranks_weakest_pair_last():
    import numpy as np
    r = broad.crank(np.array([[0.9, 0.9], [0.9, 0.9], [0.1, 0.2]]))
    assert r[2] == r.max()
