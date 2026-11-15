"""Figure S16 b: within/between distances and silhouettes at K = 55 are rebuilt from the released cosine matrix and
cluster labels (disatlas.rebuild.cluster_within_between) and equal the shipped table and json block."""
import json
from pathlib import Path
import numpy as np, pandas as pd
import pytest
from disatlas.rebuild import cluster_within_between as W

D = Path(__file__).resolve().parents[1] / "data/derived/si_cluster"


@pytest.fixture(scope="module")
def built():
    return W.compute()


def test_table(built):
    T, _ = built
    S = pd.read_csv(D / "within_between_K55_13021.tsv.gz", sep="\t")
    assert list(T.columns) == list(S.columns) and (T.mondo_id.to_numpy() == S.mondo_id.to_numpy()).all()
    assert (T.cluster.to_numpy() == S.cluster.to_numpy()).all()
    for c in ("within", "nearest_other", "silhouette"):  # stored at 5 decimals
        assert np.abs(T[c].to_numpy() - S[c].to_numpy()).max() <= 5e-6 + 1e-12, c


def test_internal(built):
    _, I = built
    S = json.load(open(D / "sifu_cluster_13021_values.json"))["internal"]
    assert set(I) == set(S)
    for k, v in S.items():
        assert I[k] == v if isinstance(v, int) else abs(I[k] - v) <= 1e-9, (k, I[k], v)
