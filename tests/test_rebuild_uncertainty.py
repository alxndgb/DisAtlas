"""Parity: disatlas.rebuild.uncertainty reproduces data/derived/uncertainty/{seed_vs_sample_sd,headline_cosine_ci}.tsv
(B = 2,000; about two minutes)."""
import numpy as np, pandas as pd, pytest
from disatlas.rebuild import uncertainty as M

DV = M.ROOT / "data/derived/uncertainty"


@pytest.fixture(scope="module")
def built():
    return M.compute()


def test_seed_vs_sample(built):
    ref = pd.read_csv(DV / "seed_vs_sample_sd.tsv", sep="\t"); new = built[0]
    assert new[["stratum", "metric"]].equals(ref[["stratum", "metric"]])
    np.testing.assert_allclose(new[["seed_sd", "sample_sd_cluster"]].to_numpy(float), ref[["seed_sd", "sample_sd_cluster"]].to_numpy(float), rtol=0, atol=1e-15)


def test_headline_cosine_ci(built):
    ref = pd.read_csv(DV / "headline_cosine_ci.tsv", sep="\t"); new = built[1]
    assert new[["metric", "stratum", "n", "n_pos", "B"]].astype(str).equals(ref[["metric", "stratum", "n", "n_pos", "B"]].astype(str))
    np.testing.assert_allclose(new[["value", "lo", "hi"]].to_numpy(float), ref[["value", "lo", "hi"]].to_numpy(float), rtol=0, atol=1e-12)
