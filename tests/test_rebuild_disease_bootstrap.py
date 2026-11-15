"""Parity: disatlas.rebuild.disease_bootstrap reproduces data/derived/ixidn_robustness/subgroups_shared_bootstrap_B2000.tsv
(B = 2,000, four scopes x five systems; counts and p-values exact, decimals to 1e-12; a few minutes)."""
import numpy as np, pandas as pd, pytest
from disatlas.rebuild import disease_bootstrap as M


@pytest.fixture(scope="module")
def built():
    return M.compute()


def test_subgroups_shared_bootstrap(built):
    ref = pd.read_csv(M.ROOT / "data/derived/ixidn_robustness/subgroups_shared_bootstrap_B2000.tsv", sep="\t")
    assert list(built.columns) == list(ref.columns) and len(built) == len(ref) == 20
    assert built[["scope", "system", "n_pairs", "n_pos", "B", "n_valid"]].equals(ref[["scope", "system", "n_pairs", "n_pos", "B", "n_valid"]])
    np.testing.assert_array_equal(built.p_two.to_numpy(float), ref.p_two.to_numpy(float))   # p-values exact
    num = ["auroc", "lo", "hi", "dis_minus_sys", "dms_lo", "dms_hi"]
    np.testing.assert_allclose(built[num].to_numpy(float), ref[num].to_numpy(float), rtol=0, atol=1e-12, equal_nan=True)
