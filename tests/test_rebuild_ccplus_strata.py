"""Parity: disatlas.rebuild.ccplus_strata reproduces data/derived/si_rare_common/ccplus_strata.tsv.
Full B = 2000 run takes a few minutes; set DISATLAS_SLOW=1 to run it (otherwise point estimates and counts only)."""
import os
import numpy as np, pandas as pd, pytest
from disatlas.rebuild import ccplus_strata as M

REF = M.ROOT / "data/derived/si_rare_common/ccplus_strata.tsv"
KEY = ["benchmark", "scope", "stratum"]
INT = ["n_pairs", "n_pos", "n_neg", "n_added_neg", "n_diseases"]
PT = ["auroc", "shortcut_stack_matched", "delta_vs_stack"]
CI = ["ci_lo", "ci_hi", "delta_vs_stack_lo", "delta_vs_stack_hi", "delta_vs_stack_p_le0", "n_boot"]


@pytest.fixture(scope="module")
def tables():
    B = 2000 if os.environ.get("DISATLAS_SLOW") else 10
    return M.compute(B), pd.read_csv(REF, sep="\t"), B


def test_rows_and_counts(tables):
    new, ref, _ = tables
    assert new[KEY].equals(ref[KEY])
    assert (new[INT].to_numpy() == ref[INT].to_numpy()).all()
    assert new.note.fillna("").equals(ref.note.fillna(""))


def test_point_estimates(tables):
    new, ref, _ = tables
    np.testing.assert_allclose(new[PT].to_numpy(float), ref[PT].to_numpy(float), rtol=0, atol=1e-12, equal_nan=True)


def test_intervals(tables):
    new, ref, B = tables
    if B != 2000: pytest.skip("set DISATLAS_SLOW=1 for the B = 2000 intervals")
    np.testing.assert_allclose(new[CI].to_numpy(float), ref[CI].to_numpy(float), rtol=0, atol=1e-12, equal_nan=True)
    assert new.ci_bootstrap.fillna("").equals(ref.ci_bootstrap.fillna(""))
