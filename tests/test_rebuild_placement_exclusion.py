"""Parity: disatlas.rebuild.placement_exclusion reproduces data/derived/placement_confidence/exclusion_sensitivity.tsv
(1,000 bootstrap draws and 500 null removals per rule and benchmark; set DISATLAS_SLOW=1, about 10 minutes)."""
import os, numpy as np, pandas as pd, pytest
from disatlas.rebuild import placement_exclusion as M

pytestmark = pytest.mark.skipif(not os.environ.get("DISATLAS_SLOW"), reason="slow; set DISATLAS_SLOW=1")


def test_exclusion_sensitivity():
    ref = pd.read_csv(M.PC / "exclusion_sensitivity.tsv", sep="\t"); new = M.compute()
    assert list(new.columns) == list(ref.columns) and len(new) == 20
    cat = ["benchmark", "exclusion"]; assert new[cat].equals(ref[cat])
    num = [c for c in ref.columns if c not in cat]
    np.testing.assert_allclose(new[num].to_numpy(float), ref[num].to_numpy(float), rtol=0, atol=1e-12, equal_nan=True)
