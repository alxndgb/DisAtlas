"""Parity: disatlas.rebuild.degree_distributions.summary vs data/derived/degree_distributions/degree_summary.tsv
(per-graph degree summary of the released graph: counts exact, decimals at the stored 6 significant digits)."""
import numpy as np, pandas as pd
from disatlas.rebuild import degree_distributions as M
from disatlas import paths


def test_degree_summary(tmp_path):
    ref = pd.read_csv(paths.DERIVED / "degree_distributions/degree_summary.tsv", sep="\t")
    S = M.summary(M.load()); S.to_csv(tmp_path / "s.tsv", sep="\t", index=False, float_format="%.6g")
    new = pd.read_csv(tmp_path / "s.tsv", sep="\t")
    assert list(new.columns) == list(ref.columns) and len(new) == len(ref)
    for c in ref.columns:
        if pd.api.types.is_numeric_dtype(ref[c]): np.testing.assert_allclose(new[c].astype(float), ref[c].astype(float), rtol=2e-6, atol=0, err_msg=c)
        else: assert (new[c].astype(str) == ref[c].astype(str)).all(), c
