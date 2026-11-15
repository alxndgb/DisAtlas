"""Figure 6 / Table 1 inputs: data/derived/fig6_comparison/{per_run_scores,node_disjoint_per_seed}.tsv are rebuilt from
released files by disatlas.rebuild.fig6_inputs (same rows, same columns, values equal at the stored precision)."""
from pathlib import Path
import numpy as np, pandas as pd
from disatlas.rebuild import fig6_inputs as F

D = Path(__file__).resolve().parents[1] / "data/derived/fig6_comparison"


def _same(new, old, atol):
    assert list(new.columns) == list(old.columns)
    new = new.set_index("blind_id").loc[old.blind_id].reset_index()
    for c in old.columns:
        if not pd.api.types.is_numeric_dtype(old[c]):
            assert (new[c].astype(str) == old[c].astype(str)).all(), c
        else:
            a, b = new[c].to_numpy(float), old[c].to_numpy(float)
            assert (np.isnan(a) == np.isnan(b)).all(), c
            m = ~np.isnan(b); assert np.abs(a[m] - b[m]).max() <= atol, (c, np.abs(a[m] - b[m]).max())


def test_per_run_scores():
    _same(F.per_run_scores(), pd.read_csv(D / "per_run_scores.tsv", sep="\t"), 5e-7)


def test_node_disjoint_per_seed():
    _same(F.node_disjoint_per_seed(), pd.read_csv(D / "node_disjoint_per_seed.tsv", sep="\t"), 5e-7)
