"""Supplementary Table S7 (SI_heuristic_vs_ai): disatlas.rebuild.si_heur_vs_ai against data/derived/si_heur_vs_ai/heuristic_minus_ai_B2000.tsv.
Slow (B = 2,000 disease-cluster bootstrap, 2 subsets x 9 scores): DISATLAS_SLOW=1.
Compares every row at the 6 decimals stored (AUROCs, differences, interval ends)."""
import os, numpy as np, pandas as pd, pytest
from disatlas.rebuild import si_heur_vs_ai as S

pytestmark = pytest.mark.skipif(not os.environ.get("DISATLAS_SLOW"), reason="2,000-draw bootstrap; set DISATLAS_SLOW=1")


@pytest.fixture(scope="module")
def built():
    return S.build()


def test_pairs():
    df = S.pair_table()
    assert len(df) == 4357 and (df.pool == "test").sum() == 1093


def test_table(built):
    T = pd.read_csv(S.D / "si_heur_vs_ai/heuristic_minus_ai_B2000.tsv", sep="\t")
    assert list(built.columns) == list(T.columns) and len(built) == len(T) == 58
    for c in ["subset", "heuristic", "system"]:
        assert (built[c].values == T[c].values).all(), c
    for c in ["heuristic_auroc", "system_auroc", "diff", "ci_lo", "ci_hi"]:
        np.testing.assert_allclose(built[c].astype(float).to_numpy(), T[c].astype(float).to_numpy(), rtol=0, atol=5.1e-7, equal_nan=True, err_msg=c)
