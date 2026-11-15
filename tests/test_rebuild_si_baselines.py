"""Supplementary Table S5 (si_baselines): every AUROC, interval and paired difference of the shipped table is rebuilt
from DisAtlas cosines and the staged training-free pair scores (disatlas.rebuild.si_baselines; B = 2,000, about 15 min).
Shipped values carry 6 decimals; the rebuild is rounded the same way and must then agree exactly. Slow: DISATLAS_SLOW=1."""
import os, numpy as np, pandas as pd, pytest
from disatlas.rebuild import si_baselines as SB

REF = pd.read_csv(SB.D / "si_baselines/baselines_auroc_B2000.tsv", sep="\t")


def test_pair_inputs():
    S = SB.pairs()
    assert S[SB.IG].notna().all().all() and S.loc[S.nd_c.notna(), [f"ho_{c}" for c in SB.HO]].notna().all().all()


@pytest.mark.skipif(not os.environ.get("DISATLAS_SLOW"), reason="B = 2,000 bootstrap; set DISATLAS_SLOW=1")
def test_table():
    T = SB.build()
    assert T.shape == REF.shape and (T[["block", "subset", "score"]].to_numpy() == REF[["block", "subset", "score"]].to_numpy()).all()
    for c in ["auc", "lo", "hi", "d", "dlo", "dhi"]:
        np.testing.assert_array_equal(np.round(T[c].to_numpy(float), 6), REF[c].to_numpy(float), err_msg=c)
