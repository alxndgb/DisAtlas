"""Supplementary Figure SI_degree_distributions: the power-law fit table is rebuilt from the released per-disease degrees
(disatlas.rebuild.powerlaw_csn). Needs the optional extra requirements-cluster.txt (powerlaw); skipped otherwise.
Observed fit and likelihood-ratio columns are always checked (a few minutes). The bootstrap columns are checked against
a completed sharded run: set DISATLAS_POWERLAW_BOOT=<OUTDIR of observed + boot_task array + collect>."""
import os, numpy as np, pandas as pd, pytest
pytest.importorskip("powerlaw")
from disatlas.rebuild import powerlaw_csn as P

REF = pd.read_csv(P.DD / "powerlaw_summary.tsv", sep="\t")
OBS = ["alpha", "alpha_mle", "sigma", "xmin", "D"] + [f"LR_{a}_{s}" for a in P.ALTS for s in ("R", "p")]
INT = ["n_degree0", "n_diseases", "n_fit", "n_tail"]


def _cmp(S, cols):
    for c in cols: np.testing.assert_allclose(S[c].to_numpy(float), REF[c].to_numpy(float), rtol=6e-6, atol=0, err_msg=c)


def test_observed(tmp_path):
    o = P.observed(tmp_path); S = pd.DataFrame([o[f"{g}:{m}"] for g, m in P.CONFIGS]); S["alpha_mle"] = S.alpha
    assert (S[["graph", "fit"]].to_numpy() == REF[["graph", "fit"]].to_numpy()).all()
    for c in INT: assert (S[c].to_numpy() == REF[c].to_numpy()).all(), c
    _cmp(S, OBS)


@pytest.mark.skipif(not os.environ.get("DISATLAS_POWERLAW_BOOT"), reason="set DISATLAS_POWERLAW_BOOT to a sharded run")
def test_bootstrap():
    S = pd.read_csv(os.path.join(os.environ["DISATLAS_POWERLAW_BOOT"], "powerlaw_summary.tsv"), sep="\t")
    assert list(S.columns) == list(REF.columns)
    for c in ["n_boot", "n_boot_D_ge_obs"] + INT: assert (S[c].to_numpy() == REF[c].to_numpy()).all(), c
    assert (S.best_rival == REF.best_rival).all()
    _cmp(S, OBS + ["p_gof", "p_gof_mc_se", "boot_xmin_median", "boot_alpha_median", "best_rival_R", "best_rival_p"])
