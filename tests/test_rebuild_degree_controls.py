"""Parity: disatlas.rebuild.degree_controls vs data/derived/degree_controls/degree_controls_table.tsv (Table S8 (SI_degree_controls)).

Class (c), rebuilt but not bit-identical. Everything that does not involve the DisAtlas score reproduces (counts exactly, decimals to 1e-5: tie handling/summation order in the weighted AUC) (all counts,
the degree-matched subsets, degree alone, and every AUROC and interval of the five comparison systems). The two DisAtlas rows
(full model, benchmark-diseases-held-out arm) and the differences that contain the held-out arm were computed by the authors from
internal per-pair DisAtlas scores that are not in the release (Spearman 0.99 with the released cosines); the released per-pair
scores are the 5-seed mean cosine, which gives AUROCs 0.002-0.013 higher (CI ends up to 0.016). The shipped values are kept; the test pins the exact part and bounds the deviation.
"""
import numpy as np, pandas as pd, pytest
from disatlas.rebuild import degree_controls as M

EXACT_SYSTEMS = ["Claude Chat", "Claude Sci", "Phylo Biomni", "Gemini Pro 3.1", "TxGNN (published checkpoint)"]
COUNTS = ["n_pairs", "n_pos", "n_diseases", "dm_pairs", "dm_pos"]
DEGREE = ["dm_degree_alone_auroc", "dm_degree_alone_ci_lo", "dm_degree_alone_ci_hi"]
SCORED = ["dm_auroc", "dm_ci_lo", "dm_ci_hi", "within_auroc", "within_ci_lo", "within_ci_hi"]
DELTA = [c for c in pd.read_csv(M.REF, sep="\t").columns if "delta" in c]


@pytest.fixture(scope="module")
def tables():
    return M.compute(2000), pd.read_csv(M.REF, sep="\t")


def test_rows_counts_and_degree_alone_exact(tables):
    new, ref = tables
    assert new[["system", "case_set"]].equals(ref[["system", "case_set"]])
    assert (new[COUNTS].to_numpy() == ref[COUNTS].to_numpy()).all()
    np.testing.assert_allclose(new[DEGREE].to_numpy(float), ref[DEGREE].to_numpy(float), rtol=0, atol=1e-5)


def test_comparison_systems_exact(tables):
    new, ref = tables; m = ref.system.isin(EXACT_SYSTEMS).to_numpy()
    np.testing.assert_allclose(new.loc[m, SCORED].to_numpy(float), ref.loc[m, SCORED].to_numpy(float), rtol=0, atol=1e-5)


def test_disatlas_rows_deviation_bounded(tables):
    """(c): DisAtlas / held-out rows and the deltas use the retired decoder score in the stored table."""
    new, ref = tables; m = (~ref.system.isin(EXACT_SYSTEMS)).to_numpy()
    d = np.abs(new.loc[m, SCORED].to_numpy(float) - ref.loc[m, SCORED].to_numpy(float))
    assert d.max() < 0.02
    dd = np.abs(new[DELTA].to_numpy(float) - ref[DELTA].to_numpy(float))
    assert np.nanmax(dd) < 0.02
