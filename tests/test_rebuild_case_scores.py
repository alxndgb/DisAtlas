"""Parity: disatlas.rebuild.case_scores reproduces data/derived/si_explainability/case_scores_cosine.tsv exactly (needs the
released cosine matrix, 0.7 GB; skipped when absent)."""
import pandas as pd, pytest
from disatlas.rebuild import case_scores as M

pytestmark = pytest.mark.skipif(not (M.ROOT / "data/zenodo/model_outputs/cosine_matrix/cosine_ens5_13021.npy").exists(), reason="released cosine matrix not present")


def test_case_scores():
    ref = pd.read_csv(M.ROOT / "data/derived/si_explainability/case_scores_cosine.tsv", sep="\t"); new = M.compute()
    assert list(new.columns) == list(ref.columns)
    assert new.astype(str).equals(ref.astype(str)), (new, ref)
