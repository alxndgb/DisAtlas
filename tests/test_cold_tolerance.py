"""Cold mode: the retrained ensemble (build/model, build/scores) agrees with the frozen record (data/zenodo).
The retrain repeats the frozen recipe with the same random streams; only floating-point rounding differs (thread
count, BLAS build), so scores agree closely but not bit for bit. Observed values come from cold retrains of
seeds 4707-4711 on the cluster with 32, 8 and 3 torch threads (README, "Cold run")."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from disatlas import paths

pytestmark = pytest.mark.skipif(paths.MODE != "cold", reason="cold mode only (DISATLAS_DATA_MODE=cold)")
SC = paths.BUILD / "scores"
# observed max |cold - frozen| for the 5-seed ensemble over retrains with 32, 8 and 3 torch threads: <= 0.0004,
# fda_dev 0.0008
TOL_AUROC = dict.fromkeys(["ixidn_matched", "ixidn_matched|test", "ixidn_matched|rare-rare",
                           "ixidn_matched|rare-common", "ixidn_pm", "fda_dev"], 0.005)
MIN_RHO_MATCHED = 0.99  # Spearman, ensemble scores of the 4,364 IxIDN-Matched pairs (observed >= 0.99996)
MIN_RHO_ALL = 0.99      # Spearman, ensemble scores of 10^6 fixed random disease pairs (observed >= 0.9997)


@pytest.fixture(scope="module")
def frozen(tmp_path_factory):
    from disatlas.rescore import score
    return score(paths.ZENODO / "model" / paths.MODEL_NAME, tmp_path_factory.mktemp("frozen"))


@pytest.fixture(scope="module")
def cold():
    return json.load(open(SC / "benchmark_auroc.json"))


@pytest.mark.parametrize("key", sorted(TOL_AUROC))
def test_benchmark_auroc_within_tolerance(frozen, cold, key):
    assert abs(cold[key] - frozen[key]) <= TOL_AUROC[key], (key, cold[key], frozen[key])


def test_matched_pair_scores_rank_correlated():
    a = pd.read_csv(SC / "ixidn_matched_disatlas_full_scores.tsv", sep="\t")
    b = pd.read_csv(paths.ZENODO / "scores" / "ixidn_matched_disatlas_full_scores.tsv", sep="\t")
    assert (a.pair_id.values == b.pair_id.values).all()
    assert spearmanr(a.p_ensemble, b.p_ensemble, nan_policy="omit")[0] >= MIN_RHO_MATCHED


def test_all_pairs_rank_correlated():
    a, b = (np.load(d / "disatlas_all_pairs_ensemble.npz") for d in (SC, paths.ZENODO / "scores"))
    assert (a["ids"] == b["ids"]).all()
    k = np.random.default_rng(0).choice(len(b["p_upper"]), 1_000_000, replace=False)
    assert spearmanr(a["p_upper"][k], b["p_upper"][k])[0] >= MIN_RHO_ALL
