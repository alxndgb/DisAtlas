"""Cold mode: the retrained ensemble (build/model, build/scores) agrees with the frozen record (data/zenodo)
statistically. The retrain uses different random streams from the frozen run, so scores are not bit-identical.
Tolerances come from a cold retrain of seeds 4707-4711 on the cluster (README, "Cold run")."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from disatlas import paths

pytestmark = pytest.mark.skipif(paths.MODE != "cold", reason="cold mode only (DISATLAS_DATA_MODE=cold)")
SC = paths.BUILD / "scores"
TOL_AUROC = {"ixidn_matched": 0.03, "ixidn_matched|test": 0.04, "ixidn_matched|rare-rare": 0.05,
             "ixidn_matched|rare-common": 0.05, "ixidn_pm": 0.03, "fda_dev": 0.08}
MIN_RHO_MATCHED = 0.6   # Spearman, ensemble scores of the 4,364 IxIDN-Matched pairs
MIN_RHO_ALL = 0.6       # Spearman, ensemble scores of 10^6 fixed random disease pairs


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
