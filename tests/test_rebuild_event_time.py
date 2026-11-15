"""Parity: disatlas.rebuild.event_time reproduces the stored event-time tables from the staged per-event percentile table
(data/zenodo/event_time/, 41 MB; skipped when absent).

  * point estimates (all-disease and test-disease pools) and the merge / matching census: always run, counts exact, decimals 1e-9;
  * deltas with drug-clustered bootstrap intervals and the pre-specified decision (B = 2,000 per pool = 10 shards x 200): run when
    DISATLAS_ETBOOT points to the directory written by `python -m disatlas.rebuild.event_time boot POOL SHARD 200 DIR` for
    POOL in {test, all}, SHARD 0..9 (about 20 single-core jobs); otherwise skipped. Counts exact, decimals 1e-9."""
import os, json
import numpy as np, pandas as pd, pytest
from disatlas.rebuild import event_time as M

SI = M.ROOT / "data/derived/si_event_time"; F7 = M.ROOT / "data/derived/fig7_prospective"
pytestmark = pytest.mark.skipif(not (M.ZEN / "event_percentiles.parquet").exists(), reason="staged event-time table not present")
BOOT = os.environ.get("DISATLAS_ETBOOT")


@pytest.fixture(scope="module")
def point():
    return M.summarize(None)


def _cmp(new, ref, keys, num, atol=1e-9):
    assert list(new.columns) == list(ref.columns), (list(new.columns), list(ref.columns))
    assert len(new) == len(ref)
    assert new[keys].astype(str).reset_index(drop=True).equals(ref[keys].astype(str).reset_index(drop=True))
    np.testing.assert_allclose(new[num].to_numpy(float), ref[num].to_numpy(float), rtol=0, atol=atol, equal_nan=True)


@pytest.mark.parametrize("pool", ["all", "test"])
def test_point_estimates(point, pool):
    ref = pd.read_csv(SI / f"point_estimates_{pool}.tsv", sep="\t"); new = point["point"][pool]
    _cmp(new, ref, ["pool", "scope", "analysis", "method"], ["E1", "E2"])


def test_census(point):
    ref = pd.read_csv(SI / "merge_and_matching_census.tsv", sep="\t"); new = point["census"]
    assert list(new.columns) == list(ref.columns) and len(new) == len(ref)
    assert new.astype(str).reset_index(drop=True).equals(ref.astype(str).reset_index(drop=True))


@pytest.mark.skipif(BOOT is None, reason="set DISATLAS_ETBOOT to the sharded B = 2000 bootstrap directory")
def test_deltas_and_decision():
    S = M.summarize(BOOT); ref = pd.read_csv(SI / "deltas.tsv", sep="\t"); new = S["deltas"]
    num = [c for c in ref.columns if c not in ("pool", "scope", "analysis", "comparator")]
    new["scope"] = new.scope.astype(str); ref["scope"] = ref.scope.astype(str)
    key = ["pool", "scope", "analysis", "comparator"]
    a = new.sort_values(key).reset_index(drop=True); b = ref.sort_values(key).reset_index(drop=True)
    _cmp(a, b, key, num)
    d = M.decision(S["deltas"]); st = json.load(open(F7 / "decision_event_time.json"))
    assert d["verdict"] == st["verdict"] and d["n_pass"] == st["n_pass"] and d["passes"] == st["passes"] and d["rule"] == st["rule"]
    for k, v in st["ci"].items(): np.testing.assert_allclose(d["ci"][k], v, rtol=0, atol=1e-9)


def test_coverage():
    ref = pd.read_csv(SI / "coverage.tsv", sep="\t"); new = M.coverage()
    assert list(new.columns) == list(ref.columns) and len(new) == len(ref)
    assert new.sort_values(["pool", "method"]).astype(str).reset_index(drop=True).equals(ref.sort_values(["pool", "method"]).astype(str).reset_index(drop=True))


@pytest.mark.skipif(BOOT is None, reason="set DISATLAS_ETBOOT to the sharded B = 2000 bootstrap directory")
def test_top10_precision_table():
    """Class (a). Table S13 (pooled top-10% precision) rebuilt from the released per-drug precision rows and the B = 2,000
    drug-clustered bootstrap: every cell (DisAtlas and comparator columns, intervals at 4 d.p.) equals the shipped table."""
    new = M.top10_precision_table(M.summarize(BOOT)); ref = pd.read_csv(SI / "table_pooled_top10_precision.tsv", sep="\t")
    assert list(new.columns) == list(ref.columns) and len(new) == len(ref)
    assert new.astype(str).equals(ref.astype(str))


@pytest.mark.skipif(BOOT is None, reason="set DISATLAS_ETBOOT to the sharded B = 2000 bootstrap directory")
def test_abs_percentiles():
    ref = pd.read_csv(F7 / "abs_percentiles.tsv", sep="\t"); ref = ref[ref.method != "chance_expected"].copy(); ref["scope"] = ref.scope.astype(str)
    new = M.abs_percentiles(BOOT); key = ["pool", "scope", "analysis", "method"]
    a = new.sort_values(key).reset_index(drop=True); b = ref[key + ["E1", "lo", "hi", "n_rep"]].sort_values(key).reset_index(drop=True)
    assert len(a) == len(b) == 312 and a[key + ["n_rep"]].astype(str).equals(b[key + ["n_rep"]].astype(str))
    np.testing.assert_allclose(a[["E1", "lo", "hi"]].to_numpy(float), b[["E1", "lo", "hi"]].to_numpy(float), rtol=0, atol=1e-9)
