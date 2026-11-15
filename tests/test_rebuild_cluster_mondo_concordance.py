"""Figure S16 c: agreement of the K = 55 clusters with the MONDO body-system classes, rebuilt from the staged MONDO
extract and the released cluster labels (disatlas.rebuild.cluster_mondo_concordance), against the m3_external block
of data/derived/si_cluster/sifu_cluster_13021_values.json.

Every field is compared exactly except the three within-cluster pair-share estimates (NOT_BIT_REPRODUCIBLE): their
sampler orders diseases with numpy's non-stable argsort, so the pairs drawn depend on the order of tied cluster labels,
which numpy does not specify. Rebuilt 0.62703 vs shipped 0.62678 (null mean 0.25375 vs 0.25370). The shipped values
are kept; the text prints none of them (only the permutation p, compared exactly here) and Figure S16 c reads the
shipped json. For these three the test checks only the relation the figure shows (observed above the null 97.5%).
"""
import json
from pathlib import Path
from disatlas.rebuild import cluster_mondo_concordance as M

S = json.load(open(Path(__file__).resolve().parents[1] / "data/derived/si_cluster/sifu_cluster_13021_values.json"))["m3_external"]
NOT_BIT_REPRODUCIBLE = {"pair_share_within_cluster", "pair_share_within_cluster_null_mean", "pair_share_within_cluster_null_95"}


def _eq(a, b):
    if isinstance(b, list):
        return len(a) == len(b) and all(_eq(x, y) for x, y in zip(a, b))
    if isinstance(b, float):
        return abs(a - b) <= 1e-12 * max(1.0, abs(b))
    return a == b


def test_m3_external():
    T = M.compute()
    assert set(T) == set(S)
    bad = {k: (T[k], S[k]) for k in S if k not in NOT_BIT_REPRODUCIBLE and not _eq(T[k], S[k])}
    assert not bad, bad
    for V in (T, S):
        assert V["pair_share_within_cluster"] > V["pair_share_within_cluster_null_95"][1] > V["pair_share_random"]
