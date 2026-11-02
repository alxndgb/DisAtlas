"""IxIDN-Matched QC: the benchmark is rebuilt from the stored language-model verdicts and the frozen matching
table, and the verdict counts equal the Methods / Figure 3b values. COLD QC (opt-in) is checked at a tolerance."""
import json
import os

import pandas as pd
import pytest

from disatlas.benchmarks import ixidn_matched as im

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
BM = os.path.join(ROOT, "data", "raw", "benchmarks", "ixidn_matched")
VERD = os.path.join(BM, "llm_curation_verdicts.tsv")
EDGES = os.path.join(ROOT, "data", "zenodo", "figures_main", "ixidn_graph_edges.parquet")
COLD = os.path.join(ROOT, "build", "qc", "llm_curation_verdicts.tsv")


@pytest.fixture(scope="module")
def built():
    return im.build(VERD)


def test_warm_build_equals_shipped_benchmark(built):
    M, _ = built
    ref = pd.read_csv(os.path.join(BM, "ixidn_matched_pairs.tsv"), sep="\t")
    pd.testing.assert_frame_equal(M.sort_values("pair_id").reset_index(drop=True),
                                  ref.sort_values("pair_id").reset_index(drop=True), check_dtype=False)
    assert (len(M), int(M.y.sum()), int((M.split == "test").sum())) == (4364, 2182, 1094)


def test_verdict_summary_matches_manuscript(built):
    _, S = built
    pd.testing.assert_frame_equal(S.fillna("").astype(str),
                                  pd.read_csv(im.SUMMARY, sep="\t", dtype=str, keep_default_na=False))
    n = dict(zip(S.verdict + "/" + S.decision, S.pairs))
    assert n["STRONG/kept"] == 1285 and n["PLAUSIBLE/kept"] == 897 and n["all/kept"] == 2182
    assert n["all/kept"] + n["all/removed"] == 4302
    rep = S[S.verdict == "repeat check"].iloc[0]
    assert rep.pairs == 944 and round(100 * rep.share, 1) == 84.7


@pytest.mark.skipif(not os.path.exists(EDGES), reason="Zenodo bundle not present")
def test_not_in_graph_pairs_have_no_shared_drug():
    V = pd.read_csv(VERD, sep="\t")
    E = pd.read_parquet(EDGES)
    has = {frozenset(k) for k, s in zip(zip(E.disease_1, E.disease_2), E.shared_drugs) if json.loads(s)}
    none = [frozenset(k) not in has for k in zip(V.disease_a, V.disease_b)]
    assert list(V.verdict[none].unique()) == ["NOT_IN_GRAPH"] and sum(none) == 57


@pytest.mark.skipif(not os.path.exists(COLD), reason="COLD QC not run (./run_pipeline.sh --cold-qc)")
def test_cold_qc_within_tolerance():
    """A rerun uses a reconstructed prompt and a model of the user's choice, so verdicts are not identical.
    Tolerance: keep/drop agreement with the stored verdicts >= 0.75 (the stored repeat check gave 0.847) and the
    rebuilt benchmark keeps >= 75% of the stored positives."""
    k = lambda d: {frozenset(x): v for x, v in zip(zip(d.disease_a, d.disease_b), d.verdict)}
    a, b = k(pd.read_csv(COLD, sep="\t")), k(pd.read_csv(VERD, sep="\t"))
    assert a.keys() == b.keys()
    agree = sum((a[p] in im.KEEP) == (b[p] in im.KEEP) for p in b) / len(b)
    assert agree >= 0.75, agree
    M, _ = im.build(COLD)
    ref = pd.read_csv(os.path.join(BM, "ixidn_matched_pairs.tsv"), sep="\t")
    kept = len(set(M.pair_id[M.y == 1]) & set(ref.pair_id[ref.y == 1])) / int(ref.y.sum())
    assert kept >= 0.75, kept
