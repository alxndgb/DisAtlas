"""Supplementary Figure SI_rare_common neighbour coverage (disatlas.rebuild.rc_coverage). Slow (about 5 min): DISATLAS_SLOW=1.

NOT_BIT_REPRODUCIBLE (documented exception, same class as cluster_mondo_concordance): the authors ranked neighbours on
float32 embeddings; the rebuild uses float64. Neighbour order differs only at near-ties below float32 resolution. In the
observed map this changes one disease (MONDO:0013751: ranks 15/16 are a rare and a common neighbour whose five-seed
cosines differ by 5e-9; shipped rank 15, rebuild 16). Under label permutation other same-label near-ties become visible.
Every difference is therefore bounded by one disease: k-coverage fractions by 1 / n_class, other-class-in-10 fractions
by 1 / (10 n_class). The test asserts that bound, exact equality of all p-values, observed counts and meta, and that every
value agrees at the 3 decimals the text prints."""
import json, os, numpy as np, pandas as pd, pytest
from disatlas.rebuild import rc_coverage as RC

pytestmark = pytest.mark.skipif(not os.environ.get("DISATLAS_SLOW"), reason="1,000 x 3 permutations; set DISATLAS_SLOW=1")
TIE = {"MONDO:0013751"}


@pytest.fixture(scope="module")
def built():
    return RC.build()


def test_per_disease(built):
    PD, *_ = built; P = pd.read_csv(RC.D / "per_disease_first_other_class.tsv.gz", sep="\t")
    for c in ["disease_id", "is_rare", "n_other_class_in_10", "train_degree"]:
        assert (PD[c].astype(str).values == P[c].astype(str).values).all(), c
    bad = set(P.disease_id[PD.first_other_class_rank.values != P.first_other_class_rank.values])
    assert bad == TIE


def test_tables(built):
    PD, S, E, meta = built; n = {"rare": meta["n_rare"], "common": meta["n_common"]}
    S0 = pd.read_csv(RC.D / "coverage_by_k_summary.tsv", sep="\t"); E0 = pd.read_csv(RC.D / "enrichment_other_class_frac10.tsv", sep="\t")
    S = S[S0.columns]; E = E[E0.columns]
    assert (S[["model", "anchor_class", "k"]].values == S0[["model", "anchor_class", "k"]].values).all()
    np.testing.assert_allclose(S.p_lower, S0.p_lower, rtol=0, atol=1e-12); np.testing.assert_allclose(E.p_lower, E0.p_lower, rtol=0, atol=1e-12)
    np.testing.assert_allclose(E.observed, E0.observed, rtol=0, atol=1e-12)
    bS = S0.anchor_class.map(n).rdiv(1.0).values + 1e-12; bE = E0.anchor_class.map(n).rdiv(0.1).values + 1e-12
    for c in ["mean", "lo", "hi"]:
        assert (np.abs(S[c].values - S0[c].values) <= bS).all(), c
        assert (np.round(S[c].values, 3) == np.round(S0[c].values, 3)).all(), c
    for c in ["mean_frac_other_in_10", "lo", "hi"]:
        assert (np.abs(E[c].values - E0[c].values) <= bE).all(), c
        assert (np.round(E[c].values, 3) == np.round(E0[c].values, 3)).all(), c
    m0 = json.load(open(RC.D / "rc_coverage_meta.json"))
    assert {k: m0[k] for k in meta} == meta
