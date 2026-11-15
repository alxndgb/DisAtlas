"""Parity: disatlas.rebuild.cluster_themes against data/derived/si_cluster_themes (Supplementary Table S10, K = 55).

(a) exact: body-system, MONDO top-5 and mechanism-term enrichment tables, cluster size and direct neoplasm share (rebuilt from the
    released labels and the staged MONDO 2026-02-06 inputs).
(c) rebuilt but not bit-identical: GO BP top-5 table, gene_coverage and the HPO / gene / GO pair-overlap folds. The authors' disease-gene
    table holds 610 pairs that are not redistributed and lacks 4,047 edges the release carries (2,193 relaxed-tier, 1,854 tier-1),
    and the authors' HPO table differs from the released one (data/README.md, "Table S10 (cluster themes) and the disease-gene table");
    the same code run on those tables reproduces the stored gene-share fold exactly and the HPO/GO Jaccard folds to float32
    precision (3.4e-7; checked by the maintainer, not repeatable here).
    pairs_module_jaccard_fold and graph_within_frac_specific rest on restricted module support and are shipped, not rebuilt.
    What must hold: the theme assigned to every one of the 55 clusters is the same whether the stored or the rebuilt tables are used.
"""
import numpy as np, pandas as pd, pytest
from disatlas.rebuild import cluster_themes as M

D = M.ROOT / "data/derived/si_cluster_themes"
rd = lambda f: pd.read_parquet(D / f"{f}.parquet")


@pytest.fixture(scope="module")
def T():
    return M.compute()


def same(new, ref, atol=1e-9):
    assert list(new.columns) == list(ref.columns) and len(new) == len(ref)
    for col in ref.columns:
        if ref[col].dtype.kind in "fi":
            np.testing.assert_allclose(new[col].to_numpy(float), ref[col].to_numpy(float), rtol=0, atol=atol, err_msg=col)
        else:
            assert (new[col].to_numpy() == ref[col].to_numpy()).all(), col


@pytest.mark.parametrize("name", ["bodysystem_by_cluster_long", "mondo_terms_top5_by_cluster", "mechanism_terms_by_cluster"])
def test_mondo_tables_exact(T, name):
    same(T[name], rd(name))


def test_size_and_neoplasm_exact(T):
    a, ref = T["cluster_annotation_rebuilt"], rd("cluster_annotation_k55")
    assert (a.cluster == ref.cluster).all() and (a["size"] == ref["size"]).all()
    np.testing.assert_allclose(a.neoplasm_share_direct, ref.neoplasm_share_direct, rtol=0, atol=1e-12)


def classify_with(**over):
    t = dict(bodysystem_by_cluster_long=rd("bodysystem_by_cluster_long"), mechanism_terms_by_cluster=rd("mechanism_terms_by_cluster"),
             go_bp_top5_by_cluster=rd("go_bp_top5_by_cluster"), mondo_terms_top5_by_cluster=rd("mondo_terms_top5_by_cluster"))
    t.update(over)
    return M.classify(t["bodysystem_by_cluster_long"], t["mechanism_terms_by_cluster"], t["go_bp_top5_by_cluster"],
                      t["mondo_terms_top5_by_cluster"], rd("cluster_annotation_k55"))


def test_stored_themes_match_reference():
    ref = pd.read_csv(D / "theme_reference.tsv", sep="\t").set_index("cluster").theme_class
    assert (classify_with().loc[ref.index] == ref).all()


def test_themes_unchanged_with_rebuilt_tables(T):
    ref = pd.read_csv(D / "theme_reference.tsv", sep="\t").set_index("cluster").theme_class
    new = classify_with(**{k: T[k] for k in ["bodysystem_by_cluster_long", "mechanism_terms_by_cluster",
                                             "go_bp_top5_by_cluster", "mondo_terms_top5_by_cluster"]})
    assert (new.loc[ref.index] == ref).all()
    assert new.value_counts().to_dict() == {"anatomy": 21, "unresolved": 13, "process": 11, "mixed": 6, "neoplasm": 4}


def test_go_and_pair_columns_known_deviation(T):
    """Pins the documented (c) deviation so a change of inputs is noticed; values are not expected to match the stored ones."""
    go, ref = T["go_bp_top5_by_cluster"], rd("go_bp_top5_by_cluster")
    assert go.shape == ref.shape
    assert int((go[go["rank"] == 1].term.to_numpy() == ref[ref["rank"] == 1].term.to_numpy()).sum()) == 42
    a, r = T["cluster_annotation_rebuilt"], rd("cluster_annotation_k55")
    for col in M.STORED_PAIR_COLS + ["gene_coverage"]:
        assert int((np.abs(a[col].to_numpy(float) - r[col].to_numpy(float)) > 0.003).sum()) >= 50, col

