"""Parity tests for code/00_ontology against the shipped tables and the manuscript.

Run:  python tests/test_ontology_parity.py        (or: pytest tests/test_ontology_parity.py)
Set ONTOLOGY_PARITY_FULL=1 to include the 13.2M-pair HPO disease-distance check (~2 min).
The GO/HPO disease-distance checks compare against data/02_mechdist/ and
data/zenodo/network/ (download from Zenodo first); they are skipped when those files are absent.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "code" / "00_ontology"))
import build_hpo_baseline as HB  # noqa: E402
import disease_distance as DD  # noqa: E402
import ontology as O  # noqa: E402

RAW = REPO_ROOT / "data" / "00_raw"
D_MAX = {"go_bp_semi": 23, "go_mf": 18, "go_cc": 18, "hpo": 25}


def _canon(df):
    return df.sort_values(list(df.columns)).reset_index(drop=True)


def _skip(msg):
    try:
        import pytest
    except ImportError:          # plain-python runner: report and move on
        print("SKIP", msg)
        return
    pytest.skip(msg)


def test_ancestor_tables_rebuild_from_obo():
    if not O.GO_OBO.exists():
        return _skip("go-basic.obo not downloaded")
    built = dict(O.build_go())
    built["hpo"] = O.build_hpo()
    for key, (nodes, anc, remap) in built.items():
        ref_nodes, ref_anc, ref_remap = O.load_ancestors(key)
        assert nodes == sorted(ref_nodes), key
        assert _canon(anc).equals(_canon(ref_anc)), key
        assert dict(zip(remap.source_id, remap.canonical_id)) == ref_remap, key


def test_d_max():
    for key, d in D_MAX.items():
        nodes, anc, _ = O.load_ancestors(key)
        assert int(O.distance_matrix(nodes, anc).max()) == d, key


def _compare(mine, ref, cols):
    ref = ref.copy()
    sw = ref.disease_a > ref.disease_b
    ref.loc[sw, ["disease_a", "disease_b"]] = ref.loc[sw, ["disease_b", "disease_a"]].to_numpy()
    m = mine.merge(ref, on=["disease_a", "disease_b"], suffixes=("", "_r"), how="outer", indicator=True)
    assert (m._merge == "both").all()
    for c in cols:
        a, b = m[c].to_numpy(), m[c + "_r"].to_numpy()
        assert (np.isnan(a) == np.isnan(b)).all(), c
        assert (a[~np.isnan(a)] == b[~np.isnan(b)]).all(), c


def test_go_disease_distance_bit_exact():
    ref_p = REPO_ROOT / "data" / "02_mechdist" / "disease_go_distance.parquet.gzip"
    if not ref_p.exists():
        return _skip("staged GO distance table absent")
    dg = pd.read_parquet(RAW / "associations" / "disease_gene.parquet")
    mine = DD.go_disease_distances(DD.pairs_sharing_entity(dg, "gene_id"), dg,
                                   pd.read_parquet(RAW / "genes" / "gene_go.parquet"))
    assert len(mine) == 114961
    _compare(mine, pd.read_parquet(ref_p), ["dist_BP", "dist_MF", "dist_CC", "dist_GO"])


def test_hpo_disease_distance_bit_exact():
    ref_p = REPO_ROOT / "data" / "zenodo" / "network" / "disease_hpo_distance.parquet.gzip"
    if os.environ.get("ONTOLOGY_PARITY_FULL") != "1" or not ref_p.exists():
        return _skip("HPO disease distance (set ONTOLOGY_PARITY_FULL=1; needs Zenodo network files)")
    ds = pd.read_parquet(RAW / "associations" / "disease_symptom.parquet")
    mine = DD.hpo_disease_distances(DD.pairs_sharing_entity(ds, "hpo_id"), ds)
    assert len(mine) == 13207984
    _compare(mine, pd.read_parquet(ref_p), ["dist_HPO"])


def test_hpo_baseline_matches_manuscript():
    """Appendix SI_ixidn_matched_ontology_baseline: 35% coverage, subset AUROC 0.69-0.72,
    DisAtlas higher by 0.03-0.06 (CI includes 0), all pairs 0.54."""
    pairs = pd.read_csv(HB.PAIRS, sep="\t")
    sc = HB.pair_scores(pairs)
    tab = HB.baseline_table(pairs[["blind_id", "a", "b", "y", "p", "pool"]].merge(sc, on="blind_id"))
    a = tab[tab.scope == "all"]
    assert round(a.coverage_pct.iloc[0]) == 35
    assert (round(a.auroc_scorable.min(), 2), round(a.auroc_scorable.max(), 2)) == (0.69, 0.72)
    assert (round(a.delta_model_minus.min(), 2), round(a.delta_model_minus.max(), 2)) == (0.03, 0.06)
    assert (a.d_ci_lo < 0).all()
    assert set(a.auroc_all_unscorable_last.round(2)) == {0.54}
    exp = {"hpo_resnik_bma": 0.7132, "hpo_lin_bma": 0.7193, "hpo_simgic": 0.6925, "mechdist_hpo": 0.7059}
    assert dict(zip(a.column, a.auroc_scorable)) == exp


# ---------------------------------------------------------------- ORDO / MONDO max-LCA and realism
OUT = REPO_ROOT / "data" / "00_rebuilt_ontology"
SI_SRC = REPO_ROOT / "figures" / "SI" / "source_data"


def _r(x, k):
    return round(float(x) + 1e-12, k)


def test_ordo_mondo_distance_matches_stored():
    import max_lca
    s = pd.read_csv(REPO_ROOT / "data" / "08S_eval" / "ixidn4364_scores.tsv", sep="\t")
    got = max_lca.pair_distances(s.a, s.b)
    ref = pd.read_csv(REPO_ROOT / "data" / "09S_ontology_baseline" / "ixidn_matched_ordo_mondo_distance.tsv", sep="\t")
    assert (ref.blind_id.values == s.blind_id.values).all()
    assert (got.ordo_dist.values == ref.ordo_dist.values).all()
    assert (got.mondo_dist.values == ref.mondo_dist.values).all()


def test_ontology_dependency_matches_manuscript():
    import build_ontology_dependency as BD
    BD.main()
    import hashlib
    # byte-identical to the stored production tables behind the appendix (and pinned by MD5)
    stored = REPO_ROOT / "data" / "09S_ontology_baseline"
    pairs = {"baselines": ("ontology_only_baselines.csv", "e0a51eb8182fff04df0b96d8ac9db837"),
             "correlations": ("ordo_hpo_correlations.csv", "922b3a984ae9620362aaa695ed8b8d31"),
             "discordant": ("ordo_hpo_least_concordant.csv", "7dedaba44e101ae075ee9bf795725ae3"),
             "orphanet": ("outside_orphanet_control.csv", "cd862c86bc48689687f115332358336f"),
             "annotation_source": ("hpo_annotation_orphanet_share.csv", None)}
    for k, (ref, h) in pairs.items():
        mine = (OUT / f"ontology_dependency_{k}.csv").read_bytes()
        assert mine == (stored / ref).read_bytes(), k
        assert h is None or hashlib.md5(mine).hexdigest() == h, k
    b = pd.read_csv(OUT / "ontology_dependency_baselines.csv")
    a = pd.read_csv(OUT / "ontology_dependency_annotation_source.csv").iloc[0]
    assert a.pct_annotations_on_orphanet_diseases == 99.1                       # "99.1% ... trace back to Orphanet"
    c = pd.read_csv(OUT / "ontology_dependency_correlations.csv")
    rho = c[(c.x == "neg_ordo_dist") & c.y.isin(["mechdist_hpo", "hpo_simgic"]) & (c.subset == "pooled")].spearman_rho
    assert (_r(rho.min(), 2), _r(rho.max(), 2)) == (0.51, 0.62)                  # Spearman 0.51-0.62
    d = pd.read_csv(OUT / "ontology_dependency_discordant.csv"); dd = d[d.hpo != "-"]
    assert (_r(dd.model.min(), 2), _r(dd.model.max(), 2)) == (0.70, 0.76)        # DisAtlas 0.70-0.76
    assert (_r(dd.hpo_auroc.min(), 2), _r(dd.hpo_auroc.max(), 2)) == (0.67, 0.73)  # HPO 0.67-0.73
    o = dd[dd.ontology == "neg_ordo_dist"].ontology_auroc                       # ORDO "falls to chance"
    assert o.min() > 0.44 and o.max() < 0.54
    r = pd.read_csv(OUT / "ontology_dependency_orphanet.csv").set_index("subset")
    out = r.loc["at least one disease outside ORDO"]
    # manuscript prints 62% and 0.786; exact values are 62.51% (2,728 / 4,364) and 0.78652, i.e. the
    # text truncates rather than rounds (rounded: 63%, 0.787). Asserted at the stored precision.
    assert out.n == 2728 and int(100 * out.n / 4364) == 62 and out.model == 0.7865
    assert _r(r.loc["outside ORDO, held-out test-test"].model, 3) == 0.810       # 0.810 test split
    h = b[(b.scope == "all") & b.column.str.startswith(("hpo_", "mechdist_hpo"))]
    assert set(h.coverage_pct.round().astype(int)) == {35}                         # 35% of pairs
    assert (_r(h.auroc_scorable.min(), 2), _r(h.auroc_scorable.max(), 2)) == (0.69, 0.72)
    assert (_r(h.delta_model_minus.min(), 2), _r(h.delta_model_minus.max(), 2)) == (0.03, 0.06)
    assert (h.d_ci_lo < 0).all()                                                   # not significant
    assert set(h.auroc_all_unscorable_last.map(lambda x: _r(x, 2))) == {0.54}     # all pairs 0.54


def test_negative_realism_matches_figure_source_and_manuscript():
    import build_negative_realism as NR
    NR.main()
    for f in ("negative_realism_hard_negative_auroc.tsv", "negative_realism_ordo_distance.tsv.gz"):
        assert pd.read_csv(OUT / f, sep="\t").equals(pd.read_csv(SI_SRC / f, sep="\t")), f
    m = pd.read_csv(OUT / "negative_realism_summary.tsv", sep="\t").set_index("set")
    assert m.loc["IxIDN-Matched negatives", "median"] == 8 and m.loc["IxIDN-Matched positives", "median"] == 6
    assert round(m.loc["IxIDN-Matched negatives", "pct_share_organ_system"]) == 37
    assert round(m.loc["IxIDN-Matched positives", "pct_share_organ_system"]) == 37
    h = pd.read_csv(OUT / "negative_realism_hard_negative_auroc.tsv", sep="\t")
    a = h[h.scope == "All 4,364 pairs"].set_index("subset")
    assert _r(a.loc["Share an organ system with each other", "auroc_disatlas"], 2) == 0.72
    assert _r(a.loc["MONDO distance <= positive median (ontology-close)", "auroc_disatlas"], 2) == 0.72
    assert _r(a.loc["ORDO distance <= positive median (ontology-close)", "auroc_disatlas"], 2) == 0.66
    assert _r(a.loc["GO Resnik gene-BMA >= positive median (function-similar)", "auroc_disatlas"], 2) == 0.67
    hp = a.loc["HPO Resnik-BMA >= positive median (phenotype-similar)"]
    assert round(hp.pct_neg_kept) == 7 and _r(hp.auroc_disatlas, 2) == 0.54
    org = a.loc["Both classes share an organ system"]
    assert (_r(org.auroc_disatlas, 3), _r(org.auroc_hpo_resnik_bma, 3), _r(org.auroc_go_resnik_gene_bma, 3)) == \
        (0.780, 0.585, 0.641)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
