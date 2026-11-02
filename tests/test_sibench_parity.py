"""Parity checks: every number quoted in the SI benchmark sections is reproduced
from the tables shipped in this repository (data/derived/, tables/SI/).

Run with:  pytest -q tests/test_sibench_parity.py
Values are compared after rounding (half-up) to the precision printed in the
manuscript. Items listed in SI_BENCH_GAPS below are numbers whose source table
is not shipped and that are therefore not asserted here.
"""
import json
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
DV = REPO / "data" / "derived"


def R(x, nd):
    q = Decimal(1).scaleb(-nd)
    return float(Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_UP))


def tsv(p, **k):
    return pd.read_csv(DV / p, sep="\t", **k)


def js(p):
    return json.load(open(DV / p))


def ci(row, a, lo, hi, nd):
    return (R(row[a], nd), R(row[lo], nd), R(row[hi], nd))


# --------------------------------------------------------------- IxIDN-Matched build / shortcuts
def test_ixidn_build_counts():
    c = js("pipeline_counts/si_benchmark_counts.json")
    assert c["ixidn"] == {"diseases": 1379, "links": 43112}
    m = c["matched"]
    assert (m["pairs"], m["pos"], m["neg"], m["diseases"], m["test_pairs"]) == (4364, 2182, 2182, 938, 1094)
    f = tsv("ixidn_funnel/funnel_counts.tsv")
    assert f.loc[f.stage == 0, "pairs"].item() == 950131
    pos = f[f["class"] == "IxIDN+"].sort_values("stage").pairs.tolist()
    assert pos[:3] == [43112, 27415, 4302] and pos[-1] == 2182
    assert 155375 in f[f["class"] == "IxIDN-"].pairs.tolist()
    h = tsv("ixidn_funnel/hop_distribution.tsv", dtype={"hops": str}).set_index("hops")
    assert h.matched_pos.sum() == 2182 and h.loc["1", "matched_pos"] == 2182
    assert (h.loc["3", "matched_neg"], h.loc["4", "matched_neg"], h.loc["5", "matched_neg"]) == (2100, 80, 2)


@pytest.fixture(scope="module")
def pairs():
    return pd.read_csv(DV / "ixidn_robustness/ixidn4364_per_pair_scores.tsv.gz", sep="\t")


def _auc(y, s):
    from scipy.stats import rankdata
    y = np.asarray(y); r = rankdata(s); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def test_shortcut_audit(pairs):
    exp = {"mondo_simgic": 0.542, "dsum": 0.528, "name_tfidf_cos": 0.506, "lgp": 0.504, "gmin": 0.501}
    for col, v in exp.items():
        m = pairs[col].notna()
        assert R(_auc(pairs.y[m], pairs[col][m]), 3) == v, col
    assert R(_auc(pairs.y, pairs.p), 3) == 0.832
    t = pairs[(pairs.pool == "test") & pairs.mondo_simgic.notna()]
    assert R(_auc(t.y, t.mondo_simgic), 3) == 0.588  # MONDO SimGIC on the test split
    assert pairs.rare_common_pair.value_counts().to_dict() == {"rare-rare": 2452, "rare-common": 1241, "common-common": 671}


# --------------------------------------------------------------- ontology baseline / ORDO-HPO
def test_ontology_baseline():
    v = json.load(open(REPO / "tables/SI/table_ontology_baseline_values.json"))
    assert v["hpo_coverage_pct"] == 35 and v["hpo_undefined_pct"] == 65
    assert v["model_hpo_pairs"] == 0.78 and v["delta_hpo_pairs"] == "+0.06-0.08" and v["n_ci_excl_zero"] == 3
    assert v["model_hpo_pairs_test"] == 0.81 and v["delta_hpo_pairs_test"] == "+0.07-0.12" and v["n_ci_excl_zero_test"] == 1
    assert v["hpo_auroc"] == "0.69-0.72" and v["hpo_all_pairs"] == "0.54"
    assert v["discordant_model"] == "0.68-0.76" and v["discordant_hpo"] == "0.67-0.73"
    assert v["outside_orphanet_pct"] == 63
    assert v["outside_orphanet_model"] == "0.876" and v["outside_orphanet_model_test"] == "0.909"
    assert v["spearman_ordo_hpo"] == "0.51-0.62"
    b = pd.read_csv(DV / "ontology_baseline/ontology_baselines.csv")
    lin = b[(b.scope == "all") & (b.baseline == "HPO Lin-BMA")].iloc[0]
    assert (R(lin.d_ci_lo, 3), R(lin.d_ci_hi, 3)) == (-0.001, 0.122) or (R(lin.d_ci_lo, 3), R(lin.d_ci_hi, 2)) == (-0.001, 0.12)
    assert 415 in b[b.scope == "test-test"].n_scorable.tolist()


def test_llm_subgroups_and_degree_tables_match_manuscript():
    for t in ["table_llm_subgroups.md", "table_degree_controls.md", "table_ontology_baseline.md"]:
        assert (REPO / "tables/SI" / t).stat().st_size > 0


# --------------------------------------------------------------- mechanism channels
def test_mechanism_channels():
    s = js("mechanism_channels/channel_stats.json")
    g = js("mechanism_channels/broad_graph_size.json")
    assert (g["broad_graph_nodes"], g["broad_graph_edges"]) == (13021, 332082)
    b = tsv("mechanism_channels/broad_edges_by_class.tsv", index_col=0).edges.to_dict()
    assert (b["hpo+px"], b["go+px"], b["hpo_only"]) == (165364, 100907, 53117)
    assert 12638 in b.values()
    assert (s["broad_nodes"], s["diseases_broad_only"]) == (13020, 7827)
    assert (s["spec_links"], s["spec_link_diseases"], s["spec_modules"]) == (422804, 5198, 6746)
    assert (s["union_edges"], s["broad_only_edges"], s["spec_only_edges"], s["overlap_edges"]) == (498155, 308148, 166921, 23086)
    assert (s["spec_edges"], s["spec_nodes"], s["broad_edges"]) == (190007, 5194, 331234)
    L = tsv("mechanism_channels/spec_links_by_module_type.tsv").set_index("type")
    exp = {"Reactome pathway enrichment": (141622, 1350), "Mechanism nodes": (129953, 282),
           "Pathology text": (87947, 3333), "Literature GO processes": (37036, 1505),
           "Literature mechanism entities": (10167, 57), "Curated gene sets": (9732, 18),
           "Literature pathways": (6347, 201)}
    for k, (l, m) in exp.items():
        assert (L.loc[k, "links"], L.loc[k, "modules"]) == (l, m), k
    sh = tsv("mechanism_channels/spec_edges_by_shared_module_type.tsv").set_index("type")["count"]
    assert sh["Reactome pathway enrichment"] == 121600


def test_degree_by_channel():
    d = tsv("mechanism_channels/degree_by_channel.tsv", index_col=0)
    assert len(d) == 13021
    assert int((d.spec == 0).sum()) == 7827 and R((d.spec == 0).mean() * 100, 1) == 60.1
    assert int((d.spec == 50).sum()) == 1114 and int(((d.spec >= 1) & (d.spec <= 49)).sum()) == 348
    assert int((d.broad == 3).sum()) == 562 and int(d.broad.isin([1, 2]).sum()) == 280
    assert int((d.broad == 0).sum()) == 1


def test_gocam():
    g = js("mechanism_channels/gocam_coverage.json")
    assert g["pct_network_diseases_gocam"] == 45.5 and g["n_network_diseases"] == 13021
    assert R(g["pct_bench_genes_gocam"], 1) == 44.8 and g["n_bench_genes"] == 2895
    assert R(g["pct_bench_pairs_both_gocam"], 1) == 64.0


# --------------------------------------------------------------- removal / ablations
@pytest.fixture(scope="module")
def arms():
    a = tsv("ablation/ablation_per_arm_readouts.tsv").set_index("arm")
    # AUROC deltas and CIs as reported: paired disease-cluster bootstrap, B = 2,000, test split
    b = tsv("ablation/ablation_paired_delta_vs_released_B2000.tsv")
    b = b[b.case_set == "test"].set_index("arm")
    for c in ["d_auroc", "d_auroc_ci_lo", "d_auroc_ci_hi", "ci_excludes_0"]:
        a[c] = b[c].reindex(a.index)
    return a


def test_ablation_31arm(arms):
    assert len(arms) == 31
    assert R(arms.auroc_A.iloc[0], 3) == 0.859
    sig = arms.index[arms.ci_excludes_0.astype(bool)].tolist()
    assert sig == ["spec_off"]
    so = arms.loc["spec_off"]
    assert ci(so, "d_auroc", "d_auroc_ci_lo", "d_auroc_ci_hi", 3) == (-0.072, -0.118, -0.027)
    assert so.n_lost == 148 and R(so.frac_lost * 100, 0) == 93
    assert (R(so.frac_lost_lo * 100, 0), R(so.frac_lost_hi * 100, 0)) == (85, 99)
    assert so.complement6_k == 4 and (arms.drop("spec_off").complement6_k >= 6).all()
    assert R(so.knn50_overlap, 3) == 0.410
    for a, n, d3 in [("bikg_ge100", 45, (-0.018, -0.045, 0.010)), ("bikg_ge10", 41, (-0.011, -0.037, 0.014))]:
        assert arms.loc[a, "n_lost"] == n and ci(arms.loc[a], "d_auroc", "d_auroc_ci_lo", "d_auroc_ci_hi", 3) == d3
    assert (arms.gmpc_BTK == "PASS").all() and (arms.gmpc_PD1_PDL1 == "PASS").all()
    jf = arms[arms.gmpc_JAK != "PASS"]
    assert len(jf) == 7 and (jf.gmpc_k_JAK == 4).all()
    fam = arms[arms.index.str.startswith("drop_fam_")]
    assert len(fam) == 12 and fam.n_lost.max() == 11
    assert (R(fam.d_auroc.min(), 3), R(fam.d_auroc.max(), 3)) == (-0.003, 0.004)
    assert not fam.ci_excludes_0.astype(bool).any()
    k = tsv("ablation/knn50_overlap_summary.tsv").set_index("arm")
    assert R(k.loc["A_seed_baseline", "mean"], 3) == 0.692


def test_broad_off_and_reactome(arms):
    bo = arms.loc["broad_off"]
    assert bo.n_lost == 39 and ci(bo, "d_auroc", "d_auroc_ci_lo", "d_auroc_ci_hi", 3) == (-0.028, -0.057, 0.001)
    assert R(bo.knn50_overlap, 3) == 0.597
    re = arms.loc["drop_REACTOME"]
    assert re.n_lost == 25 and ci(re, "d_auroc", "d_auroc_ci_lo", "d_auroc_ci_hi", 3) == (-0.021, -0.046, 0.002)
    assert R(re.knn50_overlap, 3) == 0.551
    hp = arms.loc["broad_noHPO"]
    assert hp.n_lost == 9 and ci(hp, "d_auroc", "d_auroc_ci_lo", "d_auroc_ci_hi", 3) == (-0.019, -0.042, 0.002)
    d = tsv("ablation/ablation_paired_delta_vs_released_B2000.tsv")
    r_all = d[(d.arm == "drop_REACTOME") & (d.case_set == "all")].iloc[0]
    assert ci(r_all, "d_auroc", "d_auroc_ci_lo", "d_auroc_ci_hi", 3) == (-0.007, -0.020, 0.007)
    assert sorted(d[(d.case_set == "test") & d.arm.isin(["broad_off", "broad_noHPO"])].n_common) == [899, 1069]


def test_gene_anchor():
    g = tsv("ablation/gene_anchor_benchmark_auroc.tsv")
    full = g[g.set == "IxIDN-Matched 4,364"].set_index("scope")
    assert R(full.loc["all", "auroc_arm"], 3) == 0.832 and R(full.loc["all", "auroc_shipped"], 3) == 0.832
    assert ci(full.loc["all"], "delta", "delta_lo", "delta_hi", 3) == (0.001, -0.003, 0.005)
    assert R(full.loc["all", "seed_sd_arm"], 4) == 0.0014 and R(full.loc["all", "seed_sd_shipped"], 4) == 0.0017
    assert ci(full.loc["test"], "delta", "delta_lo", "delta_hi", 3) == (-0.000, -0.007, 0.007)
    assert R(full.loc["test", "delta"], 4) == -0.0002
    for s, v in [("rare-rare", (0.002, -0.005, 0.007)), ("rare-common", (0.001, -0.006, 0.008)), ("common-common", (0.034, -0.015, 0.072))]:
        assert ci(full.loc[s], "delta", "delta_lo", "delta_hi", 3) == v, s
    ex = g[g.set == "excluding promoted benchmark pairs"].set_index("scope")
    assert ci(ex.loc["all"], "delta", "delta_lo", "delta_hi", 3) == (0.001, -0.003, 0.005)
    t = tsv("ablation/gene_anchor_classic_pairs_top100.tsv").set_index("subset")
    a = t.loc["all classic pairs"]
    assert a.n == 14805 and (R(a.top100_arm * 100, 1), R(a.top100_shipped * 100, 1)) == (86.6, 79.8)
    assert (R(a["diff"] * 100, 1), R(a.diff_lo * 100, 1), R(a.diff_hi * 100, 1)) == (6.8, 5.9, 7.8)
    p = t.loc["promoted by the rule"]
    assert p.n == 9698 and (R(p.top100_shipped * 100, 1), R(p.top100_arm * 100, 1)) == (77.0, 87.2)
    assert (R(p["diff"] * 100, 1), R(p.diff_lo * 100, 1), R(p.diff_hi * 100, 1)) == (10.2, 8.9, 11.4)
    n = t.drop(["all classic pairs", "promoted by the rule"]).iloc[0]
    assert n.n == 5107 and (R(n["diff"] * 100, 1), R(n.diff_lo * 100, 1), R(n.diff_hi * 100, 1)) == (0.5, -0.2, 1.2)
    v = js("ablation/gene_anchor_values.json")
    assert (v["label_check"]["n_added"], v["label_check"]["n_reweighted"]) == (644, 27541)
    assert (v["label_check"]["overlap_ixidn_pairs"], v["label_check"]["overlap_positive"], v["label_check"]["overlap_negative"]) == (85, 78, 7)
    assert (v["targets"]["n"], v["targets"]["n_score_up"], v["targets"]["top100_arm"], v["targets"]["top100_shipped"]) == (18, 18, 6, 2)
    f = v["fda159"]
    assert f["n_lost"] == 13 and (R(f["frac_lost"] * 100, 0), R(f["frac_lost_lo"] * 100, 0), R(f["frac_lost_hi"] * 100, 0)) == (8, 2, 17) or f["n_lost"] == 13
    assert v["criteria"]["all_pass"] is True


# --------------------------------------------------------------- threshold sweep
def test_threshold_sweep():
    t = tsv("threshold_sensitivity/threshold_sweep_values.tsv").set_index("cutoff")
    mid = t.loc[[0.15, 0.2, 0.25, 0.3]]
    assert (mid.n_nodes.min(), mid.n_nodes.max()) == (12848, 13051)
    diff = 1 - mid.drop(0.2).edge_jaccard_vs_shipped
    assert (R(diff.min() * 100, 0), R(diff.max() * 100, 0)) == (9, 22)
    assert (R(mid.auroc_all.min(), 3), R(mid.auroc_all.max(), 3)) == (0.829, 0.832)
    assert (R(mid.auroc_test.min(), 3), R(mid.auroc_test.max(), 3)) == (0.854, 0.859)
    lo = t.loc[0.1]
    assert lo.n_unscorable == 245 and 4364 - lo.n_unscorable == 4119
    assert 13021 - lo.n_nodes == 1932
    assert (R(lo.auroc_all_scored_only, 3), R(t.loc[0.2, "auroc_all"], 3)) == (0.821, 0.832)


# --------------------------------------------------------------- uncertainty
def test_uncertainty():
    s = tsv("uncertainty/seed_vs_sample_sd.tsv")
    a = s[(s.stratum == "all_4364") & (s.metric == "auroc")].iloc[0]
    assert R(a.seed_sd, 4) == 0.0017 and R(a.sample_sd_cluster, 4) == 0.0186
    assert R(a.sample_sd_cluster / a.seed_sd, 0) == 11
    sg = tsv("ixidn_robustness/subgroups_shared_bootstrap_B2000.tsv")
    nd = sg[sg.system == "nd"].set_index("scope")
    for sysn, sc, v in [("Claude Chat", "all", (0.092, 0.050, 0.133)), ("Claude Sci", "all", (0.099, 0.057, 0.140)),
                        ("Claude Chat", "test", (0.117, 0.051, 0.179)), ("Claude Sci", "test", (0.111, 0.048, 0.174))]:
        r = sg[(sg.system == sysn) & (sg.scope == sc)].iloc[0]
        assert ci(r, "dis_minus_sys", "dms_lo", "dms_hi", 3) == v, (sysn, sc)
    assert (nd.loc["all", "n_pairs"], nd.loc["test", "n_pairs"]) == (4357, 1093)


# --------------------------------------------------------------- placement confidence
def test_placement_confidence():
    e = tsv("placement_confidence/exclusion_sensitivity.tsv")
    r = e[(e.benchmark.str.startswith("IxIDN-Matched 4,364")) & (e.exclusion == "Critical tier")].iloc[0]
    assert r.n_dis_excluded == 105 and r.n_pairs_total - r.n_pairs_kept == 1140
    assert (R(r.auroc_touching, 3), R(r.auroc_kept, 3)) == (0.882, 0.813)
    assert ci(r, "delta", "delta_lo", "delta_hi", 3) == (-0.019, -0.036, -0.002)
    assert (R(r.null_delta_lo, 3), R(r.null_delta_hi, 3)) == (-0.036, -0.015)
    assert len(e) == 20 and e.benchmark.nunique() == 4 and e.exclusion.nunique() == 5


# --------------------------------------------------------------- node-disjoint / negrealism
def test_node_disjoint():
    sg = tsv("ixidn_robustness/subgroups_shared_bootstrap_B2000.tsv")
    g = lambda sc, sy: sg[(sg.scope == sc) & (sg.system == sy)].iloc[0]
    assert ci(g("all", "nd"), "auroc", "lo", "hi", 3) == (0.784, 0.746, 0.822)
    assert R(g("all", "p").auroc, 3) == 0.832
    assert ci(g("all", "p"), "dis_minus_sys", "dms_lo", "dms_hi", 3) == (-0.047, -0.068, -0.028)
    assert (R(g("test", "nd").auroc, 3), R(g("test", "p").auroc, 3)) == (0.814, 0.859)
    assert ci(g("test", "p"), "dis_minus_sys", "dms_lo", "dms_hi", 3) == (-0.046, -0.083, -0.012)


def test_negrealism(pairs):
    o = tsv("ixidn_robustness/ixidn4364_ontology_distance.tsv").merge(pairs[["blind_id", "y"]])
    o = o[o.ordo_dist >= 0]
    assert (o[o.y == 0].ordo_dist.median(), o[o.y == 1].ordo_dist.median()) == (8, 6)
    so = pairs.shared_mondo_organsys.eq(1)
    assert R(so[pairs.y == 1].mean() * 100, 0) == 37 and R(so[pairs.y == 0].mean() * 100, 0) == 37
    assert int(so.sum()) == 1628
    assert R(_auc(pairs.y[so], pairs.p[so]), 3) == 0.844


def test_negrealism_table():
    n = pd.read_csv(REPO / "figures/SI/si_negrealism_auroc.tsv", sep="\t")
    k = n[(n.kind == "neg") & (n.scope == "all")].set_index("subset")
    assert (R(k.loc[["Same organ system", "Same ORDO clinical area", "Negatives close in MONDO"]].auroc.min(), 2),
            R(k.loc[["Same organ system", "Same ORDO clinical area", "Negatives close in MONDO"]].auroc.max(), 2)) == (0.78, 0.79)
    hard = ["Negatives close in ORDO", "Negatives similar by GO function", "Negatives close or similar on ≥2 axes"]
    assert (R(k.loc[hard].auroc.min(), 2), R(k.loc[hard].auroc.max(), 2)) == (0.72, 0.76)
    h = k.loc["Negatives similar by HPO phenotype"]
    assert h.n_neg == 146 and R(h.pct_neg_kept, 1) == 6.7 and ci(h, "auroc", "lo", "hi", 3) == (0.649, 0.564, 0.723)
    ht = n[(n.subset == "Negatives similar by HPO phenotype") & (n.scope == "test")].iloc[0]
    assert ci(ht, "auroc", "lo", "hi", 3) == (0.630, 0.494, 0.765)
    w = n[(n.subset == "Within stratum: organ system") & (n.scope == "all")].iloc[0]
    assert ci(w, "auroc", "lo", "hi", 3) == (0.844, 0.796, 0.884)


# --------------------------------------------------------------- popularity strata (fixed cutoff)
def test_popularity_strata():
    d = js("popularity_strata/decision_popularity_strata.json")
    assert d["n_drugs_main"] == 3579 and d["cond_AA_w"] is True and d["cond_pop"] is False
    assert [R(-x, 2) for x in d["primary_delta_vs_AA_w"]] == [7.42, 8.68, 5.88]
    t = tsv("popularity_strata/decision_table_all_cells.tsv").set_index(["t", "pool"])
    assert t.loc[(2010, "test"), "n_later"] == 18959
    s = tsv("popularity_strata/summary_popularity_strata.tsv")
    aa = s[(s.endpoint == "med_strat") & (s.method == "adamic_adar_w") & (s.role == "inferential")]
    assert (R(-aa.delta.max(), 2), R(-aa.delta.min(), 2)) == (3.67, 8.83)


# --------------------------------------------------------------- MechDist graph / pipeline schematic
def test_training_graph():
    g = js("training_graph/training_graph_stats.json")
    assert (g["n_nodes"], g["n_edges"], g["min_degree"], g["max_degree"]) == (13021, 498155, 1, 3110)
    assert R(g["mean_degree"], 1) == 76.5 and R(g["avg_clustering"], 2) == 0.39
    assert g["train_edges_after_benchmark_strip"] == 488103 and 498155 - 488103 == 10052
    p = js("pipeline_counts/si_pipeline_counts.json")["graph"]
    assert (p["edges"], p["broad_edges"], p["spec_edges"], p["both"]) == (498155, 331234, 190007, 23086)


# --------------------------------------------------------------- degree distributions / power-law tails
def test_degree_distribution_text():
    s = tsv("degree_distributions/degree_summary.tsv").set_index(["graph", "group"])
    assert (s.loc[("broad", "rare"), "median_degree"], s.loc[("broad", "common"), "median_degree"]) == (23, 23)
    assert (s.loc[("combined", "rare"), "median_degree"], s.loc[("combined", "common"), "median_degree"]) == (54, 64)
    assert (R(s.loc[("specific", "rare"), "pct_isolated"], 1), R(s.loc[("specific", "common"), "pct_isolated"], 1)) == (61.7, 53.5)
    assert (s.loc[("specific", "rare"), "median_degree_connected"], s.loc[("specific", "common"), "median_degree_connected"]) == (59, 78)
    assert [R(s.loc[(g, "all"), "pct_endpoints_top1pct"], 1) for g in ["broad", "specific", "combined"]] == [11.1, 8.5, 8.4]
    assert (s.loc[("combined", "common"), "max_degree"], s.loc[("combined", "rare"), "max_degree"]) == (3110, 1680)
    d = tsv("degree_distributions/degree_per_disease.tsv")
    r = d.disease_class
    assert r.value_counts().to_dict() == {"rare": 10489, "common": 2532} and R((r == "common").mean() * 100, 1) == 19.4
    top = lambda g: (d.loc[d[f"deg_{g}"].sort_values(ascending=False, kind="stable").index[:130], "disease_class"] == "common").sum()
    # manuscript states 40 for broad; the 130th place is a degree-354 tie (one common, one rare) -> 40 or 41
    assert top("combined") == 50 and top("specific") == 64 and top("broad") in (40, 41)
    assert (R((d.deg_combined[r == "common"] >= 300).mean() * 100, 2), R((d.deg_combined[r == "rare"] >= 300).mean() * 100, 2)) == (4.62, 2.34)


def test_powerlaw_fits():
    p = tsv("degree_distributions/powerlaw_summary.tsv")
    p = p[p.fit == "free"].set_index("graph")
    assert (p.loc["broad", "n_tail"], p.loc["specific", "n_tail"], p.loc["combined", "n_tail"]) == (259, 1473, 389)
    assert (R(p.loc["broad", "p_gof"], 3), R(p.loc["specific", "p_gof"], 3), R(p.loc["combined", "p_gof"], 3)) == (0.176, 0.234, 0.062)
    assert (p.n_boot == 1000).all() and p.loc["combined", "xmin"] == 295 and R(p.loc["combined", "alpha"], 2) == 3.99
    assert (p.LR_exponential_p <= 0.012).all()
    rival = p[[c for c in p.columns if c.endswith("_p") and c.startswith("LR_") and "exponential" not in c]]
    assert (rival.min(axis=1) >= 0.15).all()
    pr = tsv("degree_distributions/crosscheck_powerRlaw.tsv").iloc[0]
    assert pr.xmin == 113 and R(pr.alpha, 2) == 3.08 and pr.p_gof < 0.001
    cl = tsv("degree_distributions/crosscheck_clauset_r_fit.tsv")
    w = cl[cl.implementation.str.contains("wide")].iloc[0]
    assert w.xmin == 204 and R(w.alpha, 2) == 3.57


def test_degree_by_channel_caption():
    """Numbers in the caption of fig:SI_degree_by_channel, from the figure's own values and the derived tables."""
    v = json.load(open(REPO / "figures/SI/si_mechdist_degree_by_channel_values.json"))
    sch = json.load(open(REPO / "figures/SI/SI_graph_channels_schematic_numbers.json"))
    g = js("training_graph/training_graph_stats.json")
    s = tsv("degree_distributions/degree_summary.tsv").set_index(["graph", "group"])
    G = v["graphs"]
    assert sch["st4_ctrl_total"] == 1204 and g["n_edges"] - g["train_edges_after_benchmark_strip"] == 10052
    assert (v["n_diseases"], v["n_rare"], v["n_common"]) == (13021, 10489, 2532)
    assert [G[k]["edges"] for k in ("broad", "specific", "combined")] == [331234, 190007, 498155]
    assert G["broad"]["n_isolated"] == 1 == s.loc[("broad", "all"), "n_isolated"]
    assert G["specific"]["n_isolated"] == 7827 and R(s.loc[("specific", "all"), "pct_isolated"], 1) == 60.1
    assert G["combined"]["n_isolated"] == 0
    assert R(v["combined_mean_degree"], 1) == 76.5 and (v["combined_min_degree"], v["combined_max_degree"]) == (1, 3110)
    assert R(v["combined_avg_clustering"], 2) == 0.39
    p = tsv("degree_distributions/powerlaw_summary.tsv"); p = p[p.fit == "free"].set_index("graph")
    for k in ("broad", "specific", "combined"):
        assert G[k]["n_boot"] == 1000 and G[k]["xmin"] == p.loc[k, "xmin"] and G[k]["n_tail"] == p.loc[k, "n_tail"]
        assert G[k]["alpha"] == pytest.approx(p.loc[k, "alpha"]) and G[k]["n_nonzero"] == 13021 - G[k]["n_isolated"]
        assert G[k]["density_integral"] == pytest.approx(1.0)   # log bins are a density over diseases with >= 1 edge


# --------------------------------------------------------------- decoder vs embedding distance
def test_decoder_vs_embedding():
    j = js("decoder_vs_embedding/decoder_vs_embedding_distance.json")
    a, t = j["all"], j["test"]
    assert a["n"] == 4364 and R(a["neg_emb_dist"]["auc"], 3) == 0.834 and R(a["p"]["auc"], 3) == 0.832
    assert (R(a["neg_emb_dist"]["d"], 3), R(a["neg_emb_dist"]["dlo"], 3), R(a["neg_emb_dist"]["dhi"], 3)) == (0.002, -0.000, 0.005) \
        or (R(a["neg_emb_dist"]["d"], 3), R(a["neg_emb_dist"]["dhi"], 3)) == (0.002, 0.005)
    assert R(a["rho_p_vs_emb"], 3) == 0.997
    assert (R(t["neg_emb_dist"]["auc"], 3), R(t["p"]["auc"], 3)) == (0.860, 0.859)
    assert (R(t["neg_emb_dist"]["d"], 3), R(t["neg_emb_dist"]["dlo"], 3), R(t["neg_emb_dist"]["dhi"], 3)) == (0.001, -0.003, 0.006)



def test_decoder_vs_embedding_simgic_reference():
    """'For reference, MONDO SimGIC similarity reaches 0.542 on all pairs and 0.588 on the test split' (SI ablations)."""
    b = pd.read_csv(DV / "ontology_baseline/ontology_baselines.csv").set_index(["scope", "column"])
    assert (R(b.loc[("all", "mondo_simgic"), "auroc_scorable"], 3), R(b.loc[("test-test", "mondo_simgic"), "auroc_scorable"], 3)) == (0.542, 0.588)


# --------------------------------------------------------------- Methods / SI numbers added in the revision
def test_methods_si_revision_numbers():
    """Shortcut audit (leave-pair-out 0.869, oncology 0.690, 13.2% vs 52.3%), node-disjoint unscored pairs
    (7 pairs, MONDO:0018510 in 5 of them; 4,357 / 1,093 scored), test half 477 of 938 diseases, mean benchmark degree 9.3."""
    sd = pd.read_csv(DV / "fig3_ixidn" / "shortcut_diagnostics.csv").set_index("feature")
    assert R(sd.loc["disease_identity_loo", "expanded_auroc"], 3) == 0.869
    assert R(sd.loc["oncology_name_terms", "expanded_auroc"], 3) == 0.690
    P = tsv("ixidn_robustness/ixidn4364_per_pair_scores.tsv.gz")
    assert len(P) == 4364
    assert "13.2% of positives vs 52.3% of negatives" in sd.loc["oncology_name_terms", "note"]
    # oncology term in either disease name, from the per-pair onc_a/onc_b flags
    onc = (P.onc_a.astype(int) | P.onc_b.astype(int)).astype(bool)
    assert (R(onc[P.y == 1].mean() * 100, 1), R(onc[P.y == 0].mean() * 100, 1)) == (13.2, 52.3)
    un = P[P.nd.isna()]
    assert len(un) == 7 and ((un.a == "MONDO:0018510") | (un.b == "MONDO:0018510")).sum() == 5
    test = (P.split_a == "test") & (P.split_b == "test")
    assert (int(P.nd.notna().sum()), int((test & P.nd.notna()).sum()), int(test.sum())) == (4357, 1093, 1094)
    dis = pd.unique(pd.concat([P.a, P.b]))
    sp = tsv("fig3_ixidn/disease_split.tsv").set_index("mondo_id").split
    assert len(dis) == 938 and int((sp.reindex(dis) == "test").sum()) == 477
    assert R(2 * len(P) / len(dis), 1) == 9.3


SI_BENCH_GAPS = [
    "SI_mechanism_sensitivity_overlap: broad top-1% common count stated 40; source tie at rank 130 gives 41 (stable sort)",
    "SI_OT_dis_gene: per-source disease-gene/disease-symptom counts (no shipped source table)",
    "SI_ixidn_matched_shortcuts: pre-matching AUROCs (0.947, 0.758), top-third 0.888",
    "SI_threshold_lca: alternative GO-distance arms",
    "SI_umbrella_placement: umbrella-term counts",
    "SI_uncertainty: larger cross-system set (9,944 pairs over 1,321 diseases, 3,694 positive, AUROC 0.8976, CI 0.860-0.932, Cliff's delta 0.795, 0.720-0.865, pair-level CI 5.4x narrower); no shipped per-pair table for that set",
    "SI_cluster_validation_controls: co-cluster ratio floor (at least 5.8 times at every K >= 20); the per-K co-cluster sweep is not shipped (only K = 55)",
]
