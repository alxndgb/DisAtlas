"""Maintainer-only: stage the stored inputs of the SI figures whose final producer is the gene-anchored cosine audit
(Supplementary Figures on uncertainty, placement confidence, popularity strata and FDA channels, and the
language-model subgroup table).  Nothing is
recomputed except the two small uncertainty input tables, which are column selections of stored results.
Source locations are passed as KEY=PATH and are never tracked.

usage: python scripts/maintainer/stage_figD.py <data dir> KEY=PATH ...
keys:  unc_svs      seed-vs-sample uncertainty table of the gene-anchored cosine (uncertainty_seed_vs_sample_cos_B2000.tsv)
       unc_g1       headline AUROC table, B = 2,000 disease-cluster bootstrap (headline_auroc_g1_B2000_geneanchored.tsv; column c)
       placement    exclusion-sensitivity table of the placement-confidence analysis (dq_exclusion_sensitivity.tsv)
       pop_dir      directory with summary_popularity_strata.tsv, strata_boundaries.tsv, decision_popularity_strata.json
       fda_dir      directory with broad_only_fda_dev_results.tsv and broad_only_fda_dev_pair_scores_ens5.tsv
       fda_ref_raw  an earlier unblinded per-pair FDA table (same pairs, same order as the shipped blinded table) that
                    lets the shipped tokens be reused: data/derived/si_fda/broad_only_fda_dev_pair_scores_ens5.tsv is
                    the blinded counterpart.  Diseases carrying a token keep the same token and an empty label.
       fda_prior    directory with the dated FDA prioritisation tables endpoints.tsv and paired_differences.tsv
       llm_sub      shared-draw subgroup bootstrap of DisAtlas against the language-model systems (llm_subgroups_cosine_B2000.tsv)
       baselines_dir directory with blockA_baselines_g1_B2000_geneanchored.tsv, heldout_baselines_g1_B2000_geneanchored.tsv,
                    own_edges_contrasts_g1_B2000_geneanchored.tsv and OA_direct_measures.tsv (SI baselines table)
       heur_ai      paired differences of training-free scores minus each AI system (heuristic_minus_ai.tsv, from
                    revision/response/_ai_comparison/heur_vs_ai.py in the manuscript repository)
Keys are optional: only the files whose keys are given are restaged.
"""
import json
import os
import shutil
import sys

import pandas as pd

data = sys.argv[1]
src = dict(a.split("=", 1) for a in sys.argv[2:])
D = os.path.join(data, "derived")


def wr(df, sub, name):
    os.makedirs(os.path.join(D, sub), exist_ok=True)
    df.to_csv(os.path.join(D, sub, name), sep="\t", index=False)
    print(sub, name, df.shape)


def cp(path, sub, name=None):
    os.makedirs(os.path.join(D, sub), exist_ok=True)
    dst = os.path.join(D, sub, name or os.path.basename(path)); shutil.copyfile(path, dst); print(sub, os.path.basename(dst))


if "unc_svs" in src and "unc_g1" in src:
    U = pd.read_csv(src["unc_svs"], sep="\t")
    G = pd.read_csv(src["unc_g1"], sep="\t"); G = G[G.column == "c"].set_index("qid")
    sv, hc = [], []
    for sc, st, q in (("all", "all_4364", "all_4364"), ("test", "test_test_heldout", "test_1094")):
        a = U[(U.scope == sc) & (U.metric == "AUROC")].iloc[0]; p = U[(U.scope == sc) & (U.metric == "AUPRC")].iloc[0]
        g = G.loc[q]; assert int(g.n_pairs) == int(a.n_pairs) and int(g.n_boot) == 2000
        sv.append(dict(stratum=st, metric="auroc", seed_sd=a.seed_sd, sample_sd_cluster=a.cluster_sd))
        hc.append(["AUROC p", sc, int(g.n_pairs), int(g.n_pos), g.auc, g.lo, g.hi, 2000])
        hc.append(["AUPRC p", sc, int(p.n_pairs), int(p.n_pos), p.ensemble, p.cluster_lo, p.cluster_hi, 2000])
    wr(pd.DataFrame(sv), "uncertainty", "seed_vs_sample_sd.tsv")
    wr(pd.DataFrame(hc, columns=["metric", "stratum", "n", "n_pos", "value", "lo", "hi", "B"]), "uncertainty", "headline_cosine_ci.tsv")

if "placement" in src:
    cp(src["placement"], "placement_confidence", "exclusion_sensitivity.tsv")

if "pop_dir" in src:
    for f in ("summary_popularity_strata.tsv", "strata_boundaries.tsv", "decision_popularity_strata.json"):
        cp(os.path.join(src["pop_dir"], f), "popularity_strata")

if "fda_prior" in src:
    for f in ("endpoints.tsv", "paired_differences.tsv", "flow.tsv"):
        cp(os.path.join(src["fda_prior"], f), "si_fda")

if "fda_dir" in src:
    cp(os.path.join(src["fda_dir"], "broad_only_fda_dev_results.tsv"), "si_fda")
    P = pd.read_csv(os.path.join(src["fda_dir"], "broad_only_fda_dev_pair_scores_ens5.tsv"), sep="\t")
    shipped = os.path.join(D, "si_fda", "broad_only_fda_dev_pair_scores_ens5.tsv")
    B = pd.read_csv(shipped, sep="\t"); R = pd.read_csv(src["fda_ref_raw"], sep="\t")
    assert len(B) == len(R) and (B.pair_class.values == R.pair_class.values).all()
    tok = {}
    for c in ("mondo_a", "mondo_b"):
        for raw, b in zip(R[c], B[c]):
            if raw != b:
                assert tok.get(raw, b) == b; tok[raw] = b
    assert set(P.mondo_a) | set(P.mondo_b) <= set(R.mondo_a) | set(R.mondo_b), "new pairs outside the reference"
    for c, l in (("mondo_a", "label_a"), ("mondo_b", "label_b")):
        m = P[c].isin(tok); P.loc[m, l] = None; P[c] = P[c].map(lambda x: tok.get(x, x))
    assert not P.astype(str).apply(lambda s: s.str.contains("holdout", case=False)).any().any()
    wr(P, "si_fda", "broad_only_fda_dev_pair_scores_ens5.tsv")

if "llm_sub" in src:
    # gene-anchored cosine rows of the shared-draw subgroup bootstrap; the schema of the shipped table is kept
    # (system "nd" = node-disjoint DisAtlas cosine, "p" = in-graph DisAtlas cosine)
    L = pd.read_csv(src["llm_sub"], sep="\t"); L = L[(L.readout == "cosine") & (L.B == 2000)].drop(columns="readout")
    L["system"] = L.system.replace({"nd_c": "nd", "c": "p"}); assert len(L) == 20 and (L.system == "nd").sum() == 4
    wr(L, "ixidn_robustness", "subgroups_shared_bootstrap_B2000.tsv")

if "baselines_dir" in src:
    # SI baselines table: AUROCs and paired differences (DisAtlas minus score) from the gene-anchored cosine audit,
    # B = 2,000 disease-cluster bootstrap. Only the score rows that the table prints are kept.
    b = src["baselines_dir"]
    def norm(f, block, keep):
        t = pd.read_csv(os.path.join(b, f), sep="\t")
        t = t.rename(columns={"d_ndc_minus_this": "d", "d_ref_minus_this": "d", "column": "score"})
        t["subset"] = t.qid.str.contains("test").map({True: "test", False: "all"}); t["block"] = block
        t = t[t.score.isin(keep)]
        return t[["block", "subset", "score", "auc", "lo", "hi", "d", "dlo", "dhi"]]
    IG = ["c", "tfidf_graph_svd_cosine__full_graph", "adjacency_svd_binary_cosine", "adjacency_svd_weighted_cosine",
          "adamic_adar__full_graph", "common_neighbours__full_graph", "weighted_adamic_adar__full_graph",
          "degree_product__full_graph", "allfeat_jaccard", "gene_jaccard", "hpo_jaccard"]
    HO = ["nd_c", "tfidf_svd_place55", "aa_place55", "waa_place55", "cn_place55", "tfidf_graph_svd_cosine",
          "adamic_adar", "common_neighbours", "allfeat_jaccard"]
    OE = ["common_neighbours", "adamic_adar", "tfidf_graph_svd_cosine"]
    T = pd.concat([norm("blockA_baselines_g1_B2000_geneanchored.tsv", "in_graph", IG),
                   norm("heldout_baselines_g1_B2000_geneanchored.tsv", "held_out", HO),
                   norm("own_edges_contrasts_g1_B2000_geneanchored.tsv", "held_out_own_edges", OE)])
    O = pd.read_csv(os.path.join(b, "OA_direct_measures.tsv"), sep="\t"); O = O[O.measure == "mondo_simgic"]
    O = pd.DataFrame({"block": "in_graph", "subset": O.scope.values, "score": "mondo_simgic", "auc": O.meas_auroc.values})
    T = pd.concat([T, O], ignore_index=True)
    assert len(T) == 2 * len(IG) + 2 * len(HO) + len(OE) + 2, len(T)
    wr(T, "si_baselines", "baselines_auroc_B2000.tsv")
if "heur_ai" in src:
    # TF-IDF graph SVD, Adamic-Adar and common neighbours (held-out diseases placed from 55 linked retained diseases)
    # minus each AI system, paired disease-cluster bootstrap (Poisson(1) per disease), B = 2,000, one set of draws
    # shared by every score; "(self)" rows give each score's own interval
    H = pd.read_csv(src["heur_ai"], sep="\t")
    # plus TF-IDF graph SVD with held-out diseases keeping their own annotation edges (same draws; SI text after the table)
    HS = ["TF-IDF graph SVD, placed", "Adamic-Adar, placed", "Common neighbours, placed", "TF-IDF graph SVD, own edges"]
    SY = ["Claude Chat", "Claude Sci", "Gemini Pro 3.1", "Phylo Biomni", "TxGNN"]
    H = H[(H.heuristic.isin(HS) & H.system.isin(SY)) | ((H.system == "(self)") & H.heuristic.isin(HS + SY))]
    assert len(H) == 2 * (len(HS) * len(SY) + len(HS) + len(SY)), len(H)
    wr(H, "si_heur_vs_ai", "heuristic_minus_ai_B2000.tsv")
