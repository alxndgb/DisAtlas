"""Maintainer-only: stage the rare/common inputs of the gene-anchored model (Supplementary Figures SI_rare_common_stratified
and SI_rare_common_neighbour_coverage, Supplementary Table SI_rare_common_strata) into data/derived/si_rare_common.
    python scripts/maintainer/stage_si_rare_common.py STRATA=<dir> COV=<dir>
STRATA: ccplus_strata.tsv, r4mc3.tsv, four.tsv and the four <group>.tsv per-anchor tables (embedding cosine);
COV: coverage_by_k_summary.tsv, enrichment_other_class_frac10.tsv, per_disease_first_other_class.tsv.gz, meta.json."""
import sys, os, json, shutil
A = dict(a.split("=", 1) for a in sys.argv[1:])
DST = os.path.join(os.path.dirname(__file__), "..", "..", "data", "derived", "si_rare_common")
G = ["rare_anchor_rank_rare_partners", "rare_anchor_rank_common_partners", "common_anchor_rank_rare_partners",
     "common_anchor_rank_common_partners"]
for f in ["ccplus_strata", "r4mc3", "four"] + G:
    shutil.copyfile(os.path.realpath(os.path.join(A["STRATA"], f + ".tsv")), os.path.join(DST, f + ".tsv"))
for f in ["coverage_by_k_summary.tsv", "enrichment_other_class_frac10.tsv", "per_disease_first_other_class.tsv.gz"]:
    shutil.copyfile(os.path.join(A["COV"], f), os.path.join(DST, f))
M = json.load(open(os.path.join(A["COV"], "meta.json")))
keep = ["n", "n_rare", "n_common", "nperm", "median_first_common_rank_for_rare", "median_first_rare_rank_for_common",
        "n_rare_no_common_in_50", "check_k10_other_frac"]
M = {k: M[k] for k in keep}; M["score"] = "gene-anchored model, five-seed ensemble embedding cosine"
json.dump(M, open(os.path.join(DST, "rc_coverage_meta.json"), "w"), indent=1)
print("staged", len(G) + 6, "files")
