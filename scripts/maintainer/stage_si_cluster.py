"""Maintainer-only: stage the SI cluster-validation inputs of the gene-anchored model (Figures SI_cluster_validation and
SI_cluster_map_vs_model, Tables SI_cluster_controls and SI_cluster_themes) from the analysis run directories.
Not needed to reproduce the figures; documents where the shipped tables in data/derived come from.
    python scripts/maintainer/stage_si_cluster.py THEMES=<dir> KSEL=<dir> SIFU=<dir> CTRL=<dir>
THEMES: per-cluster annotation tables and themes_k55.csv of the theme audit on the K = 55 cosine k-medoids labels;
KSEL: k_selection_boot.tsv.gz; SIFU: within_between_K55_13021.tsv.gz and sifu_cluster_13021_values.json;
CTRL: SI_controls_table.tsv (the 17 control pairs with names, roles, MONDO ids, clusters and cosine scores)."""
import sys, os, json, shutil
import pandas as pd
A = dict(a.split("=", 1) for a in sys.argv[1:])
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DV = os.path.join(REPO, "data", "derived")
# --- cluster themes: same columns as the shipped parquet files
TH = os.path.join(DV, "si_cluster_themes")
for name, src in [("bodysystem_by_cluster_long", "bodysystem_by_cluster_long.tsv"),
                  ("mechanism_terms_by_cluster", "mechanism_terms_by_cluster.tsv"),
                  ("go_bp_top5_by_cluster", "go_bp_top5_by_cluster.tsv"),
                  ("mondo_terms_top5_by_cluster", "mondo_terms_top5_by_cluster.tsv"),
                  ("cluster_annotation_k55", "cluster_annotation_k55_raw.tsv")]:
    cols = pd.read_parquet(os.path.join(TH, name + ".parquet")).columns.tolist()
    T = pd.read_csv(os.path.join(A["THEMES"], src), sep="\t")
    miss = [c for c in cols if c not in T.columns]; assert not miss, (name, miss)
    T[cols].to_parquet(os.path.join(TH, name + ".parquet"), index=False)
    print(name, T[cols].shape, "dropped columns", len(T.columns) - len(cols))
R = pd.read_csv(os.path.join(A["THEMES"], "themes_k55.csv"))[["cluster", "theme"]].rename(columns={"theme": "theme_class"})
assert len(R) == 55; R.to_csv(os.path.join(TH, "theme_reference.tsv"), sep="\t", index=False)
# --- cluster validation inputs
SC = os.path.join(DV, "si_cluster")
shutil.copyfile(os.path.join(A["KSEL"], "k_selection_boot.tsv.gz"), os.path.join(SC, "k_selection_boot.tsv.gz"))
shutil.copyfile(os.path.join(A["SIFU"], "within_between_K55_13021.tsv.gz"), os.path.join(SC, "within_between_K55_13021.tsv.gz"))
V = json.load(open(os.path.join(A["SIFU"], "sifu_cluster_13021_values.json")))
V["layout"]["layout_file"] = "data/zenodo/model_outputs/clusters/disease_clusters_k55.parquet (isomap_x, isomap_y)"
V["bench_file"] = "IxIDN-Matched benchmark, 4,364 pairs"
V.pop("secs", None)
json.dump(V, open(os.path.join(SC, "sifu_cluster_13021_values.json"), "w"), indent=1)
# --- control pairs (cosine score, K = 55 cosine clusters); Table SI_cluster_controls recomputes scores, percentiles and
# ranks from the model outputs and uses only pair, role, tag and the ids from this file
C = pd.read_csv(os.path.join(A["CTRL"], "SI_controls_table.tsv"), sep="\t")
assert len(C) == 17; C.to_csv(os.path.join(SC, "SI_controls_table.tsv"), sep="\t", index=False)
print("staged", len(C), "controls")
