"""Maintainer-only: stage data/derived/si_baselines/baseline_pair_scores.tsv.gz, the per-pair training-free scores
behind Supplementary Table S5 (si_baselines) for the 4,364 IxIDN benchmark pairs. Read-only sources:
  * annotation Jaccards (all-feature, gene, HPO; from restricted-source annotations, so not recomputable from the release):
    /scratch/kwfh644/reframe_v1/si/promote/r4mc1_baselines_pair_scores_full_setting_cosine.tsv.gz
    (the in-graph graph scores are recomputed by si_baselines from data/raw/training_graph and are not staged),
  * held-out block: /scratch/kwfh644/learning_smoke/test6_features_only/out_6b/pair_scores_6b.tsv.gz (same-placement
    "_place55" and own-edges scores, columns prefixed ho_),
  * MONDO simGIC row: the authors' OA_direct_measures.tsv (AUROC only; staged by stage_figD.py, not a per-pair score).
No DisAtlas or AI-system score is staged: DisAtlas c and nd_c are read from data/derived/fig6_comparison.
    python scripts/maintainer/stage_baseline_pair_scores.py OUT.tsv.gz
"""
import sys
import pandas as pd
R0 = "/projects/rd/us/AI/kwfh644/dis2vec_2"
S4 = pd.read_csv(f"{R0}/revision/frozen_mechdist_eval/data/ixidn4364_scores.tsv", sep="\t", usecols=["blind_id", "a", "b", "y"])
F = pd.read_csv("/scratch/kwfh644/reframe_v1/si/promote/r4mc1_baselines_pair_scores_full_setting_cosine.tsv.gz", sep="\t")
F = F[["pair_id", "y", "split", "allfeat_jaccard", "gene_jaccard", "hpo_jaccard"]]  # graph scores: computed in si_baselines from the released training graph
H = pd.read_csv("/scratch/kwfh644/learning_smoke/test6_features_only/out_6b/pair_scores_6b.tsv.gz", sep="\t")
HO = ["tfidf_svd_place55", "aa_place55", "waa_place55", "cn_place55", "tfidf_graph_svd_cosine", "adamic_adar", "common_neighbours", "allfeat_jaccard"]
H = H[["pair_id"] + HO].rename(columns={c: f"ho_{c}" for c in HO})
T = S4.merge(F.rename(columns={"y": "y_f"}), left_on="blind_id", right_on="pair_id", how="left", validate="1:1").drop(columns="pair_id")
T = T.merge(H, left_on="blind_id", right_on="pair_id", how="left", validate="1:1").drop(columns="pair_id")
assert len(T) == 4364 and (T.y == T.y_f).all(); T = T.drop(columns="y_f")
T.to_csv(sys.argv[1], sep="\t", index=False, compression="gzip"); print(T.shape, list(T.columns))
