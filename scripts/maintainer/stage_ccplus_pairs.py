"""Maintainer-only: stage the IxIDN-Matched-CC+ pair table (Supplementary Table SI_rare_common_strata) into
data/derived/si_rare_common/ixidn_ccplus_pairs.tsv.
    python scripts/maintainer/stage_ccplus_pairs.py CC=<authors' headline_swap/out/ccplus dir>
Keeps the pair list, labels, pool, source (IxIDN-Matched or CC+ added), rarity stratum and the five matched-shortcut
features of the shortcut stack. mdeg (log1p training-graph degree of a plus that of b) is staged as computed by the
authors on the tune-split training graph, which is not part of the release. No model score is staged: the rebuild
reads the DisAtlas cosine from the released cosine matrix."""
import sys, os
import pandas as pd
A = dict(a.split("=", 1) for a in sys.argv[1:])
DST = os.path.join(os.path.dirname(__file__), "..", "..", "data", "derived", "si_rare_common", "ixidn_ccplus_pairs.tsv")
X = pd.read_csv(os.path.join(A["CC"], "ixidn_ccplus_pairs.tsv"), sep="\t")
S = pd.read_csv(os.path.join(A["CC"], "ixidn_ccplus_scores.tsv"), sep="\t")
assert list(X.pair_uid) == list(S.pair_uid) and len(X) == 4673
X["mdeg"] = S.mdeg.values
C = ["pair_uid", "a", "b", "label", "pool", "source", "stratum", "dsum", "mdeg", "mondo_simgic", "name_tfidf_cos", "genecov"]
X[C].to_csv(DST, sep="\t", index=False, float_format="%.10g")
print("staged", len(X), "pairs")
