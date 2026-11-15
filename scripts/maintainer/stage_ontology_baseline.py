"""Maintainer-only: stage the per-pair ontology similarities used by disatlas.ontology_baseline into
data/derived/ontology_baseline/ixidn4364_hpo_similarity.tsv (HPO Resnik-BMA, Lin-BMA, SimGIC, direct MechDist HPO
distance and the MONDO SimGIC label-only control; empty = undefined because a disease lacks annotations).
    python scripts/maintainer/stage_ontology_baseline.py SEMSIM=<ixidn4364_pairs_scored.csv>"""
import sys, os
import pandas as pd
A = dict(a.split("=", 1) for a in sys.argv[1:])
DST = os.path.join(os.path.dirname(__file__), "..", "..", "data", "derived", "ontology_baseline", "ixidn4364_hpo_similarity.tsv")
COLS = ["blind_id", "hpo_resnik_bma", "hpo_lin_bma", "hpo_simgic", "mechdist_hpo", "mondo_simgic"]
S = pd.read_csv(A["SEMSIM"], usecols=COLS)[COLS]
assert len(S) == 4364 and S.blind_id.is_unique
S.to_csv(DST, sep="\t", index=False, float_format="%.10g"); print("staged", DST, S.shape)
