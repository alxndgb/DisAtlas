"""Maintainer-only: stage the inputs of the placement-confidence exclusion analysis (Supplementary Figure
SI_placement_confidence) so that data/derived/placement_confidence/exclusion_sensitivity.tsv can be rebuilt from the
released cosine matrix (disatlas.rebuild.placement_exclusion).

    python scripts/maintainer/stage_placement_confidence.py DQ=<disease_dq_summary.tsv.gz> PAIRS=<frozen_mechdist_eval/data> TEST_DM=<frozen_mechdist_eval/data_v2.1/test_dm_scores.tsv>

Writes (no model score is staged; the scores are read from the released cosine matrix):
  data/derived/placement_confidence/disease_placement_flags.tsv.gz  per-disease evidence flags and issue classes
      (the per-disease table the SI text says is released with the map). The stored dq_tier column is not staged:
      it is not the critical-tier rule of the analysis (>= 3 issue classes among scored diseases), which is derived
      in the rebuild module.
  data/derived/placement_confidence/test_dm_pairs.tsv.gz   MechDist link-prediction test split (a, b, y; 13,448 pairs)
  data/derived/placement_confidence/ixidn910_pairs.tsv     IxIDN tuning split (a, b, y; 906 pairs)
"""
import sys, os
import pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src = dict(a.split("=", 1) for a in sys.argv[1:]); OUT = os.path.join(ROOT, "data/derived/placement_confidence")
Q = pd.read_csv(src["DQ"], sep="\t")
keep = (["disease_id", "label", "disease_class", "in_disatlas_scored", "n_issue_classes", "issue_classes", "flags"]
        + [c for c in Q.columns if c[:2] in ("C1", "C2", "C3", "C4", "C5", "C6")] + [c for c in Q.columns if c.startswith("flag_")])
assert len(Q) == 13021 and Q.disease_id.is_unique
Q[keep].to_csv(os.path.join(OUT, "disease_placement_flags.tsv.gz"), sep="\t", index=False)
# test_dm: the 13,448-pair split of the frozen v2.1 evaluation (data_v2.1); frozen_mechdist_eval/data holds a later
# re-draw of the same split size that the shipped exclusion table does not use.
for p, o in [(src["TEST_DM"], "test_dm_pairs.tsv.gz"), (os.path.join(src["PAIRS"], "ixidn910_scores.tsv"), "ixidn910_pairs.tsv")]:
    d = pd.read_csv(p, sep="\t")[["a", "b", "y"]]
    d.to_csv(os.path.join(OUT, o), sep="\t", index=False); print(o, len(d))
print("staged", len(Q), "diseases")
