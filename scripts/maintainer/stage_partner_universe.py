"""Maintainer-only: write data/derived/si_rare_common/partner_universe.tsv, the candidate universe of the
rare-common documented-partner recovery (disease ids in the order of the score matrix used for the partner percentiles,
with their rare/common class), so the TxGNN column (scripts/maintainer/txgnn_partner_percentile.py) ranks against
exactly the candidates of partner_percentile_per_pair.tsv. Not needed to reproduce any result (the table ships).

Reads the internal score-matrix id list (JSON list, or {"ids": [...]}, path in DISATLAS_SCORE_IDS) and checks it against the released data: the ids are
the released all-pairs score ids in the same order (two diseases not in the internal matrix), and the classes equal
data/raw/benchmarks/disease_labels.tsv.

usage: DISATLAS_SCORE_IDS=<ids.json> python scripts/maintainer/stage_partner_universe.py OUT.tsv
"""
import json, os, sys
import numpy as np, pandas as pd

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
IDS = os.environ["DISATLAS_SCORE_IDS"]
uid = json.load(open(IDS)); uid = uid["ids"] if isinstance(uid, dict) else uid
lab = pd.read_csv(os.path.join(REPO, "data", "raw", "benchmarks", "disease_labels.tsv"), sep="\t").set_index("disease_id")
rel = [str(x) for x in np.load(os.path.join(REPO, "data", "zenodo", "scores", "disatlas_all_pairs_ensemble.npz"))["ids"]]
assert [d for d in rel if d in set(uid)] == list(uid), "order differs from the released score ids"
U = pd.DataFrame({"disease_id": uid, "disease_class": lab.disease_class.reindex(uid).values})
assert U.disease_class.notna().all()
U.to_csv(sys.argv[1], sep="\t", index=False)
print(len(U), U.disease_class.value_counts().to_dict(), "not in matrix:", sorted(set(rel) - set(uid)))
