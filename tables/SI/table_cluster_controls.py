"""Supplementary Table SI_cluster_controls: positive and negative control disease pairs at K = 55.

Pair names and roles are read from data/derived/si_cluster/SI_controls_table.tsv; same-cluster flag, DisAtlas pair link
score (mean of five seeds), percentile among all scored pairs and best partner rank are read from
figures/SI/SI_cluster_rarecommon_numbers.json (written by figures/SI/si_cluster_rarecommon_numbers.py).
Outputs next to this script: table_cluster_controls.{md,tsv}
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T = pd.read_csv(ROOT / "data/derived/si_cluster/SI_controls_table.tsv", sep="\t").set_index("tag")
N = json.load(open(ROOT / "figures/SI/SI_cluster_rarecommon_numbers.json"))
cols = ["Control pair", "Role", "Same cluster at K = 55", "DisAtlas pair link score", "Percentile among all pairs",
        "Best partner rank"]
rows = [[T.loc[c["tag"], "pair"], c["role"], c["same"], f"{c['score']:.3f}", f"{c['pct']:.1f}", f"{c['best_rank']:,}"]
        for c in N["controls"]]
pd.DataFrame(rows, columns=cols).to_csv(HERE / "table_cluster_controls.tsv", sep="\t", index=False)
md = ["| " + " | ".join(cols) + " |", "|:--|:--|:--|--:|--:|--:|"] + ["| " + " | ".join(r) + " |" for r in rows]
(HERE / "table_cluster_controls.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print(len(rows), N["n_pairs_all"], N["median_link_score_all_pairs"])
