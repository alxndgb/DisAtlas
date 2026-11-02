"""Supplementary Table SI_degree_controls: degree-matched and within-disease AUROC of every system on the 4,357
scorable IxIDN-Matched pairs and on their 1,093-pair test split.

Degree-matched: one subset shared by all systems, with equal positives and negatives in each of ten quantile bins of
the summed log training-graph degree of the two diseases. Within-disease: mean per-disease AUROC over diseases with
both a positive and a negative partner. Intervals are disease-cluster bootstrap 95% CIs (2,000 draws, the same draws
for every system); differences are DisAtlas with benchmark diseases held out minus the system.
Third-party system scores are stored outputs; nothing is called here.

Input:  data/derived/degree_controls/degree_controls_table.tsv
Output: table_degree_controls.tsv and table_degree_controls.md next to this script
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T = pd.read_csv(ROOT / "data/derived/degree_controls/degree_controls_table.tsv", sep="\t")

SYSTEMS = [("DisAtlas", "DisAtlas, trained with every disease"), ("DisAtlas node-disjoint", "DisAtlas, benchmark diseases held out"),
           ("Claude Chat", "Claude Chat"), ("Claude Sci", "Claude Sci"), ("Gemini Pro 3.1", "Gemini Pro 3.1"),
           ("Phylo Biomni", "Phylo Biomni"), ("TxGNN (published checkpoint)", "TxGNN (published checkpoint)")]
SCOPES = [("all", "All pairs"), ("test", "Test split")]
NS = "\\ "  # pandoc non-breaking space, as in the manuscript source


def num(x):
    return f"{x:.3f}".replace("-", "−")


def ci(a, lo, hi):
    return f"{num(a)}{NS}({num(lo)}–{num(hi)})"


def delta(a, lo, hi):
    if pd.isna(a):
        return "—"
    s = ("+" if a >= 0 else "") + num(a)
    return f"{s}{NS}({num(lo)}{NS}to{NS}{num(hi)})"


rows, out = [], []
for sc, lab in SCOPES:
    S = T[T.case_set == sc].set_index("system")
    r0 = S.iloc[0]
    rows.append([lab, "Degree alone (training graph)", ci(r0.dm_degree_alone_auroc, r0.dm_degree_alone_ci_lo, r0.dm_degree_alone_ci_hi), "—", "—", "—"])
    for key, name in SYSTEMS:
        r = S.loc[key]
        rows.append([lab, name, ci(r.dm_auroc, r.dm_ci_lo, r.dm_ci_hi), delta(r.dm_delta_heldout_minus_system, r.dm_delta_ci_lo, r.dm_delta_ci_hi),
                     ci(r.within_auroc, r.within_ci_lo, r.within_ci_hi),
                     delta(r.within_delta_heldout_minus_system, r.within_delta_ci_lo, r.within_delta_ci_hi)])
cols = ["Pairs", "System", "Degree-matched AUROC (95% CI)", "Held out − system (95% CI)", "Within-disease AUROC (95% CI)", "Held out − system (95% CI)"]
pd.DataFrame(rows, columns=cols).to_csv(HERE / "table_degree_controls.tsv", sep="\t", index=False)
md = ["| " + " | ".join(cols) + " |", "|:---------|:--------------------------|----------------:|--------------:|---------------:|--------------:|"]
md += ["| " + " | ".join(r) + " |" for r in rows]
(HERE / "table_degree_controls.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
