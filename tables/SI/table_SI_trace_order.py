"""Supplementary Table S18 (SI_trace_order) and its companion values table (SI_trace_order_values, not printed in the final
SI; kept as a README output): display order of the Figure 4 evidence traces.

Table S18 lists, per pair and list, the intermediates that ordering by walk-mass share x specificity weight removes from and adds
to the displayed list, compared with ordering by walk-mass share alone (no score changes). The values table lists every displayed
intermediate in display order with its walk-mass share, specificity weight ln(N / degree) (N = 13,021 diseases) and, where
Figure 4 shows a shortened name, the display label.
Input: data/derived/fig4_evidence_trace/trace_displayed_items.tsv (the displayed intermediates under either ordering).
Outputs next to this script: table_SI_trace_order.{md,tsv}, table_SI_trace_order_values.{md,tsv}
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
D = ROOT / "data/derived/fig4_evidence_trace"
T = pd.read_csv(D / "trace_displayed_items.tsv", sep="\t", keep_default_na=False, na_values=[""])
assert json.load(open(D / "fig4_values.json"))["idf_n_diseases"] == 13021
CASE = ["aHUS-PNH", "NMOSD-PNH", "Crohn-UC", "RA-polymyositis"]
PN = {"aHUS-PNH": "aHUS–PNH", "NMOSD-PNH": "NMOSD–PNH", "Crohn-UC": "Crohn–UC", "RA-polymyositis": "RA–polymyositis"}
LN = {"module": "evidence items", "gene": "genes", "neighbour": "shared neighbor diseases"}

# Table S18: lists whose displayed set changes (list order as in the trace summary: gene, module, neighbour)
c22 = ["Pair (panel)", "List", "Shown by walk-mass share only", "Shown with specificity weighting"]
r22 = []
for p in CASE:
    for lst in ("gene", "module", "neighbour"):
        g = T[(T.pair == p) & (T.list == lst)].sort_values("rank_share_only", kind="stable")
        drop = g[g.shown_share_only & ~g.shown_specificity].label.tolist(); add = g[~g.shown_share_only & g.shown_specificity].label.tolist()
        assert len(drop) == len(add)
        if drop:
            r22.append([f"{PN[p]} ({g.panel.iloc[0]})", LN[lst], ", ".join(drop), ", ".join(add)])
pd.DataFrame(r22, columns=c22).to_csv(HERE / "table_SI_trace_order.tsv", sep="\t", index=False)
md = ["| " + " | ".join(c22) + " |", "|:-----------|:-----------|:-----------------------|:-----------------------|"] + ["| " + " | ".join(r) + " |" for r in r22]
(HERE / "table_SI_trace_order.md").write_text("\n".join(md) + "\n", encoding="utf-8")

# values table: displayed intermediates in display order
same = lambda a, b: a == b or (a[:1].upper() + a[1:]) == b
c23 = ["Pair (panel)", "List", "Rank", "Displayed intermediate", "Walk-mass share", "Specificity weight (ln N/degree)", "Display label in Figure 4"]
r23 = []
for p in CASE:
    for lst in ("module", "gene", "neighbour"):
        g = T[(T.pair == p) & (T.list == lst) & T.shown_specificity].sort_values("rank_specificity")
        for r in g.itertuples():
            dl = "" if same(r.label, r.display_label) else r.display_label
            r23.append([f"{PN[p]} ({r.panel})", LN[lst], str(int(r.rank_specificity)), r.label, f"{r.frac_of_M2:.3f}", f"{r.idf:.2f}", dl])
pd.DataFrame(r23, columns=c23).to_csv(HERE / "table_SI_trace_order_values.tsv", sep="\t", index=False)
md = ["| " + " | ".join(c23) + " |", "|:-----------|:-----------|-----:|:-----------------------|-----:|-----:|:-----------|"]
md += ["| " + " | ".join(r[:-1]) + " | " + (f"{r[-1]} |" if r[-1] else "|") for r in r23]
(HERE / "table_SI_trace_order_values.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print(len(r22), len(r23), sum(bool(r[-1]) for r in r23))
