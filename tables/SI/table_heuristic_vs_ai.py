"""Supplementary Table SI_heuristic_vs_ai: TF-IDF graph SVD cosine (benchmark diseases held out and placed from 55 linked
retained diseases) against the AI systems on the 4,357 node-disjoint IxIDN-Matched pairs and the 1,093-pair test split.

The AI-system scores are stored outputs (no model is called here).  AUROCs, their intervals and the paired differences
come from one disease-cluster bootstrap (Poisson(1) weight per disease, 2,000 draws, the same draws for every score);
unscored AI pairs are ranked last.

The input also holds Adamic-Adar and common neighbours (placed the same way) and TF-IDF graph SVD with held-out
diseases keeping their own annotation edges, all on the same bootstrap draws; the table shows placed TF-IDF graph SVD only, and the numbers file records, for every placed score, subset and system, the paired difference,
its interval and whether the interval includes zero (the SI text after the table states these).

Input:  data/derived/si_heur_vs_ai/heuristic_minus_ai_B2000.tsv
Output: table_heuristic_vs_ai.tsv, table_heuristic_vs_ai.md and table_heuristic_vs_ai_numbers.json next to this script
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T = pd.read_csv(ROOT / "data/derived/si_heur_vs_ai/heuristic_minus_ai_B2000.tsv", sep="\t")
HEUR = "TF-IDF graph SVD, placed"
HEURS = [HEUR, "Adamic-Adar, placed", "Common neighbours, placed", "TF-IDF graph SVD, own edges"]
SYSTEMS = ["Claude Chat", "Claude Sci", "Gemini Pro 3.1", "Phylo Biomni", "TxGNN"]
SUBSETS = [("all_4357", "All pairs", 4357), ("test_1093", "Test split", 1093)]


def ci(r):
    return f"{r.heuristic_auroc:.3f} ({r.ci_lo:.3f}–{r.ci_hi:.3f})"


def one(sub, h, system):
    r = T[(T.subset == sub) & (T.heuristic == h) & (T.system == system)]
    assert len(r) == 1, (sub, h, system); return r.iloc[0]


rows = []
for sub, lab, n in SUBSETS:
    for i, a in enumerate(SYSTEMS):
        d = one(sub, HEUR, a)
        rows.append([lab if i == 0 else "", f"{n:,}" if i == 0 else "", ci(one(sub, HEUR, "(self)")) if i == 0 else "",
                     a, ci(one(sub, a, "(self)")), f"{d['diff']:+.3f} ({d.ci_lo:.3f} to {d.ci_hi:.3f})".replace("-", "−")])
H = ["Pairs", "n", "TF-IDF graph SVD cosine (95% CI)", "System", "System AUROC (95% CI)", "Score − system (95% CI)"]
pd.DataFrame(rows, columns=H).to_csv(HERE / "table_heuristic_vs_ai.tsv", sep="\t", index=False)
md = ["| " + " | ".join(H) + " |", "|:--|--:|:--|:--|:--|:--|"] + ["| " + " | ".join(r) + " |" for r in rows]
(HERE / "table_heuristic_vs_ai.md").write_text("\n".join(md) + "\n")
num = {}
for h in HEURS:
    for sub, lab, n in SUBSETS:
        s = one(sub, h, "(self)")
        num[f"{h}|{sub}|auroc"] = dict(value=round(s.heuristic_auroc, 6), lo=round(s.ci_lo, 6), hi=round(s.ci_hi, 6))
        for a in SYSTEMS:
            d = one(sub, h, a)
            num[f"{h}|{sub}|{a}"] = dict(diff=round(d["diff"], 6), lo=round(d.ci_lo, 6), hi=round(d.ci_hi, 6),
                                         includes_zero=bool(d.ci_lo <= 0 <= d.ci_hi))
(HERE / "table_heuristic_vs_ai_numbers.json").write_text(json.dumps(num, indent=1, ensure_ascii=False) + "\n")
print(len(rows), len(num))
