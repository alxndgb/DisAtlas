"""Supplementary Table SI_llm_subgroups: DisAtlas (benchmark diseases held out) against language-model-based systems
by subgroup of the 4,357 scorable IxIDN-Matched pairs.

The language-model scores are stored outputs (no model is called here); the AUROCs and paired differences come from
one disease-cluster bootstrap (2,000 draws, the same draws for every system).

Input:  data/derived/ixidn_robustness/subgroups_shared_bootstrap_B2000.tsv
Output: table_llm_subgroups.tsv and table_llm_subgroups.md next to this script
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
B = pd.read_csv(ROOT / "data/derived/ixidn_robustness/subgroups_shared_bootstrap_B2000.tsv", sep="\t")

SCOPES = [("all", "All pairs"), ("test", "Test split"), ("non_oncology", "Non-oncology pairs"), ("oncology", "Oncology pairs")]
SYSTEMS = ["Claude Chat", "Claude Sci", "Gemini Pro 3.1"]


def ci(a, lo, hi):
    return f"{a:.3f} ({lo:.3f}–{hi:.3f})"


def d(a, lo, hi):
    f = lambda x: f"{x:+.3f}".replace("-", "−") if x < 0 else f"{x:+.3f}"
    g = lambda x: f"{x:.3f}".replace("-", "−")
    return f"{f(a)} ({g(lo)} to {g(hi)})"


rows = []
for sc, lab in SCOPES:
    nd = B[(B.scope == sc) & (B.system == "nd")].iloc[0]
    for s in SYSTEMS:
        r = B[(B.scope == sc) & (B.system == s)].iloc[0]
        rows.append({"Pairs": lab, "n": f"{int(nd.n_pairs):,}", "DisAtlas, held out (95% CI)": ci(nd.auroc, nd.lo, nd.hi),
                     "System": s, "System AUROC (95% CI)": ci(r.auroc, r.lo, r.hi),
                     "DisAtlas − system (95% CI)": d(r.dis_minus_sys, r.dms_lo, r.dms_hi)})
T = pd.DataFrame(rows)
T.to_csv(HERE / "table_llm_subgroups.tsv", sep="\t", index=False)
md = ["| " + " | ".join(T.columns) + " |", "|" + "|".join([":--", "--:", ":--", ":--", ":--", ":--"]) + "|"]
md += ["| " + " | ".join(map(str, r)) + " |" for r in T.itertuples(index=False)]
(HERE / "table_llm_subgroups.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
