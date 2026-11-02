"""Supplementary Table SI_fig3a_trials: the trials behind the disease-drug lines of Figure 3a.

One row per trial linking one of the four Figure 3a diseases to one of the four drugs, as recorded in the trial data
used to build IxIDN. Trial is the public registry identifier (– when the trial has none). Phase is the highest phase in
the trial data; † marks trials whose ClinicalTrials.gov record (ctgov_phase, a stored snapshot) gives a different phase.
Start is the trial start year (– when not recorded).

Input:  data/derived/fig3a_trials/fig3a_trials.csv
Output: table_fig3a_trials.tsv and table_fig3a_trials.md next to this script
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T = pd.read_csv(ROOT / "data/derived/fig3a_trials/fig3a_trials.csv", dtype=str, keep_default_na=False)


def dash(x):
    return x if x.strip() else "–"


def phase(r):
    p = "other" if r.phase == "OTHER" else r.phase
    return p + ("†" if r.ctgov_phase and r.ctgov_phase != r.phase else "")


rows = [[r.disease, r.drug, dash(r.public_id), phase(r), dash(r.start_year), r.n_drugs, r.title.replace("|", "\\|")]
        for r in T.itertuples(index=False)]
cols = ["Disease", "Drug", "Trial", "Phase", "Start", "Drugs in trial", "Title"]
out = pd.DataFrame(rows, columns=cols)
out.to_csv(HERE / "table_fig3a_trials.tsv", sep="\t", index=False)
md = ["| " + " | ".join(cols) + " |", "|:--|:--|:--|:--|:--|:--|:----------|"] + ["| " + " | ".join(r) + " |" for r in rows]
(HERE / "table_fig3a_trials.md").write_text("\n".join(md) + "\n")
summary = {"trials": len(out), "single_drug": int((T.n_drugs == "1").sum()), "with_public_id": int((out.Trial != "–").sum()),
           "phase_mismatch": int(out.Phase.str.endswith("†").sum())}
print(summary)
