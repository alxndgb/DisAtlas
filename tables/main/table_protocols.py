"""Main text, Methods "Comparator systems" -- Table tbl:protocols.

Descriptive cells are fixed protocol text; every count in the table (runs per LLM system,
Phylo seeds, unscored / unmapped pairs) is computed from the stored per-run comparator
scores on the 4,364 IxIDN-Matched pairs (tables/inputs/comparator_run_scores.tsv; LLM
scores are the saved outputs of the runs used in the paper, no model is called here).
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tablib import INPUTS, write_and_check  # noqa: E402

runs = pd.read_csv(INPUTS / "comparator_run_scores.tsv", sep="\t")
assert len(runs) == 4364
systems = {}
for col in runs.columns:
    if "|" in col:
        sysname, run = col.split("|", 1)
        systems.setdefault(sysname, []).append(run)

n_runs = {k: len(v) for k, v in systems.items()}
llm = {"Claude-Chat Opus 5.5": "Claude Chat (Opus 5.5, {} runs)",
       "Gemini Pro 3.1": "Gemini Pro 3.1 ({} runs)"}
llm_examples = ", ".join(t.format(n_runs[k]) for k, t in llm.items())
agentic_examples = f"Claude Sci (Claude Science, Opus 5.5, {n_runs['Claude-Sci Opus 5.5']} runs); Phylo Biomni"

seeds = sorted(int(re.search(r"seed(\d+)", r).group(1)) for r in systems["Phylo"])
assert seeds == list(range(seeds[0], seeds[-1] + 1)) and len(seeds) == 5
txgnn_unmapped = int(runs["TxGNN|DistMult"].isna().sum())

rows = [
    ["1", "DisAtlas-UPNA-DDAL", "DisAtlas (unseen diseases)",
     "Trained on MechDist with every benchmark disease withheld from embedding and decoder "
     "training; evaluated on the blinded IxIDN-Matched pairs"],
    ["2", "Large language model", llm_examples,
     "Structured prompt with the blinded pair file; mean score over runs"],
    ["3", "Agentic AI system", agentic_examples, "Mean score over runs or seeds"],
    ["4", "Knowledge-graph foundation model (published checkpoint)", "TxGNN [@TxGNN]",
     "Published checkpoint, DistMult decoder; not retrained on MechDist or on IxIDN-Matched; "
     f"{txgnn_unmapped} unmapped pairs ranked last"],
]
df = pd.DataFrame(rows, columns=["Protocol", "System type", "Examples", "Evaluation approach"])

EXPECTED = [
    ["1", "DisAtlas-UPNA-DDAL", "DisAtlas (unseen diseases)",
     "Trained on MechDist with every benchmark disease withheld from embedding and decoder "
     "training; evaluated on the blinded IxIDN-Matched pairs"],
    ["2", "Large language model", "Claude Chat (Opus 5.5, 5 runs), Gemini Pro 3.1 (3 runs)",
     "Structured prompt with the blinded pair file; mean score over runs"],
    ["3", "Agentic AI system", "Claude Sci (Claude Science, Opus 5.5, 5 runs); Phylo Biomni",
     "Mean score over runs or seeds"],
    ["4", "Knowledge-graph foundation model (published checkpoint)", "TxGNN [@TxGNN]",
     "Published checkpoint, DistMult decoder; not retrained on MechDist or on IxIDN-Matched; "
     "213 unmapped pairs ranked last"],
]

if __name__ == "__main__":
    write_and_check(df, __file__, "tbl_protocols", EXPECTED,
                    "Summary of cross-evidence generalization protocols applied across system categories.")
